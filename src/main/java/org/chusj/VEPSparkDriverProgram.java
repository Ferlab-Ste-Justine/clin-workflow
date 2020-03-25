package org.chusj;

import com.whalin.MemCached.MemCachedClient;
import com.whalin.MemCached.SockIOPool;
import org.apache.http.HttpHost;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.elasticsearch.action.DocWriteRequest;
import org.elasticsearch.action.bulk.BulkRequest;
import org.elasticsearch.action.bulk.BulkResponse;
import org.elasticsearch.action.index.IndexRequest;
import org.elasticsearch.action.update.UpdateRequest;

import org.elasticsearch.client.RequestOptions;
import org.elasticsearch.client.RestClient;
import org.elasticsearch.client.RestHighLevelClient;
import org.elasticsearch.common.xcontent.XContentType;
import org.elasticsearch.script.Script;
import org.elasticsearch.script.ScriptType;
import org.json.JSONArray;
import org.json.JSONObject;

import javax.xml.bind.DatatypeConverter;
import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.util.*;

import static org.chusj.PatientHelper.loadPedigree;
import static org.chusj.VepHelper.*;

public class VEPSparkDriverProgram {

    public static RestHighLevelClient client;
    private static SockIOPool pool;
    private static MemCachedClient mcc;
    private static String MUTATION_INDEX_NAME = "mutations";
    private static String GENE_INDEX_NAME = "genes";

    public static int TOTAL_COUNT = 0;
    private static List<GeneVariants> geneVariants = new ArrayList<>();


    public static void main(String[] args) throws Exception {

        if (args.length < 9 ) {
            throw new Exception("Missing params; need extractFile ETLPropertiesFile ES_UPSERT (true|false) " +
                    "sparkMaster local[nbWorkers] 8g nbPartitions bulkSize pedFile");
        }

        String extractFile = args[0];
        String pedigreePropsFile = args[1];

        boolean esUpsert = Boolean.parseBoolean(args[2]);
        String sparkMaster = args[3];
        String localThread = args[4];
        String memory = args[5];
        int nbPartitions = Integer.parseInt(args[6]);
        int bulkOpsQty = Integer.parseInt(args[7]);
        String pedFile = args[8];
        //spliting gene into a join relation is currently disabled as it's not helping in our use cases
        // The feature is kept for now as the simplification of the ES indexing was done at same time
        boolean splitGene = false; //  Boolean.parseBoolean(args[8]);
//        if (splitGene) {
//            INDEX_NAME = "mutations_genes";
//            System.err.println("$$$$$$$\n$$$$$$$$\nswitching to mutation_gene\n\n$$$$$$$$$");
//        }

        JavaSparkContext sc;
        try (RestHighLevelClient clientTry = new RestHighLevelClient(
                RestClient.builder(
                        new HttpHost("localhost", 9200, "http")))) {


            client = clientTry;
            PatientHelper.client = clientTry;

            Properties pedigreeProps = VepHelper.getPropertiesFromFile(pedigreePropsFile);
            //pedigreeProps.forEach( (k,v) -> System.out.println(k+"="+v) );

            List<Pedigree> pedigrees = loadPedigree(pedFile);
//            System.out.println("pedigrees");
            pedigrees.forEach(System.out::println);
            List<String> specimenList = getSpecimenList(extractFile);
            Map<String, Patient> patientMap = PatientHelper.preparePedigree(specimenList);
//            System.out.println("patientMap");
//            patientMap.forEach((k, v) -> System.out.println(k + "\n\t" + v));
            Map<String,Family> familyMap = getFamilyMap(specimenList, pedigrees);
//            System.out.println("familyMap");
//            familyMap.forEach((k, v) -> System.out.println(k + "\n\t" + v));

            /* Define Spark Configuration */
            SparkConf conf = new SparkConf().setAppName("ExtractTLoad").setMaster(sparkMaster)
                    .setMaster(localThread).set("spark.executor.memory",memory)
                    .set("spark.driver.allowMultipleContexts", "true");

            /* Create Spark Context with configuration */
            sc = new JavaSparkContext(conf);

            /* Create a Resilient Distributed Dataset for a long file
             * Each line in long file become a record in RDD
             * */
            JavaRDD<String> lines = sc.textFile(extractFile, nbPartitions);

            lines.foreachPartition(partitionOfRecords -> {
                List<JSONObject> mutationPayloadList = new ArrayList<>();
                List<JSONObject> genePayloadList = new ArrayList<>();
                while (partitionOfRecords.hasNext()) {
                    Variant variant = VepHelper.processVcfDataLine(partitionOfRecords.next(), pedigreeProps,
                            patientMap, pedigrees, specimenList, familyMap);

                    if (variant != null) {
                        for (Gene gene: variant.getGenes()) {
                            //GeneVariants geneVariants = new GeneVariants();
                            //geneVariants.setGene(gene);
                            //geneVariants.setVariant(variant);
                            //VEPSparkDriverProgram.geneVariants.add(geneVariants);
                            JSONObject genePayload = new JSONObject();
                            String mutationId = variant.getMutationId();
                            JSONArray variants = new JSONArray();
                            variants.put(mutationId);
                            String id = gene.getEnsemblId();
                            if (id != null && !id.isEmpty()) {
                                genePayload.put("id", id);
                                genePayload.put("geneSymbol", gene.getGeneSymbol());
                                genePayload.put("biotype", gene.getBiotype());
                                genePayload.put("ensemblId", gene.getEnsemblId());
                                genePayload.put("mutationId", mutationId);
                                genePayload.put("variants", variants);
                                genePayload.put("donors", new JSONArray());
                                genePayload.put("frequencies", new JSONArray());
                                JSONArray alias = new JSONArray();
                                // Add newSymbol and newAlias to track/update the gene index
                                if (gene.getNewAlias() != null) {
                                    alias.put(gene.getNewAlias());
                                }
                                genePayload.put("alias", alias);
                                genePayload.put("newSymbol", gene.getNewSymbol());
                                genePayload.put("newAlias", gene.getNewAlias());

                                genePayloadList.add(genePayload);
                            } else {
                                System.out.println("empty ensemblid for "+gene.toString());
                            }

                        }

                        JSONObject payload = new JSONObject(variant.getJsonObjInString());
                        mutationPayloadList.add(payload);

                    }
                    if (mutationPayloadList.size() >= bulkOpsQty) {
                        //System.out.println("Bulk Items in partition-" + mutationPayloadList.size());
                        bulkStoreJsonObj(mutationPayloadList, esUpsert, false, false);
                        mutationPayloadList = new ArrayList<>();
                    }
                    if (genePayloadList.size() >= bulkOpsQty) {
//                        System.out.println("Bulk gene Items in partition-" + genePayloadList.size());
                        bulkStoreJsonObj(genePayloadList, esUpsert, false, true);
                        genePayloadList = new ArrayList<>();
                    }
                }
//                System.err.println("Bulk Items left in partition-" + mutationPayloadList.size());
                bulkStoreJsonObj(mutationPayloadList, esUpsert,false, false);
                bulkStoreJsonObj(genePayloadList, esUpsert, false, true);


            });
            // time to split by genes
            //System.out.println("variant count = " + variants.size());
            System.out.println("//////////////");
            System.out.println("//////////////");
            System.out.println("//////////////");
            System.out.println("//////////////");


//            sc = new JavaSparkContext(conf);
//
//            JavaRDD<GeneVariants> blob = sc.parallelize(geneVariants, nbPartitions);
//            //System.out.println("GeneVariants count = " + blob.count());
//            blob.foreachPartition(partitionOfRecords -> {
//                List<JSONObject> jsonObjectList = new ArrayList<>();
//                while (partitionOfRecords.hasNext()) {
//                    JSONObject payload = new JSONObject();
//                    GeneVariants geneVariants = partitionOfRecords.next();
//                    String mutationId = geneVariants.getVariant().getMutationId();
//                    JSONArray variants = new JSONArray();
//                    variants.put(mutationId);
//                    payload.put("id", geneVariants.getGene().getEnsemblId());
//                    payload.put("geneSymbol", geneVariants.getGene().getGeneSymbol());
//                    payload.put("biotype", geneVariants.getGene().getBiotype());
//                    payload.put("ensemblId", geneVariants.getGene().getEnsemblId());
//                    payload.put("mutationId", mutationId);
//                    payload.put("variants", variants);
//                    payload.put("donors", new JSONArray());
//                    payload.put("frequencies", new new JSONArray());
//                    jsonObjectList.add(payload);
//
//                    if (jsonObjectList.size() >= bulkOpsQty) {
//                        //System.out.println("Bulk Items in partition-" + jsonObjectList.size());
//                        bulkStoreJsonObj(jsonObjectList, esUpsert, false, true);
//                        jsonObjectList = new ArrayList<>();
//                    }
//                }
//                bulkStoreJsonObj(jsonObjectList, esUpsert, false, true);
//            });

        }
        sc.close();
        client.close();
    }

    public static Boolean bulkStoreJsonObj(List<JSONObject> payloads,
                                           boolean esUpsert, boolean exoUpsert,
                                           boolean geneInsert) {

        if (payloads == null || payloads.isEmpty()) {
//            System.err.println("-");
            return null;
        }

        BulkRequest request = new BulkRequest();
//        List<JSONObject> geneList = new ArrayList<>();

        for (JSONObject payload: payloads) {
            if (payload == null ) {
                System.out.print("-");
                continue;
            }

            JSONArray labName = (JSONArray) payload.remove("labos");
            String uid = (String) payload.remove("id");

            if (geneInsert) {
                if (esUpsert) {
                    // bulk upsert gene with variants...
                    String mutationId = (String) payload.remove("mutationId");
                    //payload.put("uid", uid);
                    String newSymbol = (String) payload.remove("newSymbol");
                    String newAlias = (String) payload.remove("newAlias");
                    request.add(
                            upsertGeneRequest(payload.toString(0), uid, mutationId, newAlias)
                    );
                } else {
                    request.add(new IndexRequest(GENE_INDEX_NAME, "_doc", uid)
                            .source(payload.toString(0), XContentType.JSON));
                }
            } else if (exoUpsert) {
                payload.put("uid", uid);
                String specimenId = (String) payload.remove("specimenId");
                request.add(
                        upsertExomiserRequest(payload,uid,specimenId));
            } else if (esUpsert) {
                JSONArray donorArray = (JSONArray) payload.get("donors");
                JSONObject frequencies = (JSONObject) payload.get("frequencies");
                JSONArray specimenList = (JSONArray) payload.get("specimenList");

                request.add(
                        upsertMutationRequest(payload.toString(0), uid, donorArray, specimenList, frequencies, labName)
                );
            } else {
                request.add(new IndexRequest(MUTATION_INDEX_NAME, "_doc", uid)
                        .source(payload.toString(0), XContentType.JSON));
            }
        }
        boolean success = true;
        if (request.numberOfActions() > 0 ) {
            success = sendToES(request, (esUpsert||exoUpsert), false);
            if (!success) {
                // try splitting in third
                int qtyByThird = request.numberOfActions() / 3;
                int reminder = request.numberOfActions() % 3;
                int pos =0;
                boolean[] successArray = {false, false, false };
                for (int i=0; i<3; i++) {
                    BulkRequest thirdRequest = new BulkRequest();
                    List<DocWriteRequest<?>> requestList = request.requests();
                    if (reminder > 0 && i==2) {
                        qtyByThird += reminder;
                    }
                    for (int j=0; j<qtyByThird; j++) {
                        thirdRequest.add(requestList.get(pos++));
                    }
                    successArray[i] = sendToES(thirdRequest, (esUpsert||exoUpsert), true);

                }
                boolean retrySuccess = true;
                for (int i=0; i<3; i++) {
                    if (!successArray[i]) {
                        retrySuccess = false;
                        System.err.print("Unable to bulk " + ((esUpsert) ? "upsert " : "insert "));
                        for (JSONObject propertiesOneMutation : payloads) {
                            if (propertiesOneMutation == null) {
                                continue;
                            }
                            String mutationId = null;
                            if (propertiesOneMutation.isNull("mutationId")) {
                                JSONArray variants = (JSONArray) propertiesOneMutation.get("variants");
                                mutationId = variants.toString(0);
                            } else {
                                mutationId = (String) propertiesOneMutation.get("mutationId");
                            }
                            String uid = (String) propertiesOneMutation.remove("uid");
                            if (uid == null) {
                                uid = (String) propertiesOneMutation.get("ensemblId");
                            }
                            System.err.println(uid +"->" + mutationId);
                        }
                    }
                }
                if (retrySuccess) {
                    System.err.print("bulk " + ((esUpsert) ? "upsert " : "insert ") + "retried ok");
                }
            }
        }

        return success;
    }

    private static boolean sendToES(BulkRequest request, boolean esUpsert, boolean secondTry) {
        boolean indexingSuccess = false;
        BulkResponse bulkResponse = null;
        for (int i=0; i< 10; i++) {
            try {
                bulkResponse = client.bulk(request, RequestOptions.DEFAULT);
                if (bulkResponse != null && !bulkResponse.hasFailures()) {
                    indexingSuccess = true;
                }
                break;
            } catch (Exception e) {
                System.err.println("Error :"+e);
                if (e.getCause().getMessage().contains("Request Entity Too Large")) {
                    System.err.println("3x1/3");
                    break;
                }
                System.err.println("*******Bulk "+((esUpsert) ? "upsert":"insert")+" try #"+(i+1)+" failed...");
                if (secondTry) {
                    //System.err.println("bulResponse=" + bulkResponse.buildFailureMessage());
                    System.err.println("Error :"+e);
                }
            }
        }
        if (indexingSuccess) {
            if (secondTry) {
                System.out.print(".");
            } else {
                System.out.print("|");
            }
            return true;
        } else {
            return false;
        }
    }

    private static UpdateRequest upsertMutationRequest(String object, String uid,
                                                       JSONArray donors, JSONArray specimenList, JSONObject frequencies,
                                                       JSONArray labosName) {



        Map<String, Object> parameters = new HashMap<>();
        Map<String, Object> donorMap = null;//= new HashMap<>();
        Map<String, Object> laboMap = null;
        List<Map<String, Object>> labFreqList = new ArrayList<>();
        List<String> specimenLst = new ArrayList<>();
        List<String> labNameLst = new ArrayList<>();
        List<Map<String, Object>> donorsLst = new ArrayList<>();
        for (int i=0; i<donors.length(); i++) {
            donorMap = new HashMap<>();
            specimenLst.add( (String) specimenList.get(i) );
            JSONObject donor = (JSONObject) donors.get(i);
//            donorMap.put("depth", donor.get("depth"));
//            donorMap.put("mq", donor.get("mq"));
            donorMap.put("filter", donor.get("filter"));
            donorMap.put("specimenId", specimenList.get(i));
            donorMap.put("patientId", donor.get("patientId"));
            donorMap.put("familyId", donor.get("familyId"));
            donorMap.put("relation", donor.get("relation"));
            donorMap.put("sequencingStrategy", donor.get("sequencingStrategy"));
            donorMap.put("studyId", donor.get("studyId"));
            donorMap.put("zygosity", donor.get("zygosity"));
            //donorMap.put("ad", donor.get("ad"));
            donorMap.put("adFreq", donor.get("adFreq"));
            donorMap.put("adAlt", donor.get("adAlt"));
            donorMap.put("adTotal", donor.get("adTotal"));
            //donorMap.put("exomiserScore", donor.get("exomiserScore"));
//            donorMap.put("af", donor.get("af"));
//            donorMap.put("dp", donor.get("dp"));
            donorMap.put("gt", donor.get("gt"));
            if (!donor.isNull("qd")) {
                donorMap.put("qd", donor.get("qd"));
            }
            if (!donor.isNull("gq")) {
                donorMap.put("gq", donor.get("gq"));
            }
            donorMap.put("practitionerId", donor.get("practitionerId"));
            donorMap.put("organizationId", donor.get("organizationId"));
            if (!donor.isNull("genotypeFamily")) {
                donorMap.put("genotypeFamily", donor.get("genotypeFamily"));
            }
            if (!donor.isNull("transmission")) {
                donorMap.put("transmission", donor.get("transmission"));
            }
//            if (!donor.isNull("dn")) {
//                donorMap.put("dn", donor.get("dn"));
//            }
//            if (!donor.isNull("dq")) {
//                donorMap.put("dq", donor.get("dq"));
//            }
            if (!donor.isNull("nbHpoTerms")) {
                donorMap.put("nbHpoTerms", donor.get("nbHpoTerms"));
            }

            donorMap.put("lastUpdate", donor.get("lastUpdate"));

            donorsLst.add(donorMap);
        }

        for (int i=0; i<labosName.length(); i++) {


            String labName = (String) labosName.get(i);
            labNameLst.add(labName);
            JSONObject freqLabo = (JSONObject) frequencies.get(labName);
            laboMap = new HashMap<>();

            laboMap.put("AC", freqLabo.get("AC"));
            laboMap.put("AN", freqLabo.get("AN"));
            laboMap.put("AF", freqLabo.get("AF"));
            laboMap.put("HC", freqLabo.get("HC"));
            laboMap.put("PN", freqLabo.get("PN"));
            laboMap.put("labName", labName);
            labFreqList.add(laboMap);

        }


        JSONObject freqInterne = (JSONObject) frequencies.get("interne");


        parameters.put("specimen", specimenList);
        parameters.put("donorsLst",  donorsLst );
        //parameters.put("donorMap", donorMap);
        parameters.put("AC", freqInterne.get("AC"));
        parameters.put("AN", freqInterne.get("AN"));
        parameters.put("HC", freqInterne.get("HC"));
        parameters.put("PN", freqInterne.get("PN"));
        parameters.put("labNameLst", labNameLst);
        parameters.put("labFreqList", labFreqList);
        //parameters.put("")


//        Map<String, Object> freqLaboMap = new HashMap<>();
//        freqLaboMap.put("AC", freqLabo.get("AC"));
//        freqLaboMap.put("AN", freqLabo.get("AN"));
//        freqLaboMap.put("AF", freqLabo.get("AF"));
//        freqLaboMap.put("HC", freqLabo.get("HC"));
//        freqLaboMap.put("PN", freqLabo.get("PN"));
        parameters.put("freqLabList", labFreqList);


        Script inline = new Script(ScriptType.INLINE, "painless",

                "boolean toUpdateFreq = false; " +
                        "for (int i=0; i<params.specimen.size(); i++) {" +
                            "String s = params.specimen.get(i); " +
                            "Map d = params.donorsLst.get(i); " +
                            "if (!ctx._source.specimenList.contains(s) && d != null) {" +
                                "ctx._source.donors.add(d);" +
                                "ctx._source.specimenList.add(s);" +
                                "toUpdateFreq = true " +
                            "} " +
                        "} " +
                        "if (toUpdateFreq) {" +
                            "ctx._source.frequencies.interne.AC += params.AC; " +
                            "ctx._source.frequencies.interne.AN += params.AN; " +
                            "ctx._source.frequencies.interne.HC += params.HC; " +
                            "ctx._source.frequencies.interne.PN += params.PN; " +
                            "if (ctx._source.frequencies.interne.AN > 0) {" +
                                " ctx._source.frequencies.interne.AF = (float) ctx._source.frequencies.interne.AC / ctx._source.frequencies.interne.AN " +
                            "}" +
                        "} " +
                        "for (int i=0; i<params.freqLabList.size(); i++) {" +
                            "Map freqLabo = params.freqLabList.get(i);" +
                            "String labNameToUpdate = freqLabo.get(\"labName\");" +
                            "boolean freqLaboExist = false; " +
                            "if (ctx._source.frequencies.containsKey(labNameToUpdate)) " +
                                "{ freqLaboExist = true } " +
                            "if (toUpdateFreq && freqLaboExist) { " +
                                "Map theFreqToUpdate = ctx._source.frequencies.get(labNameToUpdate); " +
                                "theFreqToUpdate.AC += freqLabo.AC; " +
                                "theFreqToUpdate.AN += freqLabo.AN; " +
                                "theFreqToUpdate.HC += freqLabo.HC; " +
                                "theFreqToUpdate.PN += freqLabo.PN; " +
                                "if (theFreqToUpdate.AN > 0) {" +
                                    " theFreqToUpdate.AF = " +
                                        "(float) theFreqToUpdate.AC / theFreqToUpdate.AN }" +
                                "} " +
                                "if (toUpdateFreq && !freqLaboExist) {" +
                                    "ctx._source.frequencies.put(labNameToUpdate,freqLabo)" +
                            "}" +
                        "}"
                , parameters);


        UpdateRequest request = new UpdateRequest(MUTATION_INDEX_NAME, "_doc", uid);
        request.script(inline);

        request.upsert(object, XContentType.JSON);
        return request;
    }

    protected static UpdateRequest upsertExomiserRequest(JSONObject object, String uid, String specimenId) {


        Map<String, Object> parameters = new HashMap<>();

        parameters.put("transmission", object.get("transmission"));
        parameters.put("exomiserScore", object.get("combinedScore"));
        parameters.put("lastUpdate", object.get("lastUpdate"));
        //parameters.put("transmissionKeyName", "transmission");
        parameters.put("exomiserScoreKeyName", "exomiserScore");
        parameters.put("overwriteKey", "lastUpdate");
        parameters.put("specimen", specimenId);


        Script inline = new Script(ScriptType.INLINE, "painless",

        "" +

                "String specimen = params.specimen; " +
                "if (ctx._source.specimenList != null && ctx._source.specimenList.contains(specimen) ) { " +
                    "for (Map donor : ctx._source.donors) { " +
                        "if (donor.specimenId == specimen) { " +
                            "donor.put(params.exomiserScoreKeyName,params.exomiserScore); " +
//                            "donor.put(params.transmissionKeyName,params.transmission); " +
                            "donor.put(params.overwriteKey,params.lastUpdate); " +
                        "} " +
                    "} " +
                "} "
                , parameters);


        UpdateRequest request = new UpdateRequest(MUTATION_INDEX_NAME, "_doc", uid);
        request.script(inline);

        request.upsert(object.toString(0), XContentType.JSON);
        return request;
    }

    private static UpdateRequest upsertGeneRequest(String object, String uid,
                                               JSONArray donors, JSONArray specimenList, JSONObject frequencies,
                                               JSONArray labosName, JSONArray mutations) {



        Map<String, Object> parameters = new HashMap<>();
        Map<String, Object> donorMap = null;//= new HashMap<>();
        Map<String, Object> laboMap = null;
        List<Map<String, Object>> labFreqList = new ArrayList<>();
        List<String> specimenLst = new ArrayList<>();
        List<String> mutationList = new ArrayList<>();
        List<String> labNameLst = new ArrayList<>();
        List<Map<String, Object>> donorsLst = new ArrayList<>();
        for (int i=0; i<donors.length(); i++) {
            donorMap = new HashMap<>();
            specimenLst.add( (String) specimenList.get(i) );
            JSONObject donor = (JSONObject) donors.get(i);

            donorMap.put("filter", donor.get("filter"));
            donorMap.put("specimenId", specimenList.get(i));
            donorMap.put("patientId", donor.get("patientId"));
            donorMap.put("familyId", donor.get("familyId"));
            donorMap.put("relation", donor.get("relation"));
            donorMap.put("sequencingStrategy", donor.get("sequencingStrategy"));
            donorMap.put("studyId", donor.get("studyId"));
            donorMap.put("zygosity", donor.get("zygosity"));

            donorMap.put("adFreq", donor.get("adFreq"));
            donorMap.put("adAlt", donor.get("adAlt"));
            donorMap.put("adTotal", donor.get("adTotal"));

            donorMap.put("gt", donor.get("gt"));
            if (!donor.isNull("qd")) {
                donorMap.put("qd", donor.get("qd"));
            }
            if (!donor.isNull("gq")) {
                donorMap.put("gq", donor.get("gq"));
            }
            donorMap.put("practitionerId", donor.get("practitionerId"));
            donorMap.put("organizationId", donor.get("organizationId"));
            if (!donor.isNull("genotypeFamily")) {
                donorMap.put("genotypeFamily", donor.get("genotypeFamily"));
            }

            if (!donor.isNull("nbHpoTerms")) {
                donorMap.put("nbHpoTerms", donor.get("nbHpoTerms"));
            }
            donorMap.put("lastUpdate", donor.get("lastUpdate"));

            donorsLst.add(donorMap);
        }

        for (int i=0; i<labosName.length(); i++) {

            String labName = (String) labosName.get(i);
            labNameLst.add(labName);
            JSONObject freqLabo = (JSONObject) frequencies.get(labName);
            laboMap = new HashMap<>();

            laboMap.put("AC", freqLabo.get("AC"));
            laboMap.put("AN", freqLabo.get("AN"));
            laboMap.put("AF", freqLabo.get("AF"));
            laboMap.put("HC", freqLabo.get("HC"));
            laboMap.put("PN", freqLabo.get("PN"));
            laboMap.put("labName", labName);
            labFreqList.add(laboMap);
        }

        JSONObject freqInterne = (JSONObject) frequencies.get("interne");

        parameters.put("specimen", specimenList);
        parameters.put("donorsLst",  donorsLst );
        parameters.put("mutation", mutations);

        //parameters.put("donorMap", donorMap);
        parameters.put("AC", freqInterne.get("AC"));
        parameters.put("AN", freqInterne.get("AN"));
        parameters.put("HC", freqInterne.get("HC"));
        parameters.put("PN", freqInterne.get("PN"));
        parameters.put("labNameLst", labNameLst);
        parameters.put("labFreqList", labFreqList);

        parameters.put("freqLabList", labFreqList);


        Script inline = new Script(ScriptType.INLINE, "painless",

                ""+
                        "boolean toUpdateFreq = false; " +
                        "for (int i=0; i<params.specimen.size(); i++) {" +
                            "String s = params.specimen.get(i); " +
                            "Map d = params.donorsLst.get(i); " +
                            "if (!ctx._source.specimenList.contains(s) && d != null) {" +
                                "ctx._source.donors.add(d);" +
                                "ctx._source.specimenList.add(s);" +
                                "toUpdateFreq = true " +
                            "} " +
                        "} " +
                        "if (toUpdateFreq) {" +
                            "ctx._source.frequencies.interne.AC += params.AC; " +
                            "ctx._source.frequencies.interne.AN += params.AN; " +
                            "ctx._source.frequencies.interne.HC += params.HC; " +
                            "ctx._source.frequencies.interne.PN += params.PN; " +
                                "if (ctx._source.frequencies.interne.AN > 0) {" +
                                " ctx._source.frequencies.interne.AF = (float) ctx._source.frequencies.interne.AC / ctx._source.frequencies.interne.AN " +
                                "}" +
                        "} " +
                        "for (int i=0; i<params.freqLabList.size(); i++) {" +
                        "Map freqLabo = params.freqLabList.get(i);" +
                        "String labNameToUpdate = freqLabo.get(\"labName\");" +
                        "boolean freqLaboExist = false; " +
                        "if (ctx._source.frequencies.containsKey(labNameToUpdate)) " +
                            "{ freqLaboExist = true } " +
                        "if (toUpdateFreq && freqLaboExist) { " +
                            "Map theFreqToUpdate = ctx._source.frequencies.get(labNameToUpdate); " +
                            "theFreqToUpdate.AC += freqLabo.AC; " +
                            "theFreqToUpdate.AN += freqLabo.AN; " +
                            "theFreqToUpdate.HC += freqLabo.HC; " +
                            "theFreqToUpdate.PN += freqLabo.PN; " +
                            "if (theFreqToUpdate.AN > 0) {" +
                                " theFreqToUpdate.AF = " +
                                "(float) theFreqToUpdate.AC / theFreqToUpdate.AN }" +
                            "} " +
                            "if (toUpdateFreq && !freqLaboExist) {" +
                                "ctx._source.frequencies.put(labNameToUpdate,freqLabo)" +
                            "}" +
                        "}"
                , parameters);

        UpdateRequest request = new UpdateRequest(GENE_INDEX_NAME, "_doc", uid);
        request.script(inline);

        request.upsert(object, XContentType.JSON);
        return request;
    }

    private static UpdateRequest upsertGeneRequest(String object, String uid, String mutationId, String newAlias) {


        Map<String, Object> parameters = new HashMap<>();
        Map<String, Object> aliasMap = new HashMap<>();

        parameters.put("mutationId", mutationId);
        parameters.put("newAlias", newAlias);
        parameters.put("aliasKeyName", "alias");
        String[] aliasArray = new String[]{newAlias};

        aliasMap.put("alias", aliasArray);
        parameters.put("alias", aliasMap);


        Script inline = new Script(ScriptType.INLINE, "painless",

            "" +
                    "String mutationId = params.mutationId; " +
                    "if (!ctx._source.variants.contains(mutationId)) {" +
                           " ctx._source.variants.add(mutationId);" +
                    "} " +
                    "if (params.newAlias != null && !ctx._source.alias.contains(params.newAlias)) { " +
                        "ctx._source.alias.add(params.newAlias);" +
                    "} "
                , parameters);

        UpdateRequest request = new UpdateRequest(GENE_INDEX_NAME, "_doc", uid);
        request.script(inline);

        request.upsert(object, XContentType.JSON);
        return request;
    }


    public static String getSHA256Hash(String data) {
        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            byte[] hash = digest.digest(data.getBytes(StandardCharsets.UTF_8));
            return bytesToHex(hash);
        } catch(Exception ex) {
            ex.printStackTrace();
        }
        return null;
    }

    static String getSHA1Hash(String data) {
        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-1");
            byte[] hash = digest.digest(data.getBytes(StandardCharsets.UTF_8));
            return bytesToHex(hash);
        } catch(Exception ex) {
            ex.printStackTrace();
        }
        return null;
    }

    static String getMD5Hash(String data) {
        try {
            MessageDigest digest = MessageDigest.getInstance("MD5");
            byte[] hash = digest.digest(data.getBytes(StandardCharsets.UTF_8));
            return bytesToHex(hash);
        } catch(Exception ex) {
            ex.printStackTrace();
        }
        return null;
    }

    private static String bytesToHex(byte[] hash) {
        return DatatypeConverter.printHexBinary(hash);
    }

    private static JSONArray toJson(String[] bob) {
        JSONArray array = new JSONArray();
        for (String line: bob) {
            JSONObject obj = new JSONObject(line);
            array.put(obj);
        }
        return array;

    }

    private void testFreqLab(List<Map<String, Object>> labFreqList, Map<String, Object> theMapToUpdate) {
        for (int i=0; i< labFreqList.size(); i++) {
            Map<String, Object> freqLabo = labFreqList.get(i);
            String labName = (String) freqLabo.get("labName");
            if (theMapToUpdate.containsKey(labName)) {
                JSONObject theObj = (JSONObject) theMapToUpdate.get(labName);
                theObj.put("AN", (String) freqLabo.get("labName"));
                theObj.remove("labName");
            }
        }
    }
}
