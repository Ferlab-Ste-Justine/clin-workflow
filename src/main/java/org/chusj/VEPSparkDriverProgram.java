package org.chusj;

import com.whalin.MemCached.MemCachedClient;
import com.whalin.MemCached.SockIOPool;
import org.apache.http.HttpHost;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.ForeachPartitionFunction;
import org.apache.spark.api.java.function.Function2;
import org.elasticsearch.action.bulk.BulkRequest;
import org.elasticsearch.action.bulk.BulkResponse;
import org.elasticsearch.action.get.GetRequest;
import org.elasticsearch.action.get.GetResponse;
import org.elasticsearch.action.index.IndexRequest;
import org.elasticsearch.action.update.UpdateRequest;
import org.elasticsearch.action.update.UpdateResponse;
import org.elasticsearch.client.RequestOptions;
import org.elasticsearch.client.RestClient;
import org.elasticsearch.client.RestHighLevelClient;
import org.elasticsearch.common.xcontent.XContentType;
import org.elasticsearch.script.Script;
import org.elasticsearch.script.ScriptType;
import org.json.JSONArray;
import org.json.JSONObject;

import javax.xml.bind.DatatypeConverter;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.util.*;

import static org.chusj.VepHelper.addDonorArrayToDonorArray;

public class VEPSparkDriverProgram {

    public static RestHighLevelClient client;
    private static SockIOPool pool;
    private static MemCachedClient mcc;
    private static String MEMCACHED_SERVER = "localhost:11212";
    private static String INDEX_NAME = "mutations";


    public static void main(String[] args) throws Exception {

        if (args.length < 7 ) {
            throw new Exception("Missing params; need extractFile pedigreePropertiesFile ES_UPSERT (true|false) " +
                    "sparkMaster local[nbWorkers] 8g nbPartitions");
        }

        String extractFile = args[0];
        String pedigreePropsFile = args[1];

        boolean esUpsert = Boolean.parseBoolean(args[2]);
        String sparkMaster = args[3];
        String localThread = args[4];
        String memory = args[5];
        int nbPartitions = Integer.parseInt(args[6]);
        int bulkOpsQty = Integer.parseInt(args[7]);
//        int localTi = Integer.valueOf(localThread);

        /* Define Spark Configuration */
        SparkConf conf = new SparkConf().setAppName("ExtractTLoad").setMaster(sparkMaster)
                .setMaster(localThread).set("spark.executor.memory",memory);

        /* Create Spark Context with configuration */
        JavaSparkContext sc = new JavaSparkContext(conf);

        /* Create a Resilient Distributed Dataset for a long file
         * Each line in long file become a record in RDD
         * */
        JavaRDD<String> lines = sc.textFile(extractFile, nbPartitions);

        Properties pedigreeProps = VepHelper.getPropertiesFromFile(pedigreePropsFile);
        pedigreeProps.forEach( (k,v) -> System.out.println(k+"="+v) );

        try (RestHighLevelClient clientTry = new RestHighLevelClient(
                RestClient.builder(
                        new HttpHost("localhost", 9200, "http")))) {

            client = clientTry;
//            connectToMemCached();

            //lines.foreach((a) -> buildAndStoreJsonObj(a, esUpsert, pedigreeProps));

            //JavaRDD<JSONObject> listJsonObjs = lines.map((a) -> buildJsonObj(a,pedigreeProps));
            //JavaRDD<List<JSONObject>> listOfJsonObj = lines.map((a) -> buildJsonObj(a,pedigreeProps));
//            JavaRDD<JSONObject> collection = lines.mapPartitions( (Iterator<String> iter1) -> {
//                ArrayList<JSONObject> out = new ArrayList<>();
//                while(iter1.hasNext()) {
//                    String value = iter1.next();
//                    JSONObject current = buildJsonObj(value, pedigreeProps);
//                    out.add(current);
//                }
//                return out.iterator();
//            });

//            collection.foreach((listOfJsonObjs) -> {
//                if (listOfJsonObjs != null) {
////                    System.out.print("-"+listOfJsonObjs.length()+"-");
//                    //bulkStoreJsonObj(listOfJsonObjs,esUpsert,pedigreeProps);
//                }
//            });


            lines.foreachPartition(partitionOfRecords -> {
                List<JSONObject> jsonObjectList = new ArrayList<>();
                while (partitionOfRecords.hasNext()) {
                    if (jsonObjectList.size() >= bulkOpsQty) {
                        // empty and bulk updates
                        System.out.println("Bulk Items in partition-" + jsonObjectList.size());
                        bulkStoreJsonObj(jsonObjectList,esUpsert,pedigreeProps);
                        jsonObjectList = new ArrayList<>();
                    }
                    jsonObjectList.add(buildJsonObj(partitionOfRecords.next(), pedigreeProps));
                }
                if (jsonObjectList.size() > 0) {
                    System.out.println("Bulk Items in partition-" + jsonObjectList.size());
                    bulkStoreJsonObj(jsonObjectList,esUpsert,pedigreeProps);
                }
            });
        }
        sc.close();
        client.close();
    }




    private static JSONObject buildJsonObj(String extractedLine, Properties pedigreeProps) {
        return  VepHelper.processVcfDataLine(extractedLine, "dn,dq", pedigreeProps);
    }


    private static void buildAndStoreJsonObj(String extractedLine, boolean esUpsert, Properties pedigreeProps) throws IOException {


        JSONObject propertiesOneMutation = VepHelper.processVcfDataLine(extractedLine, "dn,dq", pedigreeProps);

        // extract donor info if not found
        //String qual = (String) propertiesOneMutation.remove("qual");
        //String filter = (String) propertiesOneMutation.remove("filter");
        if (propertiesOneMutation == null) return;
        //JSONArray donorArray = (JSONArray) propertiesOneMutation.remove("donors");
        JSONArray donorArray = (JSONArray) propertiesOneMutation.get("donors");

        JSONArray newDonorArray;
        // verify if variant already exist
        boolean toMapping = false;
        boolean toIndex = true;
        boolean toMemcached = false;
        boolean checkES = false;
        boolean checkMemcached = false;
        String msg= "";



        String uid = (String) propertiesOneMutation.remove("id");
//        if ("B5C46555D00EC3C0638D5B84535185A1".equalsIgnoreCase(uid)) {
//            System.err.println(propertiesOneMutation.toString(0));
//        }
        String dnaChanges = (String) propertiesOneMutation.get("mutationId");

        String specimenIdList = pedigreeProps.getProperty("pedigree");
        String laboName = pedigreeProps.getProperty("laboName");

        //String uid = getSHA1Hash(familyId+ "@" + dnaChanges);

        JSONArray previousDonorArray = null;
        JSONArray specimentList;


        if (checkMemcached) {


            String donorArrayStr = (String) mcc.get(uid);
            if (donorArrayStr != null) {
                previousDonorArray = new JSONArray(donorArrayStr);
                if (previousDonorArray != null && previousDonorArray.length() > 0) {

                    boolean donorFound = checkForSpecimen(donorArray, specimenIdList);

                    if (!donorFound) {
                        // add new donor(s) to previous one
                        msg += "m0";
                        //donorArray.put(newDonor);
                        //previousDonorArray.put(donorArray);
                        addDonorArrayToDonorArray(previousDonorArray, donorArray);
                        toIndex = true;
                        toMemcached = true;

                    } else {
                        msg += "m1"; //already present - no update needed - no check ES

                    }
                    System.out.print(msg);
                } else {
                    //donorArray = new JSONArray();
                    //donorArray.put(newDonor);
                    toMemcached = true;
                    checkES = true;
                }
            } else {
                // might need to check ES
                toMemcached = true;
                checkES = true;
            }
        }

        if (checkES) {
            GetRequest getRequest = new GetRequest("mutations", "_doc", uid);

            if (client.exists(getRequest, RequestOptions.DEFAULT)) {

                GetResponse getResponse = client.get(getRequest, RequestOptions.DEFAULT);

                if (getResponse.isExists()) {

                    JSONObject obj = new JSONObject(getResponse.getSourceAsString());

                    //JSONObject props = (JSONObject) obj.get("_source");

                    String mutationId = (String) obj.get("mutationId");
                    //msg += mutationId;

                    if (!mutationId.equalsIgnoreCase(dnaChanges)) {
                        // collision detected
                        msg = "(" + dnaChanges + ";" + uid + ";";
                        System.err.println("collision="+msg);
                    }
                    try {
                        previousDonorArray = (JSONArray) obj.get("donors");
                    } catch (Exception e) {
                        previousDonorArray = new JSONArray();
                    }

                    boolean donorFound = checkForSpecimen(previousDonorArray, specimenIdList);

                    if (!donorFound) {
                        // add new donor to previous one
                        msg += "e0";
                        //System.out.print("0)");
                        addDonorArrayToDonorArray(previousDonorArray, donorArray);
                        toIndex = true;

                    } else {
                        msg += "e1";
                        // nothing - might have to update
                        if (toMemcached) {
                            addDonorArrayToDonorArray(previousDonorArray, donorArray);
                        }
                    }
                    System.out.print(msg);
                } else { // mutation not found

                    previousDonorArray = donorArray;
                    toIndex = true;
                }
            } else { // index not created yet
                toMapping = true;
                toIndex = true;
                previousDonorArray = donorArray;

            }
        }

        if (toMapping) {
            // have to put mapping TODO (currently done manually on Kibana)
        }
        boolean indexingSuccess = false;

        if (toIndex) {
//            if (previousDonorArray == null) {
//                previousDonorArray = donorArray;
//            }

//            previousDonorArray = (JSONArray) propertiesOneMutation.get("donors");
//            JSONObject frequencies = (JSONObject) propertiesOneMutation.get("frequencies");
//            specimentList = (JSONArray) propertiesOneMutation.get("specimenList");
            for (int i=0; i< 3; i++) {
                try {
                    if (esUpsert) {
                        previousDonorArray = (JSONArray) propertiesOneMutation.get("donors");
                        JSONObject frequencies = (JSONObject) propertiesOneMutation.get("frequencies");
                        specimentList = (JSONArray) propertiesOneMutation.get("specimenList");
                        upsert( upsertRequest(propertiesOneMutation.toString(), uid, previousDonorArray, specimentList, frequencies, laboName), uid);
                    } else {
                        index(propertiesOneMutation.toString(), client, uid, INDEX_NAME); // overwrite index
                    }
                    indexingSuccess = true;
                    break;
                } catch (Exception e) {
                    System.err.println("*******Indexation try #"+i+" failed..."+e);
                }
            }

            if (!indexingSuccess) System.err.println("#########\n\n\n######### Unable to index " + uid + "\n############");

        }
        if (toMemcached && indexingSuccess) {
//            if (previousDonorArray == null) {
//                previousDonorArray = donorArray;
//            }
            mcc.set(uid, previousDonorArray.toString());
        }

    }


    private static void storeJsonObj(JSONObject propertiesOneMutation, boolean esUpsert, Properties pedigreeProps) throws IOException {


//        JSONObject propertiesOneMutation = VepHelper.processVcfDataLine(extractedLine, "dn,dq", pedigreeProps);
//
        // extract donor info if not found
        //String qual = (String) propertiesOneMutation.remove("qual");
        //String filter = (String) propertiesOneMutation.remove("filter");
        if (propertiesOneMutation == null) return;
        //JSONArray donorArray = (JSONArray) propertiesOneMutation.remove("donors");
        JSONArray donorArray = (JSONArray) propertiesOneMutation.get("donors");

        JSONArray newDonorArray;
        // verify if variant already exist
        boolean toMapping = false;
        boolean toIndex = true;
        boolean toMemcached = false;
        boolean checkES = false;
        boolean checkMemcached = false;
        String msg= "";



        String uid = (String) propertiesOneMutation.remove("id");
//        if ("B5C46555D00EC3C0638D5B84535185A1".equalsIgnoreCase(uid)) {
//            System.err.println(propertiesOneMutation.toString(0));
//        }
        String dnaChanges = (String) propertiesOneMutation.get("mutationId");

        String specimenIdList = pedigreeProps.getProperty("pedigree");
        String laboName = pedigreeProps.getProperty("laboName");

        //String uid = getSHA1Hash(familyId+ "@" + dnaChanges);

        JSONArray previousDonorArray = null;
        JSONArray specimenList;


        if (checkMemcached) {


            String donorArrayStr = (String) mcc.get(uid);
            if (donorArrayStr != null) {
                previousDonorArray = new JSONArray(donorArrayStr);
                if (previousDonorArray != null && previousDonorArray.length() > 0) {

                    boolean donorFound = checkForSpecimen(donorArray, specimenIdList);

                    if (!donorFound) {
                        // add new donor(s) to previous one
                        msg += "m0";
                        //donorArray.put(newDonor);
                        //previousDonorArray.put(donorArray);
                        addDonorArrayToDonorArray(previousDonorArray, donorArray);
                        toIndex = true;
                        toMemcached = true;

                    } else {
                        msg += "m1"; //already present - no update needed - no check ES

                    }
                    System.out.print(msg);
                } else {
                    //donorArray = new JSONArray();
                    //donorArray.put(newDonor);
                    toMemcached = true;
                    checkES = true;
                }
            } else {
                // might need to check ES
                toMemcached = true;
                checkES = true;
            }
        }

        if (checkES) {
            GetRequest getRequest = new GetRequest("mutations", "_doc", uid);

            if (client.exists(getRequest, RequestOptions.DEFAULT)) {

                GetResponse getResponse = client.get(getRequest, RequestOptions.DEFAULT);

                if (getResponse.isExists()) {

                    JSONObject obj = new JSONObject(getResponse.getSourceAsString());

                    //JSONObject props = (JSONObject) obj.get("_source");

                    String mutationId = (String) obj.get("mutationId");
                    //msg += mutationId;

                    if (!mutationId.equalsIgnoreCase(dnaChanges)) {
                        // collision detected
                        msg = "(" + dnaChanges + ";" + uid + ";";
                        System.err.println("collision="+msg);
                    }
                    try {
                        previousDonorArray = (JSONArray) obj.get("donors");
                    } catch (Exception e) {
                        previousDonorArray = new JSONArray();
                    }

                    boolean donorFound = checkForSpecimen(previousDonorArray, specimenIdList);

                    if (!donorFound) {
                        // add new donor to previous one
                        msg += "e0";
                        //System.out.print("0)");
                        addDonorArrayToDonorArray(previousDonorArray, donorArray);
                        toIndex = true;

                    } else {
                        msg += "e1";
                        // nothing - might have to update
                        if (toMemcached) {
                            addDonorArrayToDonorArray(previousDonorArray, donorArray);
                        }
                    }
                    System.out.print(msg);
                } else { // mutation not found

                    previousDonorArray = donorArray;
                    toIndex = true;
                }
            } else { // index not created yet
                toMapping = true;
                toIndex = true;
                previousDonorArray = donorArray;

            }
        }

        if (toMapping) {
            // have to put mapping TODO (currently done manually on Kibana)
        }
        boolean indexingSuccess = false;

        if (toIndex) {
//            if (previousDonorArray == null) {
//                previousDonorArray = donorArray;
//            }

//            previousDonorArray = (JSONArray) propertiesOneMutation.get("donors");
//            JSONObject frequencies = (JSONObject) propertiesOneMutation.get("frequencies");
//            specimentList = (JSONArray) propertiesOneMutation.get("specimenList");
            for (int i=0; i< 3; i++) {
                try {
                    if (esUpsert) {
                        previousDonorArray = (JSONArray) propertiesOneMutation.get("donors");
                        JSONObject frequencies = (JSONObject) propertiesOneMutation.get("frequencies");
                        specimenList = (JSONArray) propertiesOneMutation.get("specimenList");
                        upsert( upsertRequest(propertiesOneMutation.toString(), uid, previousDonorArray, specimenList, frequencies, laboName), uid);
                    } else {
                        index(propertiesOneMutation.toString(), client, uid, INDEX_NAME); // overwrite index
                    }
                    indexingSuccess = true;
                    break;
                } catch (Exception e) {
                    System.err.println("*******Indexation try #"+i+" failed..."+e);
                }
            }

            if (!indexingSuccess) System.err.println("#########\n\n\n######### Unable to index " + uid + "\n############");

        }
        if (toMemcached && indexingSuccess) {
//            if (previousDonorArray == null) {
//                previousDonorArray = donorArray;
//            }
            mcc.set(uid, previousDonorArray.toString());
        }

    }


    private static Boolean bulkStoreJsonObj(List<JSONObject> propertiesOneMutations, boolean esUpsert, Properties pedigreeProps) {

        if (propertiesOneMutations == null) return null;

        BulkRequest request = new BulkRequest();

        for (JSONObject propertiesOneMutation: propertiesOneMutations) {
            if (propertiesOneMutation == null ) {
                continue;
            }

            String uid = (String) propertiesOneMutation.remove("id");
            if (esUpsert) {
                JSONArray donorArray = (JSONArray) propertiesOneMutation.get("donors");
                JSONObject frequencies = (JSONObject) propertiesOneMutation.get("frequencies");
                JSONArray specimenList = (JSONArray) propertiesOneMutation.get("specimenList");
                String laboName = pedigreeProps.getProperty("laboName");

                request.add(
                        upsertRequest(propertiesOneMutation.toString(0), uid, donorArray, specimenList, frequencies, laboName)
                );
            } else {
                request.add(new IndexRequest("mutations", "_doc", uid)
                        .source(propertiesOneMutation.toString(0), XContentType.JSON));
            }
        }

        boolean indexingSuccess = false;
        //BulkResponse bulkResponse;
        for (int i=0; i< 10; i++) {
            try {
                //bulkResponse =
                        client.bulk(request, RequestOptions.DEFAULT);
                indexingSuccess = true;
                break;
            } catch (Exception e) {
                System.err.println("*******Bulk "+((esUpsert) ? "upsert":"insert")+" try #"+i+" failed...\n"+e);
            }
        }
        if (indexingSuccess) {
            System.out.print(".");
            return true;
        } else {
            System.err.print("Unable to bulk "+((esUpsert) ? "upsert ":"insert "));
            for (JSONObject propertiesOneMutation: propertiesOneMutations) {
                if (propertiesOneMutation == null) {
                    continue;
                }
                System.err.println("->"+propertiesOneMutation.get("mutationId"));
            }
            return false;
        }
    }

    private static void index(String object, RestHighLevelClient client, String uid, String index) throws Exception {
        IndexRequest request = new IndexRequest("mutations", "_doc", uid);

        request.source(object, XContentType.JSON);
//        System.out.println("\n"+object);
        //IndexResponse indexResponse;
        try {
            client.index(request, RequestOptions.DEFAULT);
            //System.out.print("."+indexResponse.status().getStatus()+".");
        } catch (Exception e) {
            System.err.println("E="+e.getMessage()+"-"+uid);
            try {
                Thread.sleep(200);
            } catch (InterruptedException e1) {
                e1.printStackTrace();
            }
            //indexResponse = client.index(request, RequestOptions.DEFAULT);
            throw new Exception(e);
        }

        System.out.print(".");
    }

    private static void upsert(UpdateRequest updateRequest, String uid) throws Exception {

        UpdateResponse result = null;
        try {
            result = client.update(updateRequest, RequestOptions.DEFAULT);
            //System.out.print("."+indexResponse.status().getStatus()+".");
            //System.err.println("result="+result);
        } catch (Exception e) {
            System.err.println("E="+e.getMessage()+"-"+uid);
            System.err.println("result="+result+e);
            //indexResponse = client.index(request, RequestOptions.DEFAULT);
            throw new Exception(e);
        }

        System.out.print(".");
    }


    private static UpdateRequest upsertRequest(String object, String uid,
                               JSONArray donors, JSONArray specimenList, JSONObject frequencies,
                               String laboName) {



        Map<String, Object> parameters = new HashMap<>();
        Map<String, Object> donorMap = null;//= new HashMap<>();
        List<String> specimenLst = new ArrayList<>();
        List<Map<String, Object>> donorsLst = new ArrayList<>();
        for (int i=0; i<donors.length(); i++) {
            donorMap = new HashMap<>();
            specimenLst.add( (String) specimenList.get(i) );
            JSONObject donor = (JSONObject) donors.get(i);
            donorMap.put("depth", donor.get("depth"));
            donorMap.put("mq", donor.get("mq"));
            donorMap.put("filter", donor.get("filter"));
            donorMap.put("specimenId", specimenList.get(i));
            donorMap.put("patientId", donor.get("patientId"));
            donorMap.put("familyId", donor.get("familyId"));
            donorMap.put("relation", donor.get("relation"));
            donorMap.put("sequencingStrategy", donor.get("sequencingStrategy"));
            donorMap.put("studyId", donor.get("studyId"));
            donorMap.put("zygosity", donor.get("zygosity"));
            donorMap.put("ad", donor.get("ad"));
            donorMap.put("af", donor.get("af"));
            donorMap.put("dp", donor.get("dp"));
            donorMap.put("gt", donor.get("gt"));
            if (!donor.isNull("gq")) {
                donorMap.put("gq", donor.get("gq"));
            }
            donorMap.put("practitionerId", donor.get("practitionerId"));
            donorMap.put("organizationId", donor.get("organizationId"));
            if (!donor.isNull("genotypeFamily")) {
                donorMap.put("genotypeFamily", donor.get("genotypeFamily"));
            }
            if (!donor.isNull("dn")) {
                donorMap.put("dn", donor.get("dn"));
            }
            if (!donor.isNull("dq")) {
                donorMap.put("dq", donor.get("dq"));
            }

            donorMap.put("lastUpdate", donor.get("lastUpdate"));


            donorsLst.add(donorMap);
        }

        String labName = "labo"+laboName;
        JSONObject freqInterne = (JSONObject) frequencies.get("interne");
        JSONObject freqLabo = (JSONObject) frequencies.get(labName);


        parameters.put("specimen", specimenList);
        parameters.put("donorsLst",  donorsLst );
        parameters.put("donorMap", donorMap);
        parameters.put("AC", freqInterne.get("AC"));
        parameters.put("AN", freqInterne.get("AC"));
        parameters.put("HC", freqInterne.get("HC"));
        parameters.put("PN", freqInterne.get("PN"));
        parameters.put("lAC", freqLabo.get("AC"));
        parameters.put("lAN", freqLabo.get("AC"));
        parameters.put("lAF", freqLabo.get("AF"));
        parameters.put("lHC", freqLabo.get("HC"));
        parameters.put("lPN", freqLabo.get("PN"));
        parameters.put("labName", labName);

        Map<String, Object> freqLaboMap = new HashMap<>();
        freqLaboMap.put("AC", freqLabo.get("AC"));
        freqLaboMap.put("AN", freqLabo.get("AN"));
        freqLaboMap.put("AF", freqLabo.get("AF"));
        freqLaboMap.put("HC", freqLabo.get("HC"));
        freqLaboMap.put("PN", freqLabo.get("PN"));
        parameters.put("freqLab", freqLaboMap);


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
                        "if (ctx._source.frequencies.interne.AN > 0) " +
                        "{ ctx._source.frequencies.interne.AF = (float) ctx._source.frequencies.interne.AC / ctx._source.frequencies.interne.AN }" +
                        "} boolean freqLaboExist = false; " +
                        "if (ctx._source.frequencies.containsKey(params.labName)) { freqLaboExist = true } " +
                        "if (toUpdateFreq && freqLaboExist) { " +
                        "ctx._source.frequencies."+labName+".AC += params.lAC; " +
                        "ctx._source.frequencies."+labName+".AN += params.lAN; " +
                        "ctx._source.frequencies."+labName+".HC += params.lHC; " +
                        "ctx._source.frequencies."+labName+".PN += params.lPN; " +
                        "if (ctx._source.frequencies."+labName+".AN > 0) " +
                        "{ ctx._source.frequencies."+labName+".AF = (float) ctx._source.frequencies."+labName+".AC / ctx._source.frequencies."+labName+".AN }" +
                        "} " +
                        "if (toUpdateFreq && !freqLaboExist) {" +
                        "ctx._source.frequencies.put(params.labName,params.freqLab)" +
                        "}"
                , parameters);


        UpdateRequest request = new UpdateRequest("mutations", "_doc", uid);
        request.script(inline);

        request.upsert(object, XContentType.JSON);
        return request;
    }


    public static String getSHA256Hash(String data) {
        String result = null;
        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            byte[] hash = digest.digest(data.getBytes(StandardCharsets.UTF_8));
            return bytesToHex(hash);
        } catch(Exception ex) {
            ex.printStackTrace();
        }
        return result;
    }

    static String getSHA1Hash(String data) {
        String result = null;

        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-1");
            byte[] hash = digest.digest(data.getBytes(StandardCharsets.UTF_8));
            return bytesToHex(hash);
        } catch(Exception ex) {
            ex.printStackTrace();
        }
        return result;
    }

    static String getMD5Hash(String data) {
        String result = null;

        try {
            MessageDigest digest = MessageDigest.getInstance("MD5");
            byte[] hash = digest.digest(data.getBytes(StandardCharsets.UTF_8));
            return bytesToHex(hash);
        } catch(Exception ex) {
            ex.printStackTrace();
        }
        return result;
    }

    private static String bytesToHex(byte[] hash) {
        return DatatypeConverter.printHexBinary(hash);
    }

    private static void connectToMemCached() {
        String[] servers = {MEMCACHED_SERVER};
        pool = SockIOPool.getInstance();
        pool.setServers( servers );
        pool.setFailover( true );
        pool.setInitConn( 10 );
        pool.setMinConn( 5 );
        pool.setMaxConn( 250 );
        pool.setMaintSleep( 30 );
        pool.setNagle( false );
        pool.setSocketTO( 3000 );
        pool.setAliveCheck( true );
        pool.initialize();
        //Get the Memcached Client from SockIOPool
        mcc = new MemCachedClient();
        System.err.println("*********** Connected to Memcached");
    }

    private static boolean checkForSpecimen(JSONArray donorArray, String specimenIdList) {


        String[] specimenIds = specimenIdList.split(",");
        for (int i = 0; i < donorArray.length(); i++) {
            JSONObject currentDonor = (JSONObject) donorArray.get(i);
            String specimen = (String) currentDonor.get("specimenId");
            for (String specimenId: specimenIds) {
                if (specimenId.equalsIgnoreCase(specimen)) {
                    return true;
                }
            }
        }
        return false;
    }

    private Properties getPropertiesFromFile(String filename) {
        Properties prop = new Properties();

        try (InputStream in =
                     getClass().getResourceAsStream(filename)) {
            prop.load(in);
        } catch (IOException e) {
            e.printStackTrace();
        }


        return prop;
    }


}