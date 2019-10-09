package org.chusj;

import com.whalin.MemCached.MemCachedClient;
import com.whalin.MemCached.SockIOPool;
import org.apache.http.HttpHost;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.elasticsearch.action.get.GetRequest;
import org.elasticsearch.action.get.GetResponse;
import org.elasticsearch.action.index.IndexRequest;
import org.elasticsearch.action.update.UpdateRequest;
import org.elasticsearch.action.update.UpdateResponse;
import org.elasticsearch.client.RequestOptions;
import org.elasticsearch.client.RestClient;
import org.elasticsearch.client.RestHighLevelClient;
import org.elasticsearch.common.xcontent.XContentBuilder;
import org.elasticsearch.common.xcontent.XContentFactory;
import org.elasticsearch.common.xcontent.XContentType;
import org.elasticsearch.index.engine.Engine;
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

import static java.util.Collections.singletonMap;
import static org.chusj.VepHelper.addDonorArrayToDonorArray;

public class VEPSparkDriverProgram {

    public static RestHighLevelClient client;
    private static SockIOPool pool;
    private static MemCachedClient mcc;
    private static String MEMCACHED_SERVER = "localhost:11212";
    private static String INDEX_NAME = "mutations";


    public static void main(String[] args) throws Exception {

        if (args.length != 7 ) {
            throw new Exception("Missing params; need extractFile pedigreePropertiesFile ES_UPSERT (true|false) " +
                    "sparkMaster local[nbWorkers] 8g nbPartitions");
        }

        String extractFile = args[0];
        String pedigreePropsFile = args[1];

        boolean esUpsert = Boolean.valueOf(args[2]);
        String sparkMaster = args[3];
        String localThread = args[4];
        String memory = args[5];
        int nbPartitions = Integer.valueOf(args[6]);
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

        //System.out.println("Total lines in log file " + lines.count());

        /* Map operation -> Mapping number of characters into each line as RDD */
        //JavaRDD<Integer> lineCharacters = lines.map(s -> s.length());
        /* reduces operation -> Calculating total characters */

        Properties pedigreeProps = VepHelper.getPropertiesFromFile(pedigreePropsFile);
        pedigreeProps.forEach( (k,v) -> System.out.println(k+"="+v) );

        try (RestHighLevelClient clientTry = new RestHighLevelClient(
                RestClient.builder(
                        new HttpHost("localhost", 9200, "http")))) {

            client = clientTry;
            connectToMemCached();

            lines.foreach((a) -> buildAndStoreJsonObj(a, esUpsert, pedigreeProps));
        }
        sc.close();
        client.close();
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



        String uid = (String) propertiesOneMutation.get("id");
        String dnaChanges = (String) propertiesOneMutation.get("mutationId");

        String specimenIdList = pedigreeProps.getProperty("pedigree");
        String familyId = pedigreeProps.getProperty("familyId");

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

            previousDonorArray = (JSONArray) propertiesOneMutation.get("donors");
            specimentList = (JSONArray) propertiesOneMutation.get("specimenList");
            for (int i=0; i< 3; i++) {
                try {
                    if (esUpsert) {
                        upsert(propertiesOneMutation.toString(), client, uid, previousDonorArray, specimentList);
                    } else {
                        index(propertiesOneMutation.toString(), client, uid, INDEX_NAME);
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


    private static void index(String object, RestHighLevelClient client, String uid, String index) throws Exception {
        IndexRequest request = new IndexRequest("mutations", "_doc", uid);

        request.source(object, XContentType.JSON);
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

    private static void upsert(String object, RestHighLevelClient client, String uid, JSONArray donors, JSONArray specimenList) throws Exception {



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
            //donorMap.put("pid", donor.get("pid"));
            donorMap.put("relation", donor.get("relation"));
            donorMap.put("sequencingStrategy", donor.get("sequencingStrategy"));
            donorMap.put("studyId", donor.get("studyId"));
            donorMap.put("zygosity", donor.get("zygosity"));
            donorMap.put("ad", donor.get("ad"));
            donorMap.put("af", donor.get("af"));
            donorMap.put("dp", donor.get("dp"));
            donorMap.put("gt", donor.get("gt"));
            donorMap.put("gq", donor.get("gq"));
            donorMap.put("practitionerId", donor.get("practitionerId"));
            donorMap.put("organizationId", donor.get("organizationId"));
            if (!donor.isNull("genotypeFamily")) {
                donorMap.put("genotypeFamily", donor.get("genotypeFamily"));
            }

            donorMap.put("lastUpdate", donor.get("lastUpdate"));


            donorsLst.add(donorMap);
//            System.err.println("donorsLst@"+i+"="+donorsLst.get(i));
//            System.err.println("spIds@"+i+"="+donorMap.get("specimenId"));

        }

        parameters.put("specimen", specimenList);
        parameters.put("donorsLst",  donorsLst );
        parameters.put("donorMap", donorMap);

        Script inline = new Script(ScriptType.INLINE, "painless",
                "for (int i=0; i<params.specimen.size(); i++) {" +
                            "String s = params.specimen.get(i); " +
                            "Map d = params.donorsLst.get(i); " +
                            "if (!ctx._source.specimenList.contains(s) && d != null) {" +
                               " ctx._source.donors.add(d);" +
                               " ctx._source.specimenList.add(s)" +
                            "}" +
                          "}", parameters);

        UpdateRequest request = new UpdateRequest("mutations", "_doc", uid);
        request.script(inline);

        request.upsert(object, XContentType.JSON);
        UpdateResponse result = null;
        try {
            result = client.update(request, RequestOptions.DEFAULT);
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