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
import org.elasticsearch.client.RequestOptions;
import org.elasticsearch.client.RestClient;
import org.elasticsearch.client.RestHighLevelClient;
import org.elasticsearch.common.xcontent.XContentType;
import org.json.JSONArray;
import org.json.JSONObject;

import javax.xml.bind.DatatypeConverter;
import java.io.IOException;
import java.io.InputStream;
import java.security.MessageDigest;
import java.util.Properties;

import static org.chusj.VepHelper.addDonorArrayToDonorArray;

public class VEPSparkDriverProgram {

    public static RestHighLevelClient client;
    public static SockIOPool pool;
    public static MemCachedClient mcc;
    private static String MEMCACHED_SERVER = "localhost:11212";
    private static String INDEX_NAME = "mutations";


    public static void main(String args[]) throws Exception {

        if (args.length == 0 || args.length != 12 ) {
            throw new Exception("Missing params; need extractFile patientId familyId projectId studyId sequencingStrategy type(P|M|F) ES_LOAD(Y|N) " +
                    "sparkMaster local[nbWorkers] 8g nbPartitions");
        }

        String extractFile = args[0];
        String patientId = args[1];
        String familyId = args[2];
        String projectId = args[3];
        String studyId = args[4];
        String sequencingStrategy = args[5];
        String type = args[6];
        String esLoad = args[7];
        String sparkMaster = args[8];
        String localThread = args[9];
        String memory = args[10];
        int nbPartitions = Integer.valueOf(args[11]);
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

        Properties pedigreeProps = VepHelper.getPropertiesFromFile("pedigree.properties");

        try (RestHighLevelClient clientTry = new RestHighLevelClient(
                RestClient.builder(
                        new HttpHost("localhost", 9200, "http")))) {

            client = clientTry;
            connectToMemCached();

            lines.foreach((a) -> buildAndStoreJsonObj(a, esLoad, pedigreeProps));
        }
        sc.close();
        client.close();
    }


    public static void buildAndStoreJsonObj(String extractedLine, String esLoad, Properties pedigreeProps) throws IOException {


        JSONObject propertiesOneMutation = VepHelper.processVcfDataLine(extractedLine, "dn,dq", pedigreeProps);

        // extract donor info if not found
        //String qual = (String) propertiesOneMutation.remove("qual");
        //String filter = (String) propertiesOneMutation.remove("filter");
        if (propertiesOneMutation == null) return;
        JSONArray donorArray = (JSONArray) propertiesOneMutation.remove("donors");

        JSONArray newDonorArray;
        // verify if variant already exist
        boolean toMapping = false;
        boolean toIndex = false;
        boolean toMemcached = false;
        boolean checkES = false;
        String msg= "";



        String uid = (String) propertiesOneMutation.get("id");
        String dnaChanges = (String) propertiesOneMutation.get("mutationId");

        String specimenIdList = pedigreeProps.getProperty("pedigree");

        JSONArray previousDonorArray = null;
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
            //int retry = 0;
            //if (donorArray==null) donorArray = new JSONArray();

            propertiesOneMutation.put("donors", previousDonorArray);
            //mutationCentricIndexjson.put("properties", propertiesOneMutation);
            for (int i=0; i< 3; i++) {
                try {
                    index(propertiesOneMutation.toString(), client, uid, INDEX_NAME);
                    indexingSuccess = true;
                    break;
                } catch (Exception e) {
                    System.err.println("*********** Try #"+i+" failed...");
                    continue;
                }
            }

            if (!indexingSuccess) System.err.println("#########\n\n\n######### Unable to index " + uid + "\n############");

        }
        if (toMemcached && indexingSuccess) {
            mcc.set(uid,
                    previousDonorArray.toString());
        }

    }


    public static void index(String object, RestHighLevelClient client, String uid, String index) throws Exception {
        IndexRequest request = new IndexRequest("mutations", "_doc", uid);

        request.source(object, XContentType.JSON);
        //IndexResponse indexResponse;
        try {
            client.index(request, RequestOptions.DEFAULT);
            //System.out.print("."+indexResponse.status().getStatus()+".");
        } catch (Exception e) {
            System.err.println("E="+e.getMessage()+"-"+uid);
            //indexResponse = client.index(request, RequestOptions.DEFAULT);
            throw new Exception(e);
        }

        System.out.print(".");
    }

    public static String getSHA256Hash(String data) {
        String result = null;
        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            byte[] hash = digest.digest(data.getBytes("UTF-8"));
            return bytesToHex(hash);
        } catch(Exception ex) {
            ex.printStackTrace();
        }
        return result;
    }

    public static String getSHA1Hash(String data) {
        String result = null;

        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-1");
            byte[] hash = digest.digest(data.getBytes("UTF-8"));
            return bytesToHex(hash);
        } catch(Exception ex) {
            ex.printStackTrace();
        }
        return result;
    }

    public static String getMD5Hash(String data) {
        String result = null;

        try {
            MessageDigest digest = MessageDigest.getInstance("MD5");
            byte[] hash = digest.digest(data.getBytes("UTF-8"));
            return bytesToHex(hash);
        } catch(Exception ex) {
            ex.printStackTrace();
        }
        return result;
    }

    public static String bytesToHex(byte[] hash) {
        return DatatypeConverter.printHexBinary(hash);
    }

    public static void connectToMemCached() {
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

    public static boolean checkForSpecimen(JSONArray donorArray, String specimenIdList) {


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