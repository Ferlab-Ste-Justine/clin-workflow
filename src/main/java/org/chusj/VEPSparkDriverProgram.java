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

public class VEPSparkDriverProgram {

    public static RestHighLevelClient client;
    public static SockIOPool pool;
    public static MemCachedClient mcc;


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

        try (RestHighLevelClient clientTry = new RestHighLevelClient(
                RestClient.builder(
                        new HttpHost("localhost", 9200, "http")))) {

            client = clientTry;
            connectToMemCached();

            lines.foreach((a) -> buildAndStoreJsonObj(a, familyId, patientId, studyId, type, projectId, sequencingStrategy, esLoad));
        }
        sc.close();
        client.close();
    }


    public static void buildAndStoreJsonObj(String extractedLine, String familyId, String patientId, String studyId,
                                            String type, String projectId, String sequencingStrategy, String esLoad) throws IOException {


        String[] pedigree = {"14140,P","14141,M", "14142,F"};
        JSONObject propertiesOneMutation = VepHelper.processVcfDataLine(extractedLine, "dn,dq", pedigree);
        propertiesOneMutation.put("assemblyVersion", "GRCh38");
        propertiesOneMutation.put("annotationTool", "VEP");
        propertiesOneMutation.put("annotationToolVersion", 96.3);

        // extract donor info if not found
        //String qual = (String) propertiesOneMutation.remove("qual");
        //String filter = (String) propertiesOneMutation.remove("filter");





    }

    public static void index(String object, RestHighLevelClient client, String uid, String index) throws Exception {
        IndexRequest request = new IndexRequest("variants", "family", uid);

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

    public static String bytesToHex(byte[] hash) {
        return DatatypeConverter.printHexBinary(hash);
    }

    public static void connectToMemCached() {
        String[] servers = {"localhost:11211"};
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

    public static boolean checkForDonor(JSONArray donorArray, String patientId) {

        for (int i = 0; i < donorArray.length(); i++) {
            JSONObject currentDonor = (JSONObject) donorArray.get(i);
            String donorId = (String) currentDonor.get("donorId");
            if (patientId.equalsIgnoreCase(donorId)) {
                return true;
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