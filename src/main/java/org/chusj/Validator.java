package org.chusj;


import com.whalin.MemCached.MemCachedClient;
import com.whalin.MemCached.SockIOPool;
import org.apache.http.HttpHost;
import org.elasticsearch.action.get.GetRequest;
import org.elasticsearch.action.get.GetResponse;
import org.elasticsearch.client.RequestOptions;
import org.elasticsearch.client.RestClient;
import org.elasticsearch.client.RestHighLevelClient;
import org.json.JSONArray;
import org.json.JSONObject;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashSet;
import java.util.Set;

import static org.chusj.ESSparkDriverProgram.checkForDonor;
import static org.chusj.ESSparkDriverProgram.getSHA256Hash;


public class Validator {

    public static RestHighLevelClient client;
    public static SockIOPool pool;
    public static MemCachedClient mcc;


    public static void main(String args[]) throws Exception {

        String extractFile = args[0];
        String patientId = args[1];

        try (BufferedReader buf = new BufferedReader(new FileReader(extractFile));
             RestHighLevelClient client = new RestHighLevelClient(
                     RestClient.builder(
                             new HttpHost("localhost", 9200, "http"))) ) {

            String fetchedLine;
            connectToMemCached();

            //JSONObject projectIdO = new JSONObject();
            //projectIdO.put("pid", projectId);

            Set<String> setTest = new HashSet<>();


            // Record and prep metaData
            fetchedLine = buf.readLine();

            // Main

            while (true) {
                fetchedLine = buf.readLine();
                if (fetchedLine == null) {
                    break;
                } else {

                    System.out.print(".");
                    String[] lineValueArray = fetchedLine.split("\t");
                    String chrom = lineValueArray[0];

                    String position = lineValueArray[1];
                    String reference = lineValueArray[3].trim();
                    String alt = lineValueArray[4].replace(",<NON_REF>", ""); // CT,<NON_REF> or G,TGG,<NON_REF>
                    String chrPos = chrom.substring(3); // remove 'chr'
                    String mutation = reference + ">" + alt.split(",")[0];
                    String dnaChanges = chrPos + ":g." + position + mutation;

                    String uid = getSHA256Hash(dnaChanges);

                    if (setTest.contains(uid)) {
                        // found in memory!!!
                        System.err.println("*****  Found " + dnaChanges + " in memory - uid ="+uid);
                    } else {
                        setTest.add(uid);
                    }

                    // Check in Memcached
                    String donorArrayStr = (String) mcc.get(uid);
                    if (donorArrayStr != null) {
                        JSONArray donorArray = new JSONArray(donorArrayStr);

                        if (donorArray != null && donorArray.length() > 0) {
                            boolean donorFound = checkForDonor(donorArray, patientId);

                            if (!donorFound) {
                                System.err.println("*****  Donor not Found in Memcached:" + dnaChanges + " - uid ="+uid);
                            } else {
                                System.out.print("m1");
                            }

                        } else {
                            System.err.println("*****  No donor array Found in Memcached:" + dnaChanges + " - uid ="+uid);
                        }

                    } else {
                        System.err.println("*****  No Record not Found in Memcached:" + dnaChanges + " - uid ="+uid);
                    }

                    // check in ES
                    String msg;
                    GetRequest getRequest = new GetRequest("variants", "family", uid);

                    if (client.exists(getRequest, RequestOptions.DEFAULT)) {
                        GetResponse getResponse = client.get(getRequest, RequestOptions.DEFAULT);

                        if (getResponse.isExists()) {

                            JSONObject obj = new JSONObject(getResponse.getSourceAsString());
                            JSONObject props = (JSONObject) obj.get("properties");
                            JSONArray donorArray = new JSONArray();

                            String mutationId = (String) props.get("mutationId");

                            if (!mutationId.equalsIgnoreCase(dnaChanges)) {
                                // collision detected
                                msg = "(" + dnaChanges + ";" + uid + ";";
                                System.err.println("collision="+msg);
                            }
                            try {
                                donorArray = (JSONArray) props.get("donor");
                            } catch (Exception e) {
                                donorArray = new JSONArray();
                            }

                            boolean donorFound = ESSparkDriverProgram.checkForDonor(donorArray, patientId);
                            if (!donorFound) {
                                System.err.println("*****  Donor not Found in ES:" + dnaChanges + " - uid ="+uid);
                            } else {
                                System.out.print("e1");
                            }
                        } else {
                            System.err.println("*****  No Record not Found in ES:" + dnaChanges + " - uid ="+uid);

                        }

                    } else {
                        System.err.println("#### not exist()");
                    }
                }
            }
        }
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


}
