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
import java.util.*;

import static org.chusj.ESSparkDriverProgram.checkForDonor;
import static org.chusj.ESSparkDriverProgram.getSHA256Hash;
import static org.chusj.PatientHelper.loadPedigree;
import static org.chusj.VEPSparkDriverProgram.getMD5Hash;
import static org.chusj.VepHelper.getPropertiesFromFile;


public class Validator {

    public static RestHighLevelClient client;
    private static SockIOPool pool;
    private static MemCachedClient mcc;
    private static long FOUND = 0L;
    private static long NOT_FOUND = 0L;
    private static long PROCESS = 0L;
    private static long DISCARD = 0L;
    private static long ERROR = 0L;


    public static void main(String args[]) throws Exception {

        if (args.length != 4) {
            args = new String[]{"out10.txt","pedigree.properties", "pedigree.ped", "validate"};
        }

        String extractFile = args[0];
        String pedigrePropsFile = args[1];
        String pedFile = args[2];

        Properties pedigreeProps = getPropertiesFromFile(pedigrePropsFile);
        List<Pedigree> pedigrees = loadPedigree(pedFile);

        Map<String, Patient> patientMap = PatientHelper.preparePedigreeFromProps(pedigreeProps);
        pedigrees.forEach(System.out::println);
        patientMap.forEach((k,v)->System.out.println(k+"\n\t"+v));
        String build = pedigreeProps.getProperty("assemblyVersion");

        try (BufferedReader buf = new BufferedReader(new FileReader(extractFile));
             RestHighLevelClient client = new RestHighLevelClient(
                     RestClient.builder(
                             new HttpHost("localhost", 9200, "http"))) ) {

            String fetchedLine;
            //connectToMemCached();

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

                    PROCESS++;
                    String[] lineValueArray = fetchedLine.split("\t");
                    String chrom = lineValueArray[0];

                    String position = lineValueArray[1];
                    String reference = lineValueArray[2].trim();
                    String alt = lineValueArray[3].replace(",<NON_REF>", ""); // CT,<NON_REF> or G,TGG,<NON_REF>
                    String chrPos = chrom.substring(3); // remove 'chr'
                    if (chrPos.length() > 2) {
                        DISCARD++;
                        continue;
                    }
                    String mutation = reference + ">" + alt.split(",")[0];
                    String dnaChanges = "chr" + chrPos + ":g." + position + mutation;


//                    String uid = getSHA256Hash(dnaChanges);
                    String uid = getMD5Hash(dnaChanges +"@"+build);

                    if (setTest.contains(uid)) {
                        // found in memory!!!
                        System.err.println("*****  Found " + dnaChanges + " in memory - uid ="+uid);
                    } else {
                        setTest.add(uid);
                    }

                    // check in ES
                    String msg;
                    GetRequest getRequest = new GetRequest("mutations", "_doc", uid);

                    if (client.exists(getRequest, RequestOptions.DEFAULT)) {
                        GetResponse getResponse = client.get(getRequest, RequestOptions.DEFAULT);

                        if (getResponse.isExists()) {

                            JSONObject obj = new JSONObject(getResponse.getSourceAsString());
                            JSONArray donorArray = new JSONArray();

                            String mutationId = (String) obj.get("mutationId");

                            if (!mutationId.equalsIgnoreCase(dnaChanges)) {
                                // collision detected
                                msg = "(" + dnaChanges + ";" + uid + ";";
                                System.err.println("collision="+msg);
                            }
                            System.out.print(".");
                            FOUND++;

//                            boolean donorFound = ESSparkDriverProgram.checkForDonor(donorArray, patientId);
//                            if (!donorFound) {
//                                System.err.println("*****  Donor not Found in ES:" + dnaChanges + " - uid ="+uid);
//                            } else {
//                                System.out.print("e1");
//                            }
                        } else {
                            System.err.println("*****  Cannot get record from ES for " + dnaChanges + " - uid = "+uid);
                            ERROR++;

                        }

                    } else {
                        System.err.println("#### Record not Found in ES:" + dnaChanges + " - uid = "+uid);
                        NOT_FOUND++;
                    }
                }
            }
        }
        long total = (FOUND+NOT_FOUND+DISCARD+ERROR);
        System.out.println("");
        System.out.println("Number of Record Processed = "+ PROCESS);
        System.out.println("Number of Record Found     = "+ FOUND);
        System.out.println("Number of Record Not Found = "+ NOT_FOUND);
        System.out.println("Number of Record Discarded = "+ DISCARD);
        System.out.println("Number of error Record     = "+ ERROR);
        System.out.println("Total                      = "+ total);
        System.out.println(( NOT_FOUND == 0L ) ? "passed" : "failed" );

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
