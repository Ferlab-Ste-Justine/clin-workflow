package org.chusj;

import org.apache.http.HttpHost;
import org.elasticsearch.client.RestClient;
import org.elasticsearch.client.RestHighLevelClient;
import org.json.JSONObject;
import org.json.simple.JSONArray;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

public class RedisGeneSetIndexer {

    public static RestHighLevelClient client;
    public static int bulkOpsQty = 51;

    public static void main(String[] args) throws Exception {
        if (args.length != 2) {
            args = new String[]{"Homo_sapiens.gene_info.txt", "9200"};
        }

        String extractFile = args[0];
        int esPort = Integer.valueOf(args[1]);

        try (RestHighLevelClient clientTry = new RestHighLevelClient(
                RestClient.builder(
                        new HttpHost("localhost", esPort, "http")))) {

            client = clientTry;
            VEPSparkDriverProgram.client = clientTry;

            try (BufferedReader buf = new BufferedReader(new FileReader(extractFile))) {

                String fetchedLine;

                buf.readLine();
                List<JSONObject> jsonObjectList = new ArrayList<>();

                while (true) {
                    fetchedLine = buf.readLine();
                    if (fetchedLine == null) {
                        break;
                    } else {
                        //System.out.println(fetchedLine);
                        String[] words = fetchedLine.split("\t");
                        //val ensemblId = words(5).trim.split("[|]", 0).last.replace("Ensembl:", "")
                        String[] ensArray = words[5].trim().split("[|]", 0);
                        for (String ens: ensArray) {
                            String ensemblId = ens.replace("Ensembl:", "");
                            if (ensemblId.startsWith("ENS")) {
                                //System.out.println("id:"+ensemblId);
                                JSONObject gene = new JSONObject();
                                // private static void addGeneSetsToObjs(String ensId, JSONObject jsonObject, JSONObject availObj,
                                //                                          Map<String, Patient> patientMap, List<Pedigree> pedigrees)
                                gene.put("ensemblId", ensemblId);
                                gene.put("variants", new JSONArray());
                                gene.put("donors", new JSONArray());
                                gene.put("frequencies", new JSONArray());
                                VepHelper.addGeneSetsToObjs(ensemblId, gene, new JSONObject(), null, new ArrayList<>());
                                if (gene.isNull("geneId")) {
                                    System.err.println(" Empty gene set for id:"+ensemblId);
                                } else {
                                    gene.put("id", gene.get("ensemblId"));
                                    jsonObjectList.add(gene);
                                }
                                if (gene.isNull("alias")) {
                                    gene.put("alias", new JSONArray());
                                }
                            }
                        }
                        if (jsonObjectList.size() >= bulkOpsQty) {
                            //System.out.println("Bulk Items in partition-" + jsonObjectList.size());
                            VEPSparkDriverProgram.bulkStoreJsonObj(
                                    jsonObjectList, false,false, true);
                            jsonObjectList = new ArrayList<>();
                        }
                    }
                }
                VEPSparkDriverProgram.bulkStoreJsonObj(
                        jsonObjectList, false, false, true);
            }
            client.close();
        }
    }

}
