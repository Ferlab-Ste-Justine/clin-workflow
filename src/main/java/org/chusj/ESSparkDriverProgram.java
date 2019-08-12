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
import org.elasticsearch.action.index.IndexResponse;
import org.elasticsearch.client.RequestOptions;
import org.elasticsearch.client.RestClient;
import org.elasticsearch.client.RestHighLevelClient;
import org.elasticsearch.common.xcontent.XContentType;
import org.json.JSONArray;
import org.json.JSONObject;
import scala.Int;

import javax.xml.bind.DatatypeConverter;
import java.io.IOException;
import java.security.MessageDigest;

public class ESSparkDriverProgram {

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



        String[] lineValueArray = extractedLine.split("\t");
        String chrom = lineValueArray[0];

        if ("CHROM".equalsIgnoreCase(chrom)) {
            return; // Meta data line
        }

        JSONObject mutationCentricIndexjson = new JSONObject();
        JSONObject propertiesOneMutation = new JSONObject();

        String position = lineValueArray[1];
        String id = lineValueArray[2].trim();   // mdb_snp membership id
        String reference = lineValueArray[3].trim();
        String alt = lineValueArray[4].replace(",<NON_REF>", ""); // CT,<NON_REF> or G,TGG,<NON_REF>
        String af = lineValueArray[5].trim();

        String dpS = lineValueArray[6].trim();

        String ac = lineValueArray[7].trim();
        String an = lineValueArray[8].trim();
        float qual = Float.valueOf(lineValueArray[9].trim());
        //boolean db = Boolean.valueOf(lineValueArray[10]); // not used since we only use membership id and not membership true/false

        String[] alleleA = lineValueArray[11].split(",");
        String[] effetA = lineValueArray[12].split(",");
        String[] impactA = lineValueArray[13].split(",");
        String[] geneA = lineValueArray[14].split(",");
        String[] geneIdA = lineValueArray[15].split(",");
        String[] featureA = lineValueArray[16].split(",");
        String[] featureIdA = lineValueArray[17].split(",");
        String[] biotypeA = lineValueArray[18].split(",");
        String[] hgvs_cA = lineValueArray[19].split(",");
        String[] hgvs_pA = lineValueArray[20].split(",");
        String[] cdnaPosA = lineValueArray[21].split(",");
        String[] cdnaLenA = lineValueArray[22].split(",");
        String[] cdsPosA = lineValueArray[23].split(",");
        String[] cdsLenA = lineValueArray[24].split(",");
        String[] aaPosA = lineValueArray[25].split(",");
        String[] aaLenA = lineValueArray[26].split(",");
        String[] distA = lineValueArray[27].split(",");
        String[] rankA = lineValueArray[28].split(",");
        String[] errA = lineValueArray[29].split(",");

        String lof = lineValueArray[30].trim();
        String nmd = lineValueArray[31].trim();
        boolean snp = Boolean.valueOf(lineValueArray[32].trim());
        boolean mnp = Boolean.valueOf(lineValueArray[33].trim());
        boolean ins = Boolean.valueOf(lineValueArray[34].trim());
        boolean del = Boolean.valueOf(lineValueArray[35].trim());
        boolean mixed = Boolean.valueOf(lineValueArray[36].trim());
        boolean hom = Boolean.valueOf(lineValueArray[37].trim());
        boolean het = Boolean.valueOf(lineValueArray[38].trim());

        /* population frequencies */

        JSONObject frequencies = new JSONObject();
        JSONObject frequencyEx = new JSONObject();
        JSONObject frequency1000Gp3 = new JSONObject();
        JSONObject frequencyEsp6500 = new JSONObject();
        JSONObject frequencyGnomadEx = new JSONObject();
        JSONObject frequencyGnomadGen = new JSONObject();


        JSONObject predictions = new JSONObject();
        JSONObject prediction = new JSONObject();

        String exAC_AC = lineValueArray[39].trim().split(",")[0];
        String exAC_AF = lineValueArray[40].trim().split(",")[0];
        String exAC_SAS_AC = lineValueArray[41].trim().split(",")[0];
        String exAC_SAS_AF = lineValueArray[42].trim().split(",")[0];
        String exAC_ADJ_AC = lineValueArray[43].trim().split(",")[0];
        String exAC_ADJ_AF = lineValueArray[44].trim().split(",")[0];
        String exAC_AFR_AC = lineValueArray[45].trim().split(",")[0];
        String exAC_AFR_AF = lineValueArray[46].trim().split(",")[0];
        String ExAC_FIN_AC = lineValueArray[47].trim().split(",")[0];
        String ExAC_FIN_AF = lineValueArray[48].trim().split(",")[0];
        String ExAC_AMR_AC = lineValueArray[49].trim().split(",")[0];
        String ExAC_AMR_AF = lineValueArray[50].trim().split(",")[0];
        String ExAC_NFE_AC = lineValueArray[51].trim().split(",")[0];
        String ExAC_NFE_AF = lineValueArray[52].trim().split(",")[0];
        String ExAC_EAS_AC = lineValueArray[53].trim().split(",")[0];
        String ExAC_EAS_AF = lineValueArray[54].trim().split(",")[0];
        String GERP___RS = lineValueArray[55].trim().split(",")[0];
        String GERP___NR = lineValueArray[56].trim().split(",")[0];
        String d1000Gp3_AC = lineValueArray[57].trim().split(",")[0];
        String d1000Gp3_AF = lineValueArray[58].trim().split(",")[0];
        String d1000Gp3_AMR_AC = lineValueArray[59].trim().split(",")[0];
        String d1000Gp3_AMR_AF = lineValueArray[60].trim().split(",")[0];
        String d1000Gp3_EAS_AC = lineValueArray[61].trim().split(",")[0];
        String d1000Gp3_EAS_AF = lineValueArray[62].trim().split(",")[0];
        String d1000Gp3_EUR_AC = lineValueArray[63].trim().split(",")[0];
        String d1000Gp3_EUR_AF = lineValueArray[64].trim().split(",")[0];
        String d1000Gp3_AFR_AC = lineValueArray[65].trim().split(",")[0];
        String d1000Gp3_AFR_AF = lineValueArray[66].trim().split(",")[0];
        String d1000Gp3_SAS_AC = lineValueArray[67].trim().split(",")[0];
        String d1000Gp3_SAS_AF = lineValueArray[68].trim().split(",")[0];

        String FATHMM_pred = lineValueArray[69].trim().split(",")[0];
        String Polyphen2_HDIV_pred = lineValueArray[70].trim().split(",")[0];
        String Polyphen2_HVAR_pred = lineValueArray[71].trim().split(",")[0];

        String ESP6500_EA_AC = lineValueArray[72].trim().split(",")[0];
        String ESP6500_EA_AF = lineValueArray[73].trim().split(",")[0];
        String ESP6500_AA_AC = lineValueArray[74].trim().split(",")[0];
        String ESP6500_AA_AF = lineValueArray[75].trim().split(",")[0];

        String SIFT_pred = lineValueArray[76].trim().split(",")[0];

        String gnomAD_exomes_AC = lineValueArray[77].trim().split(",")[0];
        String gnomAD_exomes_AN = lineValueArray[78].trim().split(",")[0];
        String gnomAD_exomes_AF = lineValueArray[79].trim().split(",")[0];
        String gnomAD_exomes_AFR_AC = lineValueArray[80].trim().split(",")[0];
        String gnomAD_exomes_AFR_AN = lineValueArray[81].trim().split(",")[0];
        String gnomAD_exomes_AFR_AF = lineValueArray[82].trim().split(",")[0];
        String gnomAD_exomes_AMR_AC = lineValueArray[83].trim().split(",")[0];
        String gnomAD_exomes_AMR_AN = lineValueArray[84].trim().split(",")[0];
        String gnomAD_exomes_AMR_AF = lineValueArray[85].trim().split(",")[0];
        String gnomAD_exomes_ASJ_AC = lineValueArray[86].trim().split(",")[0];
        String gnomAD_exomes_ASJ_AN = lineValueArray[87].trim().split(",")[0];
        String gnomAD_exomes_ASJ_AF = lineValueArray[88].trim().split(",")[0];
        String gnomAD_exomes_EAS_AC = lineValueArray[89].trim().split(",")[0];
        String gnomAD_exomes_EAS_AN = lineValueArray[90].trim().split(",")[0];
        String gnomAD_exomes_EAS_AF = lineValueArray[91].trim().split(",")[0];
        String gnomAD_exomes_FIN_AC = lineValueArray[92].trim().split(",")[0];
        String gnomAD_exomes_FIN_AN = lineValueArray[93].trim().split(",")[0];
        String gnomAD_exomes_FIN_AF = lineValueArray[94].trim().split(",")[0];
        String gnomAD_exomes_NFE_AC = lineValueArray[95].trim().split(",")[0];
        String gnomAD_exomes_NFE_AN = lineValueArray[96].trim().split(",")[0];
        String gnomAD_exomes_NFE_AF = lineValueArray[97].trim().split(",")[0];
        String gnomAD_exomes_SAS_AC = lineValueArray[98].trim().split(",")[0];
        String gnomAD_exomes_SAS_AN = lineValueArray[99].trim().split(",")[0];
        String gnomAD_exomes_SAS_AF = lineValueArray[100].trim().split(",")[0];
        String gnomAD_genomes_AC = lineValueArray[101].trim().split(",")[0];
        String gnomAD_genomes_AN = lineValueArray[102].trim().split(",")[0];
        String gnomAD_genomes_AF = lineValueArray[103].trim().split(",")[0];
        String gnomAD_genomes_AFR_AC = lineValueArray[104].trim().split(",")[0];
        String gnomAD_genomes_AFR_AN = lineValueArray[105].trim().split(",")[0];
        String gnomAD_genomes_AFR_AF = lineValueArray[106].trim().split(",")[0];
        String gnomAD_genomes_AMR_AC = lineValueArray[107].trim().split(",")[0];
        String gnomAD_genomes_AMR_AN = lineValueArray[108].trim().split(",")[0];
        String gnomAD_genomes_AMR_AF = lineValueArray[109].trim().split(",")[0];
        String gnomAD_genomes_ASJ_AC = lineValueArray[110].trim().split(",")[0];
        String gnomAD_genomes_ASJ_AN = lineValueArray[111].trim().split(",")[0];
        String gnomAD_genomes_ASJ_AF = lineValueArray[112].trim().split(",")[0];
        String gnomAD_genomes_EAS_AC = lineValueArray[113].trim().split(",")[0];
        String gnomAD_genomes_EAS_AN = lineValueArray[114].trim().split(",")[0];
        String gnomAD_genomes_EAS_AF = lineValueArray[115].trim().split(",")[0];
        String gnomAD_genomes_FIN_AC = lineValueArray[116].trim().split(",")[0];
        String gnomAD_genomes_FIN_AN = lineValueArray[117].trim().split(",")[0];
        String gnomAD_genomes_FIN_AF = lineValueArray[118].trim().split(",")[0];
        String gnomAD_genomes_NFE_AC = lineValueArray[119].trim().split(",")[0];
        String gnomAD_genomes_NFE_AN = lineValueArray[120].trim().split(",")[0];
        String gnomAD_genomes_NFE_AF = lineValueArray[121].trim().split(",")[0];
        String clinvar_id = lineValueArray[122].trim();
        String clinvar_clnsig = lineValueArray[123].trim();
        String clinvar_trait = lineValueArray[124].trim();
        String clinvar_hgvs = lineValueArray[125].trim();
        String gq = lineValueArray[126].trim();
        String ac2 = lineValueArray[127].trim();
        String gt = lineValueArray[128].trim();
        String ad = lineValueArray[129].trim();


        boolean freqExAcAvailable = false;
        boolean freq1000Gp3Available = false;
        boolean freqEsp6500Avail = false;
        boolean freqGnomadExAvail = false;
        boolean freqGnomadGenAvail = false;
        boolean predAvail = false;

        if (exAC_AC.length()>0) {
            frequencyEx.put("AC", Integer.valueOf(exAC_AC));
            freqExAcAvailable = true;
        }
        if (exAC_AF.length()>0) {
            frequencyEx.put("AF", Float.valueOf(exAC_AF));
            freqExAcAvailable = true;
        }
        if (exAC_SAS_AC.length()>0) {
            frequencyEx.put("SAS_AC", Integer.valueOf(exAC_SAS_AC));
            freqExAcAvailable = true;
        }
        if (exAC_SAS_AF.length()>0) {
            frequencyEx.put("SAS_AF", Float.valueOf(exAC_SAS_AF));
            freqExAcAvailable = true;
        }
        if (exAC_ADJ_AC.length()>0) {
            frequencyEx.put("ADJ_AC", Integer.valueOf(exAC_ADJ_AC));
            freqExAcAvailable = true;
        }
        if (exAC_ADJ_AF.length()>0) {
            frequencyEx.put("ADJ_AF", Float.valueOf(exAC_ADJ_AF));
            freqExAcAvailable = true;
        }
        if (exAC_AFR_AC.length()>0) {
            frequencyEx.put("AFR_AC", Integer.valueOf(exAC_AFR_AC));
            freqExAcAvailable = true;
        }
        if (exAC_AFR_AF.length()>0) {
            frequencyEx.put("AFR_AF", Float.valueOf(exAC_AFR_AF));
            freqExAcAvailable = true;
        }
        if (ExAC_FIN_AC.length()>0) {
            frequencyEx.put("FIN_AC", Integer.valueOf(ExAC_FIN_AC));
            freqExAcAvailable = true;
        }
        if (ExAC_FIN_AF.length()>0) {
            frequencyEx.put("FIN_AF", Float.valueOf(ExAC_FIN_AF));
            freqExAcAvailable = true;
        }
        if (ExAC_AMR_AC.length()>0) {
            frequencyEx.put("AMR_AC", Integer.valueOf(ExAC_AMR_AC));
            freqExAcAvailable = true;
        }
        if (ExAC_AMR_AF.length()>0) {
            frequencyEx.put("AMR_AF", Float.valueOf(ExAC_AMR_AF));
            freqExAcAvailable = true;
        }
        if (ExAC_NFE_AC.length()>0) {
            frequencyEx.put("NFE_AC", Integer.valueOf(ExAC_NFE_AC));
            freqExAcAvailable = true;
        }
        if (ExAC_NFE_AF.length()>0) {
            frequencyEx.put("NFE_AF", Float.valueOf(ExAC_NFE_AF));
            freqExAcAvailable = true;
        }
        if (ExAC_EAS_AC.length()>0) {
            frequencyEx.put("EAS_AC", Integer.valueOf(ExAC_EAS_AC));
            freqExAcAvailable = true;
        }
        if (ExAC_EAS_AF.length()>0) {
            frequencyEx.put("EAS_AF", Float.valueOf(ExAC_EAS_AF));
            freqExAcAvailable = true;
        }

        if (d1000Gp3_AC.length()>0) {
            frequency1000Gp3.put("AC", Integer.valueOf(d1000Gp3_AC));
            freq1000Gp3Available = true;
        }
        if (d1000Gp3_AF.length()>0) {
            frequency1000Gp3.put("AF", Float.valueOf(d1000Gp3_AF));
            freq1000Gp3Available = true;
        }
        if (d1000Gp3_AMR_AC.length()>0) {
            frequency1000Gp3.put("AMR_AC", Integer.valueOf(d1000Gp3_AMR_AC));
            freq1000Gp3Available = true;
        }
        if (d1000Gp3_AMR_AF.length()>0) {
            frequency1000Gp3.put("AMR_AF", Float.valueOf(d1000Gp3_AMR_AF));
            freq1000Gp3Available = true;
        }
        if (d1000Gp3_EAS_AC.length()>0) {
            frequency1000Gp3.put("EAS_AC", Integer.valueOf(d1000Gp3_EAS_AC));
            freq1000Gp3Available = true;
        }
        if (d1000Gp3_EAS_AF.length()>0) {
            frequency1000Gp3.put("EAS_AF", Float.valueOf(d1000Gp3_EAS_AF));
            freq1000Gp3Available = true;
        }
        if (d1000Gp3_EUR_AC.length()>0) {
            frequency1000Gp3.put("EUR_AC", Integer.valueOf(d1000Gp3_EUR_AC));
            freq1000Gp3Available = true;
        }
        if (d1000Gp3_EUR_AF.length()>0) {
            frequency1000Gp3.put("EUR_AF", Float.valueOf(d1000Gp3_EUR_AF));
            freq1000Gp3Available = true;
        }
        if (d1000Gp3_AFR_AC.length()>0) {
            frequency1000Gp3.put("AFR_AC", Integer.valueOf(d1000Gp3_AFR_AC));
            freq1000Gp3Available = true;
        }
        if (d1000Gp3_AFR_AF.length()>0) {
            frequency1000Gp3.put("AFR_AF", Float.valueOf(d1000Gp3_AFR_AF));
            freq1000Gp3Available = true;
        }
        if (d1000Gp3_SAS_AC.length()>0) {
            frequency1000Gp3.put("SAS_AC", Integer.valueOf(d1000Gp3_SAS_AC));
            freq1000Gp3Available = true;
        }
        if (d1000Gp3_SAS_AF.length()>0) {
            frequency1000Gp3.put("SAS_AF", Float.valueOf(d1000Gp3_SAS_AF));
            freq1000Gp3Available = true;
        }


        if (gnomAD_exomes_AC.length()>0) {
            frequencyGnomadEx.put("AC", Integer.valueOf(gnomAD_exomes_AC));
            freqGnomadExAvail = true;
        }
        if (gnomAD_exomes_AN.length()>0) {
            frequencyGnomadEx.put("AN", Long.valueOf(gnomAD_exomes_AN));
            freqGnomadExAvail = true;
        }
        if (gnomAD_exomes_AF.length()>0) {
            frequencyGnomadEx.put("AF", Float.valueOf(gnomAD_exomes_AF));
            freqGnomadExAvail = true;
        }
        if (gnomAD_exomes_AFR_AC.length()>0) {
            frequencyGnomadEx.put("AFR_AC", Integer.valueOf(gnomAD_exomes_AFR_AC));
            freqGnomadExAvail = true;
        }
        if (gnomAD_exomes_AFR_AN.length()>0) {
            frequencyGnomadEx.put("AFR_AN", Long.valueOf(gnomAD_exomes_AFR_AN));
            freqGnomadExAvail = true;
        }
        if (gnomAD_exomes_AFR_AF.length()>0) {
            frequencyGnomadEx.put("AFR_AF", Float.valueOf(gnomAD_exomes_AFR_AF));
            freqGnomadExAvail = true;
        }
        if (gnomAD_exomes_AMR_AC.length()>0) {
            frequencyGnomadEx.put("AMR_AC", Integer.valueOf(gnomAD_exomes_AMR_AC));
            freqGnomadExAvail = true;
        }
        if (gnomAD_exomes_AMR_AN.length()>0) {
            frequencyGnomadEx.put("AMR_AN", Long.valueOf(gnomAD_exomes_AMR_AN));
            freqGnomadExAvail = true;
        }
        if (gnomAD_exomes_AMR_AF.length()>0) {
            frequencyGnomadEx.put("AMR_AF", Float.valueOf(gnomAD_exomes_AMR_AF));
            freqGnomadExAvail = true;
        }
        if (gnomAD_exomes_ASJ_AC.length()>0) {
            frequencyGnomadEx.put("ASJ_AC", Integer.valueOf(gnomAD_exomes_ASJ_AC));
            freqGnomadExAvail = true;
        }
        if (gnomAD_exomes_ASJ_AN.length()>0) {
            frequencyGnomadEx.put("ASJ_AN", Long.valueOf(gnomAD_exomes_ASJ_AN));
            freqGnomadExAvail = true;
        }
        if (gnomAD_exomes_ASJ_AF.length()>0) {
            frequencyGnomadEx.put("ASJ_AF", Float.valueOf(gnomAD_exomes_ASJ_AF));
            freqGnomadExAvail = true;
        }
        if (gnomAD_exomes_EAS_AC.length()>0) {
            frequencyGnomadEx.put("EAS_AC", Integer.valueOf(gnomAD_exomes_EAS_AC));
            freqGnomadExAvail = true;
        }
        if (gnomAD_exomes_EAS_AN.length()>0) {
            frequencyGnomadEx.put("EAS_AN", Long.valueOf(gnomAD_exomes_EAS_AN));
            freqGnomadExAvail = true;
        }
        if (gnomAD_exomes_EAS_AF.length()>0) {
            frequencyGnomadEx.put("EAS_AF", Float.valueOf(gnomAD_exomes_EAS_AF));
            freqGnomadExAvail = true;
        }
        if (gnomAD_exomes_FIN_AC.length()>0) {
            frequencyGnomadEx.put("FIN_AC", Integer.valueOf(gnomAD_exomes_FIN_AC));
            freqGnomadExAvail = true;
        }
        if (gnomAD_exomes_FIN_AN.length()>0) {
            frequencyGnomadEx.put("FIN_AN", Long.valueOf(gnomAD_exomes_FIN_AN));
            freqGnomadExAvail = true;
        }
        if (gnomAD_exomes_FIN_AF.length()>0) {
            frequencyGnomadEx.put("FIN_AF", Float.valueOf(gnomAD_exomes_FIN_AF));
            freqGnomadExAvail = true;
        }
        if (gnomAD_exomes_NFE_AC.length()>0) {
            frequencyGnomadEx.put("NFE_AC", Integer.valueOf(gnomAD_exomes_NFE_AC));
            freqGnomadExAvail = true;
        }
        if (gnomAD_exomes_NFE_AN.length()>0) {
            frequencyGnomadEx.put("NFE_AN", Long.valueOf(gnomAD_exomes_NFE_AN));
            freqGnomadExAvail = true;
        }
        if (gnomAD_exomes_NFE_AF.length()>0) {
            frequencyGnomadEx.put("NFE_AF", Float.valueOf(gnomAD_exomes_NFE_AF));
            freqGnomadExAvail = true;
        }
        if (gnomAD_exomes_SAS_AC.length()>0) {
            frequencyGnomadEx.put("SAS_AC", Integer.valueOf(gnomAD_exomes_SAS_AC));
            freqGnomadExAvail = true;
        }
        if (gnomAD_exomes_SAS_AN.length()>0) {
            frequencyGnomadEx.put("SAS_AN", Long.valueOf(gnomAD_exomes_SAS_AN));
            freqGnomadExAvail = true;
        }
        if (gnomAD_exomes_SAS_AF.length()>0) {
            frequencyGnomadEx.put("SAS_AF", Float.valueOf(gnomAD_exomes_SAS_AF));
            freqGnomadExAvail = true;
        }

        if (gnomAD_genomes_AC.length()>0) {
            frequencyGnomadGen.put("AC", Integer.valueOf(gnomAD_genomes_AC));
            freqGnomadGenAvail = true;
        }
        if (gnomAD_genomes_AN.length()>0) {
            frequencyGnomadGen.put("AN", Long.valueOf(gnomAD_genomes_AN));
            freqGnomadGenAvail = true;
        }
        if (gnomAD_genomes_AF.length()>0) {
            frequencyGnomadGen.put("AF", Float.valueOf(gnomAD_genomes_AF));
            freqGnomadGenAvail = true;
        }
        if (gnomAD_genomes_AFR_AC.length()>0) {
            frequencyGnomadGen.put("AFR_AC", Integer.valueOf(gnomAD_genomes_AFR_AC));
            freqGnomadGenAvail = true;
        }
        if (gnomAD_genomes_AFR_AN.length()>0) {
            frequencyGnomadGen.put("AFR_AN", Long.valueOf(gnomAD_genomes_AFR_AN));
            freqGnomadGenAvail = true;
        }
        if (gnomAD_genomes_AFR_AF.length()>0) {
            frequencyGnomadGen.put("AFR_AF", Float.valueOf(gnomAD_genomes_AFR_AF));
            freqGnomadGenAvail = true;
        }
        if (gnomAD_genomes_AMR_AC.length()>0) {
            frequencyGnomadGen.put("AMR_AC", Integer.valueOf(gnomAD_genomes_AMR_AC));
            freqGnomadGenAvail = true;
        }
        if (gnomAD_genomes_AMR_AN.length()>0) {
            frequencyGnomadGen.put("AMR_AN", Long.valueOf(gnomAD_genomes_AMR_AN));
            freqGnomadGenAvail = true;
        }
        if (gnomAD_genomes_AMR_AF.length()>0) {
            frequencyGnomadGen.put("AMR_AF", Float.valueOf(gnomAD_genomes_AMR_AF));
            freqGnomadGenAvail = true;
        }
        if (gnomAD_genomes_ASJ_AC.length()>0) {
            frequencyGnomadGen.put("ASJ_AC", Integer.valueOf(gnomAD_genomes_ASJ_AC));
            freqGnomadGenAvail = true;
        }
        if (gnomAD_genomes_ASJ_AN.length()>0) {
            frequencyGnomadGen.put("ASJ_AN", Long.valueOf(gnomAD_genomes_ASJ_AN));
            freqGnomadGenAvail = true;
        }
        if (gnomAD_genomes_ASJ_AF.length()>0) {
            frequencyGnomadGen.put("ASJ_AF", Float.valueOf(gnomAD_genomes_ASJ_AF));
            freqGnomadGenAvail = true;
        }
        if (gnomAD_genomes_EAS_AC.length()>0) {
            frequencyGnomadGen.put("EAS_AC", Integer.valueOf(gnomAD_genomes_EAS_AC));
            freqGnomadGenAvail = true;
        }
        if (gnomAD_genomes_EAS_AN.length()>0) {
            frequencyGnomadGen.put("EAS_AN", Long.valueOf(gnomAD_genomes_EAS_AN));
            freqGnomadGenAvail = true;
        }
        if (gnomAD_genomes_EAS_AF.length()>0) {
            frequencyGnomadGen.put("EAS_AF", Float.valueOf(gnomAD_genomes_EAS_AF));
            freqGnomadGenAvail = true;
        }
        if (gnomAD_genomes_FIN_AC.length()>0) {
            frequencyGnomadGen.put("FIN_AC", Integer.valueOf(gnomAD_genomes_FIN_AC));
            freqGnomadGenAvail = true;
        }
        if (gnomAD_genomes_FIN_AN.length()>0) {
            frequencyGnomadGen.put("FIN_AN", Long.valueOf(gnomAD_genomes_FIN_AN));
            freqGnomadGenAvail = true;
        }
        if (gnomAD_genomes_FIN_AF.length()>0) {
            frequencyGnomadGen.put("FIN_AF", Float.valueOf(gnomAD_genomes_FIN_AF));
            freqGnomadGenAvail = true;
        }
        if (gnomAD_genomes_NFE_AC.length()>0) {
            frequencyGnomadGen.put("NFE_AC", Integer.valueOf(gnomAD_genomes_NFE_AC));
            freqGnomadGenAvail = true;
        }
        if (gnomAD_genomes_NFE_AN.length()>0) {
            frequencyGnomadGen.put("NFE_AN", Long.valueOf(gnomAD_genomes_NFE_AN));
            freqGnomadGenAvail = true;
        }
        if (gnomAD_genomes_NFE_AF.length()>0) {
            frequencyGnomadGen.put("NFE_AF", Float.valueOf(gnomAD_genomes_NFE_AF));
            freqGnomadGenAvail = true;
        }

        if (ESP6500_EA_AC.length()>0) {
            frequencyEsp6500.put("EA_AC", Integer.valueOf(ESP6500_EA_AC));
            freqEsp6500Avail = true;
        }
        if (ESP6500_EA_AF.length()>0) {
            frequencyEsp6500.put("EA_AF", Float.valueOf(ESP6500_EA_AF));
            freqEsp6500Avail = true;
        }
        if (ESP6500_AA_AC.length()>0) {
            frequencyEsp6500.put("AA_AC", Integer.valueOf(ESP6500_AA_AC));
            freqEsp6500Avail = true;
        }
        if (ESP6500_AA_AF.length()>0) {
            frequencyEsp6500.put("AA_AF", Float.valueOf(ESP6500_AA_AF));
            freqEsp6500Avail = true;
        }

        if (FATHMM_pred.length()>0) {
            prediction.put("FATHMM", FATHMM_pred);
            predAvail = true;
        }
        if (Polyphen2_HDIV_pred.length()>0) {
            prediction.put("Polyphen2_HDIV", Polyphen2_HDIV_pred);
            predAvail = true;
        }
        if (Polyphen2_HVAR_pred.length()>0) {
            prediction.put("Polyphen2_HVAR", Polyphen2_HVAR_pred);
            predAvail = true;
        }
        if (SIFT_pred.length()>0) {
            prediction.put("SIFT", SIFT_pred);
            predAvail = true;
        }

        if (clinvar_id.length()>0) propertiesOneMutation.put("clinvar_id", clinvar_id);
        if (clinvar_clnsig.length()>0) propertiesOneMutation.put("clinvar_clnsig", clinvar_clnsig);
        if (clinvar_trait.length()>0) propertiesOneMutation.put("clinvar_trait", clinvar_trait);
        if (clinvar_hgvs.length()>0) propertiesOneMutation.put("clinvar_hgvs", clinvar_hgvs);
        if (gq.length()>0) propertiesOneMutation.put("gq", gq);
        if (ac2.length()>0) propertiesOneMutation.put("ac2", gq);
        if (gt.length()>0) propertiesOneMutation.put("gt", gt);
        if (ad.length()>0) propertiesOneMutation.put("ad", ad);

        frequencyEx.put("available", freqExAcAvailable);
        frequency1000Gp3.put("available", freq1000Gp3Available);
        frequencyEsp6500.put("available", freqEsp6500Avail);
        frequencyGnomadEx.put("available", freqGnomadExAvail);
        frequencyGnomadGen.put("available", freqGnomadGenAvail);

        prediction.put("available", predAvail);


        frequencies.put("ExAc", frequencyEx);
        frequencies.put("1000Gp3", frequency1000Gp3);
        frequencies.put("ESP6500", frequencyEsp6500);
        frequencies.put("gnomAD_exomes", frequencyGnomadEx);
        frequencies.put("gnomAD_genomes", frequencyGnomadGen);

        propertiesOneMutation.put("frequencies", frequencies);
        propertiesOneMutation.put("Predictions", prediction);


        if (GERP___RS.length()>0) {
            propertiesOneMutation.put("GERP___RS", Float.valueOf(GERP___RS));
            //freqExAcAvailable = true;
        }
        if (GERP___NR.length()>0) {
            propertiesOneMutation.put("GERP___NR", Float.valueOf(GERP___NR));
            //freqExAcAvailable = true;
        }

        String chrPos = chrom.substring(3); // remove 'chr'
        String mutation = reference + ">" + alt.split(",")[0];
        String dnaChanges = chrPos + ":g." + position + mutation;

        String uid = getSHA256Hash(dnaChanges);

        propertiesOneMutation.put("id", uid);
        propertiesOneMutation.put("mutationId", dnaChanges);
        propertiesOneMutation.put("mutation", mutation);
        propertiesOneMutation.put("chrom",chrPos);

        //Pattern p = Pattern.compile("-?\\d+");
        //Matcher m;
//        m = p.matcher(chrom);
//        if (m.find()) {
//            propertiesOneMutation.put("chrom",Integer.valueOf(m.group()));
//        }

        int pos = Integer.valueOf(position);
        propertiesOneMutation.put("start",pos);

        String[] altArray = alt.split(",");

        int altLen =0;
        for (String a : altArray) {
            altLen += a.length();
        }

        //propertiesOneMutation.put("end",(pos+altLen-1));
        propertiesOneMutation.put("reference", reference);
        propertiesOneMutation.put("alt", alt);


        // af is currently unused
        if (af.length() > 0 ) propertiesOneMutation.put("af", af);

        // ac is currently unused
        if (ac.length() > 0 ) propertiesOneMutation.put("ac", ac);
        // an is currently unused
        if (an.length() > 0 ) propertiesOneMutation.put("an", an);

        if (id.length()>0) propertiesOneMutation.put("dbSNP_ID",id);
        //propertiesOneMutation.put("dbSNP_membership", db);

        if (lof.length()>0) propertiesOneMutation.put("lossOfFunctionEffect", lof);
        if (nmd.length()>0) propertiesOneMutation.put("nonsenseMediatedDecayEffects", nmd);
        if (snp) propertiesOneMutation.put("type", "SNP");
        if (mnp) propertiesOneMutation.put("type", "MNP");
        if (ins) propertiesOneMutation.put("type", "INS");
        if (del) propertiesOneMutation.put("type", "DEL");
        if (mixed) propertiesOneMutation.put("type", "MIXED");


        JSONArray functionalAnnotations = new JSONArray();
        int qtyFunctionalAnnotation = alleleA.length;

        for (int i=0; i<qtyFunctionalAnnotation; i++) {
            int valueInt;
            JSONObject functionalAnnotation = new JSONObject();
            functionalAnnotation.put("allele", alleleA[i]);
            functionalAnnotation.put("effet", effetA[i]);
            functionalAnnotation.put("impact", impactA[i]);
            String gene = geneA[i].trim();
            if (gene.length()>0) functionalAnnotation.put("gene", gene);
            String geneid = geneIdA[i].trim();
            if (geneid.length()>0) functionalAnnotation.put("geneId", geneid);

            functionalAnnotation.put("featureType", featureA[i]);

            functionalAnnotation.put("featureId", featureIdA[i]);
            String biotype = biotypeA[i].trim();
            if (biotype.length()>0) functionalAnnotation.put("biotype", biotype);
            if (rankA[i].trim().length()>0) {
                valueInt = Integer.valueOf(rankA[i]);
                if(valueInt>0) functionalAnnotation.put("rank/total",valueInt);
            }

            String hgvs_c = hgvs_cA[i].trim();
            if (hgvs_c.length()>0) functionalAnnotation.put("hgvsC", hgvs_c);
            String hgvs_p = hgvs_pA[i].trim();
            if (hgvs_p.length()>0) functionalAnnotation.put("hgvsP", hgvs_p);
            if (cdnaPosA[i].trim().length()>0) {
                valueInt = Integer.valueOf(cdnaPosA[i]);
                if (valueInt > 0) {
                    functionalAnnotation.put("cdnaPos", valueInt);
                    functionalAnnotation.put("cdnaLen", Integer.valueOf(cdnaLenA[i]));
                }
            }
            if (cdsPosA[i].trim().length()>0) {
                valueInt = Integer.valueOf(cdsPosA[i]);
                if (valueInt >= 0) {
                    functionalAnnotation.put("cdsPos", valueInt);
                    functionalAnnotation.put("cdsLen", Integer.valueOf(cdsLenA[i]));
                }
            }
            if (aaPosA[i].trim().length()>0) {
                valueInt = Integer.valueOf(aaPosA[i]);
                if (valueInt >= 0) {
                    functionalAnnotation.put("aaPos", valueInt);
                    functionalAnnotation.put("aaLen", Integer.valueOf(aaLenA[i]));
                }
            }
            if (distA[i].trim().length()>0) functionalAnnotation.put("distance", Integer.valueOf(distA[i]));
            String err = errA[i].trim();
            if (err.length()>0) functionalAnnotation.put("errorsWarningsInfo", err);
            functionalAnnotations.put(functionalAnnotation);
        }
        propertiesOneMutation.put("functionalAnnotations", functionalAnnotations);
//        _summary = new JSONObject();
//        _summary.put("_affected_donor_count", 1);
//        _summary.put("_affected_project_count", 1);
//        _summary.append("_affected_project_id", projectIdO);
//        _summary.put("_tested_donor_count", 1);
//        propertiesOneMutation.put("_summary",_summary);
        propertiesOneMutation.put("assemblyVersion", "GRCh38");
        propertiesOneMutation.put("annotationTool", "snpEff/snpSift 4.3t");

        JSONObject newDonor = new JSONObject();
        newDonor.put("familyId", familyId);
        newDonor.put("donorId", patientId);
        newDonor.put("studyId", studyId);
        //newDonor.put("isProband", true);
        newDonor.put("type", type);
        if (dpS.length()>0) newDonor.put("depth", Integer.valueOf(dpS));
        newDonor.put("quality", qual);
        newDonor.put("sequencingStrategy", sequencingStrategy);
        newDonor.put("pid", projectId);
        if (hom) newDonor.put("zygosity", "Hom");
        if (het) newDonor.put("zygosity", "Het");

        JSONArray donorArray = null;//= new JSONArray();

        // verify if variant already exist
        boolean toIndex = false;
        boolean toMemcached = false;
        boolean checkES = false;
        String msg= "";

        // try in memcached first
        String donorArrayStr = (String) mcc.get(uid);
        if (donorArrayStr != null) {
            donorArray = new JSONArray(donorArrayStr);
            if (donorArray != null && donorArray.length() > 0) {

                boolean donorFound = checkForDonor(donorArray, patientId);

                if (!donorFound) {
                    // add new donor to previous one
                    msg += "m0";
                    donorArray.put(newDonor);
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
            // might need to check
            toMemcached = true;
            checkES = true;
        }

        if (checkES) {
            GetRequest getRequest = new GetRequest("variants", "family", uid);

            if (client.exists(getRequest, RequestOptions.DEFAULT)) {

                GetResponse getResponse = client.get(getRequest, RequestOptions.DEFAULT);

                if (getResponse.isExists()) {

                    JSONObject obj = new JSONObject(getResponse.getSourceAsString());

                    JSONObject props = (JSONObject) obj.get("properties");

                    String mutationId = (String) props.get("mutationId");
                    //msg += mutationId;

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

                    boolean donorFound = checkForDonor(donorArray, patientId);

                    if (!donorFound) {
                        // add new donor to previous one
                        msg += "e0";
                        //System.out.print("0)");
                        donorArray.put(newDonor);
                        toIndex = true;

                    } else {
                        msg += "e1";
                        donorArray.put(newDonor);
                        // nothing
                    }
                    System.out.print(msg);
                } else {
                    donorArray = new JSONArray();
                    donorArray.put(newDonor);
                    toIndex = true;
                }
            } else {
                toIndex = true;
                donorArray = new JSONArray();
                donorArray.put(newDonor);
            }
        }

        if (toIndex) {
            //int retry = 0;
            //if (donorArray==null) donorArray = new JSONArray();
            boolean success = false;
            propertiesOneMutation.put("donor", donorArray);
            mutationCentricIndexjson.put("properties", propertiesOneMutation);
            for (int i=0; i< 3; i++) {
                try {
                    index(mutationCentricIndexjson.toString(), client, uid, "variants");
                    success = true;
                    break;
                } catch (Exception e) {
                    System.err.println("*********** Try #"+i+" failed...");
                    continue;
                }
            }
            if (!success) System.err.println("#########\n\n\n######### Unable to index " + uid + "\n############");

        }
        if (toMemcached) {
            mcc.set(uid,
                    donorArray.toString());
        }

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

}