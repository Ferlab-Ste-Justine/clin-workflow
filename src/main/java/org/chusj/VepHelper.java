package org.chusj;


import org.json.JSONArray;
import org.json.JSONObject;

import java.io.BufferedReader;
import java.io.FileReader;
import java.text.NumberFormat;
import java.text.ParseException;

import static org.chusj.VEPSparkDriverProgram.getSHA256Hash;

public class VepHelper {

    private static float avgFuncAnnoPerMutation = 0.0f;
    private static int countFuncAnnoPerMutation = 0;
    private static int countMutation =0;
    private static NumberFormat NF = NumberFormat.getInstance();


    public static void main(String[] args) throws Exception {

        String extractFile = args[0];
        String nbPatient = args[1];

        try (BufferedReader buf = new BufferedReader(new FileReader(extractFile))) {

            String fetchedLine;


            // Record and prep metaData
            buf.readLine();

            // Main

            while (true) {
                fetchedLine = buf.readLine();
                if (fetchedLine == null) {
                    break;
                } else {
                    processVcfDataLine(fetchedLine);
                }
            }


        }
        if (countMutation>0) {
            avgFuncAnnoPerMutation = (float) countFuncAnnoPerMutation / countMutation;
        }
        System.out.println("\navgFuncAnnoPerMutation="+avgFuncAnnoPerMutation);

    }


    private static JSONObject processVcfDataLine(String extractedLine) {

        //CHROM	POS	ID	REF	ALT	QUAL	FILTER	DP	MQ	MQRankSum
        // ReadPosRankSum	LOD	FractionInformativeReads	SNP	MNP	INS	DEL	MIXED	HOM	HET
        // GEN[*].GT	GEN[*].AD	GEN[*].AF	GEN[*].F1R2	GEN[*].F2R1	GEN[*].DP	GEN[*].SB	GEN[*].MB	CSQ
        //System.out.print(".");
        String[] lineValueArray = extractedLine.split("\t");
        // System.out.println("lineValueArray.length="+lineValueArray.length);
        // dynamic positioning system -- pos counter++
        int pos = 0;

        JSONObject propertiesOneMutation = new JSONObject();

        String chrom = lineValueArray[pos++];
        if ("CHROM".equalsIgnoreCase(chrom)) {
            return null; // Meta data line
        }
        countMutation++;
        String position = lineValueArray[pos++];
        String id = lineValueArray[pos++].trim();   // mdb_snp membership id
        String reference = lineValueArray[pos++].trim();
        String alt = lineValueArray[pos++].replace(",<NON_REF>", ""); // CT,<NON_REF> or G,TGG,<NON_REF>
        String qual = lineValueArray[pos++].trim();
        String filter = lineValueArray[pos++].trim();
        String dpS = lineValueArray[pos++].trim();
        String mqS = lineValueArray[pos++].trim();
        String mqRankSum = lineValueArray[pos++].trim();
        String readPosRankSum = lineValueArray[pos++].trim();
        String lod = lineValueArray[pos++].trim();
        String fractionInformativeReads = lineValueArray[pos++].trim();

        boolean snp = Boolean.valueOf(lineValueArray[pos++].trim());
        boolean mnp = Boolean.valueOf(lineValueArray[pos++].trim());
        boolean ins = Boolean.valueOf(lineValueArray[pos++].trim());
        boolean del = Boolean.valueOf(lineValueArray[pos++].trim());
        boolean mixed = Boolean.valueOf(lineValueArray[pos++].trim());
        boolean hom = Boolean.valueOf(lineValueArray[pos++].trim());
        boolean het = Boolean.valueOf(lineValueArray[pos++].trim());


        if (snp) propertiesOneMutation.put("type", "SNP");
        if (mnp) propertiesOneMutation.put("type", "MNP");
        if (ins) propertiesOneMutation.put("type", "INS");
        if (del) propertiesOneMutation.put("type", "DEL");
        if (mixed) propertiesOneMutation.put("type", "MIXED");
        propertiesOneMutation.put("reference", reference);
        propertiesOneMutation.put("alt", alt);
        if (id.length()>0) propertiesOneMutation.put("dbSNP_ID",id);

        String gt = lineValueArray[pos++].trim();
        String ad = lineValueArray[pos++].trim();
        String af = lineValueArray[pos++].trim();
        String f1r2 = lineValueArray[pos++].trim();
        String f2r1 = lineValueArray[pos++].trim();
        String genDP = lineValueArray[pos++].trim();
        String sb = lineValueArray[pos++].trim();
        String mb = lineValueArray[pos++].trim();
        String csq = lineValueArray[pos++].trim();
        String[] csqArray = csq.split(",");  // return functionalAnnotation array
        countFuncAnnoPerMutation += csqArray.length;

        String chrPos = chrom.substring(3); // remove 'chr'
        String mutation = reference + ">" + alt.split(",")[0];
        String dnaChanges = chrPos + ":g." + position + mutation;
        String uid = getSHA256Hash(dnaChanges);

        propertiesOneMutation.put("id", uid);
        propertiesOneMutation.put("mutationId", dnaChanges);
        propertiesOneMutation.put("mutation", mutation);
        propertiesOneMutation.put("chrom",chrPos);

        propertiesOneMutation.put("assemblyVersion", "GRCh38");
        propertiesOneMutation.put("annotationTool", "snpEff/snpSift 4.3t");

//        JSONObject funcAnnotation = null;

        JSONArray functionalAnnotations = new JSONArray();

        for (String s : csqArray) {
            functionalAnnotations.put(processVepAnnotations(s, dnaChanges));
        }

        propertiesOneMutation.put("functionalAnnotations", functionalAnnotations);

        return propertiesOneMutation;

    }

    private static JSONObject processVepAnnotations(String csqLine, String id) {

        //System.out.print("\n"+csqLine);
        String[] functionalAnnotationArray = csqLine.split("[|]", -1);
        // dynamic positioning system -- pos counter++

//        for (String s : functionalAnnotationArray) {
//            String val = s.trim();
//
//            //System.out.print(i+"'"+val+"'");
//        }

        //System.out.print(id+"("+functionalAnnotationArray.length+")");

        int pos = 0;

        if (functionalAnnotationArray.length < 403) {
            System.out.println(" " + csqLine + " short");
        }

        JSONObject funcAnnotation = new JSONObject();
        JSONObject funcAnnoProperties = new JSONObject();
        JSONObject frequencies = new JSONObject();
        JSONObject frequencyExAc = new JSONObject();
        JSONObject frequency1000Gp3 = new JSONObject();
        JSONObject frequencyEsp6500 = new JSONObject();
        JSONObject frequencyGnomadEx = new JSONObject();
        JSONObject frequencyGnomadGen = new JSONObject();

        JSONObject predictions = new JSONObject();
        JSONObject prediction = new JSONObject();



        //0 -
        addStrToJsonObject("allele", functionalAnnotationArray[pos++], funcAnnotation, false);
//        String Consequence = functionalAnnotationArray[pos++].trim();
        addStrToJsonObject("consequence", functionalAnnotationArray[pos++], funcAnnotation, false);
        addStrToJsonObject("impact", functionalAnnotationArray[pos++], funcAnnotation, false);
        addStrToJsonObject("gene", functionalAnnotationArray[pos++], funcAnnotation, false);
        addStrToJsonObject("geneId", functionalAnnotationArray[pos++], funcAnnotation, false);
        addStrToJsonObject("featureType", functionalAnnotationArray[pos++], funcAnnotation, false);
        addStrToJsonObject("featureId", functionalAnnotationArray[pos++], funcAnnotation, false);
        addStrToJsonObject("biotype", functionalAnnotationArray[pos++], funcAnnotation, false);
        //String EXON = functionalAnnotationArray[pos++].trim();
        addStrToJsonObject("exon", functionalAnnotationArray[pos++], funcAnnotation, false);
        //String INTRON = functionalAnnotationArray[pos++].trim();
        addStrToJsonObject("intron", functionalAnnotationArray[pos++], funcAnnotation, false);
        // 9 -
        addStrToJsonObject("hgvsC", functionalAnnotationArray[pos++], funcAnnotation, false);
        addStrToJsonObject("hgvsP", functionalAnnotationArray[pos++], funcAnnotation, false);
        addNumberToJsonObject("cdnaPos", functionalAnnotationArray[pos++] , funcAnnotation, false);
        addNumberToJsonObject("cdsPos", functionalAnnotationArray[pos++] , funcAnnotation, false);
        //String Protein_position = functionalAnnotationArray[pos++].trim();
        addNumberToJsonObject("ProteinPos", functionalAnnotationArray[pos++] , funcAnnotation, false);
        String Amino_acids = functionalAnnotationArray[pos++].trim();
        //addNumberToJsonObject("aaLen", functionalAnnotationArray[pos++] , funcAnnotation, false);
//        String Codons = functionalAnnotationArray[pos++].trim();
        addStrToJsonObject("codons", functionalAnnotationArray[pos++], funcAnnotation, false);
        String Existing_variation = functionalAnnotationArray[pos++].trim();
        addNumberToJsonObject("distance", functionalAnnotationArray[pos++] , funcAnnotation, false);
        addNumberToJsonObject("strand", functionalAnnotationArray[pos++] , funcAnnotation, false);
        //19 -
        String FLAGS = functionalAnnotationArray[pos++].trim();
        String SYMBOL_SOURCE = functionalAnnotationArray[pos++].trim();
        String HGNC_ID = functionalAnnotationArray[pos++].trim();
        String CANONICAL = functionalAnnotationArray[pos++].trim();
        String GIVEN_REF = functionalAnnotationArray[pos++].trim();
        String USED_REF = functionalAnnotationArray[pos++].trim();
        String BAM_EDIT = functionalAnnotationArray[pos++].trim();

        addNumberToJsonObject("AC", functionalAnnotationArray[pos++] , frequency1000Gp3, true);
        addNumberToJsonObject("AF", functionalAnnotationArray[pos++] , frequency1000Gp3, true);
        addNumberToJsonObject("AFR_AC", functionalAnnotationArray[pos++] , frequency1000Gp3, true);
        // 29
        addNumberToJsonObject("AFR_AF", functionalAnnotationArray[pos++] , frequency1000Gp3, true);
        addNumberToJsonObject("AMR_AC", functionalAnnotationArray[pos++] , frequency1000Gp3, true);
        addNumberToJsonObject("AMR_AF", functionalAnnotationArray[pos++] , frequency1000Gp3, true);
        addNumberToJsonObject("EAS_AC", functionalAnnotationArray[pos++] , frequency1000Gp3, true);
        addNumberToJsonObject("EAS_AF", functionalAnnotationArray[pos++] , frequency1000Gp3, true);
        addNumberToJsonObject("EUR_AC", functionalAnnotationArray[pos++] , frequency1000Gp3, true);
        addNumberToJsonObject("EUR_AF", functionalAnnotationArray[pos++] , frequency1000Gp3, true);
        addNumberToJsonObject("SAS_AC", functionalAnnotationArray[pos++] , frequency1000Gp3, true);
        addNumberToJsonObject("SAS_AF", functionalAnnotationArray[pos++] , frequency1000Gp3, true);
        String ALSPAC_AC = functionalAnnotationArray[pos++].trim();
        //39 -
        String ALSPAC_AF = functionalAnnotationArray[pos++].trim();
        String APPRIS = functionalAnnotationArray[pos++].trim();
        String Aloft_Confidence = functionalAnnotationArray[pos++].trim();
        String Aloft_Fraction_transcripts_affected = functionalAnnotationArray[pos++].trim();
        String Aloft_pred = functionalAnnotationArray[pos++].trim();
        String Aloft_prob_Dominant = functionalAnnotationArray[pos++].trim();
        String Aloft_prob_Recessive = functionalAnnotationArray[pos++].trim();
        String Aloft_prob_Tolerant = functionalAnnotationArray[pos++].trim();
        String AltaiNeandertal = functionalAnnotationArray[pos++].trim();
        String Ancestral_allele = functionalAnnotationArray[pos++].trim();
        //49 -
        String CADD_phred = functionalAnnotationArray[pos++].trim();
        String CADD_raw = functionalAnnotationArray[pos++].trim();
        String CADD_raw_rankscore = functionalAnnotationArray[pos++].trim();
        String DANN_rankscore = functionalAnnotationArray[pos++].trim();
        String DANN_score = functionalAnnotationArray[pos++].trim();
        String DEOGEN2_pred = functionalAnnotationArray[pos++].trim();
        String DEOGEN2_rankscore = functionalAnnotationArray[pos++].trim();
        String DEOGEN2_score = functionalAnnotationArray[pos++].trim();
        String Denisova = functionalAnnotationArray[pos++].trim();
        addNumberToJsonObject("AA_AC", functionalAnnotationArray[pos++] , frequencyEsp6500, true);
        //59 -
        addNumberToJsonObject("AA_AF", functionalAnnotationArray[pos++] , frequencyEsp6500, true);
        addNumberToJsonObject("EA_AC", functionalAnnotationArray[pos++] , frequencyEsp6500, true);
        addNumberToJsonObject("EA_AF", functionalAnnotationArray[pos++] , frequencyEsp6500, true);
        String Eigen_PC_phred_coding = functionalAnnotationArray[pos++].trim();
        String Eigen_PC_raw_coding = functionalAnnotationArray[pos++].trim();
        String Eigen_PC_raw_coding_rankscore = functionalAnnotationArray[pos++].trim();
        String Eigen_pred_coding = functionalAnnotationArray[pos++].trim();
        String Eigen_raw_coding = functionalAnnotationArray[pos++].trim();
        String Eigen_raw_coding_rankscore = functionalAnnotationArray[pos++].trim();
        String Ensembl_geneid = functionalAnnotationArray[pos++].trim();
        // 69 -
        String Ensembl_proteinid = functionalAnnotationArray[pos++].trim();
        String Ensembl_transcriptid = functionalAnnotationArray[pos++].trim();
        addNumberToJsonObject("AC", functionalAnnotationArray[pos++] , frequencyExAc, true);
        addNumberToJsonObject("AF", functionalAnnotationArray[pos++] , frequencyExAc, true);
        addNumberToJsonObject("AFR_AC", functionalAnnotationArray[pos++] , frequencyExAc, true);
        addNumberToJsonObject("AFR_AF", functionalAnnotationArray[pos++] , frequencyExAc, true);
        addNumberToJsonObject("AMR_AC", functionalAnnotationArray[pos++] , frequencyExAc, true);
        addNumberToJsonObject("AMR_AF", functionalAnnotationArray[pos++] , frequencyExAc, true);
        addNumberToJsonObject("Adj_AC", functionalAnnotationArray[pos++] , frequencyExAc, true);
        addNumberToJsonObject("Adj_AF", functionalAnnotationArray[pos++] , frequencyExAc, true);
        // 79
        addNumberToJsonObject("EAS_AC", functionalAnnotationArray[pos++] , frequencyExAc, true);
        addNumberToJsonObject("EAS_AF", functionalAnnotationArray[pos++] , frequencyExAc, true);
        addNumberToJsonObject("FIN_AC", functionalAnnotationArray[pos++] , frequencyExAc, true);
        addNumberToJsonObject("FIN_AF", functionalAnnotationArray[pos++] , frequencyExAc, true);
        addNumberToJsonObject("NFE_AC", functionalAnnotationArray[pos++] , frequencyExAc, true);
        addNumberToJsonObject("NFE_AF", functionalAnnotationArray[pos++] , frequencyExAc, true);
        addNumberToJsonObject("SAS_AC", functionalAnnotationArray[pos++] , frequencyExAc, true);
        addNumberToJsonObject("SAS_AF", functionalAnnotationArray[pos++] , frequencyExAc, true);
        String ExAC_nonTCGA_AC = functionalAnnotationArray[pos++].trim();
        String ExAC_nonTCGA_AF = functionalAnnotationArray[pos++].trim();
        // 89
        String ExAC_nonTCGA_AFR_AC = functionalAnnotationArray[pos++].trim();
        String ExAC_nonTCGA_AFR_AF = functionalAnnotationArray[pos++].trim();
        String ExAC_nonTCGA_AMR_AC = functionalAnnotationArray[pos++].trim();
        String ExAC_nonTCGA_AMR_AF = functionalAnnotationArray[pos++].trim();
        String ExAC_nonTCGA_Adj_AC = functionalAnnotationArray[pos++].trim();
        String ExAC_nonTCGA_Adj_AF = functionalAnnotationArray[pos++].trim();
        String ExAC_nonTCGA_EAS_AC = functionalAnnotationArray[pos++].trim();
        String ExAC_nonTCGA_EAS_AF = functionalAnnotationArray[pos++].trim();
        String ExAC_nonTCGA_FIN_AC = functionalAnnotationArray[pos++].trim();
        String ExAC_nonTCGA_FIN_AF = functionalAnnotationArray[pos++].trim();
        // 99
        String ExAC_nonTCGA_NFE_AC = functionalAnnotationArray[pos++].trim();
        String ExAC_nonTCGA_NFE_AF = functionalAnnotationArray[pos++].trim();
        String ExAC_nonTCGA_SAS_AC = functionalAnnotationArray[pos++].trim();
        String ExAC_nonTCGA_SAS_AF = functionalAnnotationArray[pos++].trim();
        String ExAC_nonpsych_AC = functionalAnnotationArray[pos++].trim();
        String ExAC_nonpsych_AF = functionalAnnotationArray[pos++].trim();
        String ExAC_nonpsych_AFR_AC = functionalAnnotationArray[pos++].trim();
        String ExAC_nonpsych_AFR_AF = functionalAnnotationArray[pos++].trim();
        String ExAC_nonpsych_AMR_AC = functionalAnnotationArray[pos++].trim();
        String ExAC_nonpsych_AMR_AF = functionalAnnotationArray[pos++].trim();
        // 109
        String ExAC_nonpsych_Adj_AC = functionalAnnotationArray[pos++].trim();
        String ExAC_nonpsych_Adj_AF = functionalAnnotationArray[pos++].trim();
        String ExAC_nonpsych_EAS_AC = functionalAnnotationArray[pos++].trim();
        String ExAC_nonpsych_EAS_AF = functionalAnnotationArray[pos++].trim();
        String ExAC_nonpsych_FIN_AC = functionalAnnotationArray[pos++].trim();
        String ExAC_nonpsych_FIN_AF = functionalAnnotationArray[pos++].trim();
        String ExAC_nonpsych_NFE_AC = functionalAnnotationArray[pos++].trim();
        String ExAC_nonpsych_NFE_AF = functionalAnnotationArray[pos++].trim();
        String ExAC_nonpsych_SAS_AC = functionalAnnotationArray[pos++].trim();
        String ExAC_nonpsych_SAS_AF = functionalAnnotationArray[pos++].trim();
        // 119
        String FATHMM_converted_rankscore = functionalAnnotationArray[pos++].trim();
        addStrToJsonObject("FATHMM", functionalAnnotationArray[pos++].trim(), prediction, true);
        addStrToJsonObject("FATHMM_score", functionalAnnotationArray[pos++], prediction, true);
        String GENCODE_basic = functionalAnnotationArray[pos++].trim();
        String GERPpp_NR = functionalAnnotationArray[pos++].trim();
        String GERPpp_RS = functionalAnnotationArray[pos++].trim();
        String GERPpp_RS_rankscore = functionalAnnotationArray[pos++].trim();
        String GM12878_confidence_value = functionalAnnotationArray[pos++].trim();
        String GM12878_fitCons_rankscore = functionalAnnotationArray[pos++].trim();
        String GM12878_fitCons_score = functionalAnnotationArray[pos++].trim();
        // 129
        String GTEx_V7_gene = functionalAnnotationArray[pos++].trim();
        String GTEx_V7_tissue = functionalAnnotationArray[pos++].trim();
        String GenoCanyon_rankscore = functionalAnnotationArray[pos++].trim();
        String GenoCanyon_score = functionalAnnotationArray[pos++].trim();
        String Geuvadis_eQTL_target_gene = functionalAnnotationArray[pos++].trim();
        String H1_hESC_confidence_value = functionalAnnotationArray[pos++].trim();
        String H1_hESC_fitCons_rankscore = functionalAnnotationArray[pos++].trim();
        String H1_hESC_fitCons_score = functionalAnnotationArray[pos++].trim();
        String HGVSc_ANNOVAR = functionalAnnotationArray[pos++].trim();
        String HGVSc_VEP = functionalAnnotationArray[pos++].trim();
        // 139
        String HGVSc_snpEff = functionalAnnotationArray[pos++].trim();
        String HGVSp_ANNOVAR = functionalAnnotationArray[pos++].trim();
        String HGVSp_VEP = functionalAnnotationArray[pos++].trim();
        String HGVSp_snpEff = functionalAnnotationArray[pos++].trim();
        String HUVEC_confidence_value = functionalAnnotationArray[pos++].trim();
        String HUVEC_fitCons_rankscore = functionalAnnotationArray[pos++].trim();
        String HUVEC_fitCons_score = functionalAnnotationArray[pos++].trim();
        String Interpro_domain = functionalAnnotationArray[pos++].trim();
        String LINSIGHT = functionalAnnotationArray[pos++].trim();
        String LINSIGHT_rankscore = functionalAnnotationArray[pos++].trim();
        // 149
        String LRT_Omega = functionalAnnotationArray[pos++].trim();
        String LRT_converted_rankscore = functionalAnnotationArray[pos++].trim();
        String LRT_pred = functionalAnnotationArray[pos++].trim();
        String LRT_score = functionalAnnotationArray[pos++].trim();
        String M_CAP_pred = functionalAnnotationArray[pos++].trim();
        String M_CAP_rankscore = functionalAnnotationArray[pos++].trim();
        String M_CAP_score = functionalAnnotationArray[pos++].trim();
        String MPC_rankscore = functionalAnnotationArray[pos++].trim();
        String MPC_score = functionalAnnotationArray[pos++].trim();
        String MVP_rankscore = functionalAnnotationArray[pos++].trim();
        // 159
        String MVP_score = functionalAnnotationArray[pos++].trim();
        String MetaLR_pred = functionalAnnotationArray[pos++].trim();
        String MetaLR_rankscore = functionalAnnotationArray[pos++].trim();
        String MetaLR_score = functionalAnnotationArray[pos++].trim();
        String MetaSVM_pred = functionalAnnotationArray[pos++].trim();
        String MetaSVM_rankscore = functionalAnnotationArray[pos++].trim();
        String MetaSVM_score = functionalAnnotationArray[pos++].trim();
        String MutPred_AAchange = functionalAnnotationArray[pos++].trim();
        String MutPred_Top5features = functionalAnnotationArray[pos++].trim();
        String MutPred_protID = functionalAnnotationArray[pos++].trim();
        // 169
        String MutPred_rankscore = functionalAnnotationArray[pos++].trim();
        String MutPred_score = functionalAnnotationArray[pos++].trim();
        String MutationAssessor_pred = functionalAnnotationArray[pos++].trim();
        String MutationAssessor_rankscore = functionalAnnotationArray[pos++].trim();
        String MutationAssessor_score = functionalAnnotationArray[pos++].trim();
        String MutationTaster_AAE = functionalAnnotationArray[pos++].trim();
        String MutationTaster_converted_rankscore = functionalAnnotationArray[pos++].trim();
        String MutationTaster_model = functionalAnnotationArray[pos++].trim();
        String MutationTaster_pred = functionalAnnotationArray[pos++].trim();
        String MutationTaster_score = functionalAnnotationArray[pos++].trim();
        // 179
        String PROVEAN_converted_rankscore = functionalAnnotationArray[pos++].trim();
        String PROVEAN_pred = functionalAnnotationArray[pos++].trim();
        String PROVEAN_score = functionalAnnotationArray[pos++].trim();
        addStrToJsonObject("Polyphen2_HDIV", functionalAnnotationArray[pos++].trim(), prediction, true);
        String Polyphen2_HDIV_rankscore = functionalAnnotationArray[pos++].trim();
        addStrToJsonObject("Polyphen2_HDIV_score", functionalAnnotationArray[pos++], prediction, true);
        addStrToJsonObject("Polyphen2_HVAR_pred", functionalAnnotationArray[pos++].trim(), prediction, true);
        String Polyphen2_HVAR_rankscore = functionalAnnotationArray[pos++].trim();
        addStrToJsonObject("Polyphen2_HVAR_score", functionalAnnotationArray[pos++], prediction, true);
        String PrimateAI_pred = functionalAnnotationArray[pos++].trim();
        // 189
        String PrimateAI_rankscore = functionalAnnotationArray[pos++].trim();
        String PrimateAI_score = functionalAnnotationArray[pos++].trim();
        String REVEL_rankscore = functionalAnnotationArray[pos++].trim();
        String REVEL_score = functionalAnnotationArray[pos++].trim();
        String Reliability_index = functionalAnnotationArray[pos++].trim();
        String SIFT4G_converted_rankscore = functionalAnnotationArray[pos++].trim();
        String SIFT4G_pred = functionalAnnotationArray[pos++].trim();
        String SIFT4G_score = functionalAnnotationArray[pos++].trim();
        String SIFT_converted_rankscore = functionalAnnotationArray[pos++].trim();
        addStrToJsonObject("SIFT", functionalAnnotationArray[pos++], prediction, true);
        // 199
        addStrToJsonObject("SIFT_score", functionalAnnotationArray[pos++].trim(), prediction, true);
        String SiPhy_29way_logOdds = functionalAnnotationArray[pos++].trim();
        String SiPhy_29way_logOdds_rankscore = functionalAnnotationArray[pos++].trim();
        String SiPhy_29way_pi = functionalAnnotationArray[pos++].trim();
        String TSL = functionalAnnotationArray[pos++].trim();
        String TWINSUK_AC = functionalAnnotationArray[pos++].trim();
        String TWINSUK_AF = functionalAnnotationArray[pos++].trim();
        String UK10K_AC = functionalAnnotationArray[pos++].trim();
        String UK10K_AF = functionalAnnotationArray[pos++].trim();
        String Uniprot_acc = functionalAnnotationArray[pos++].trim();
        // 209
        String Uniprot_entry = functionalAnnotationArray[pos++].trim();
        String VEP_canonical = functionalAnnotationArray[pos++].trim();
        String VEST4_rankscore = functionalAnnotationArray[pos++].trim();
        String VEST4_score = functionalAnnotationArray[pos++].trim();
        String VindijiaNeandertal = functionalAnnotationArray[pos++].trim();
//        String aaalt = functionalAnnotationArray[pos++].trim();
        addStrToJsonObject("aaAlt", functionalAnnotationArray[pos++], funcAnnotation, false);
        //String aapos = functionalAnnotationArray[pos++].trim();
        addNumberToJsonObject("aaPos", functionalAnnotationArray[pos++] , funcAnnotation, false);
//        String aaref = functionalAnnotationArray[pos++].trim();
        addStrToJsonObject("aaRef", functionalAnnotationArray[pos++], funcAnnotation, false);
//        String alt = functionalAnnotationArray[pos++].trim();
        addStrToJsonObject("alt", functionalAnnotationArray[pos++], funcAnnotation, false);
        String bStatistic = functionalAnnotationArray[pos++].trim();
        // 219
        String bStatistic_rankscore = functionalAnnotationArray[pos++].trim();
        String cds_strand = functionalAnnotationArray[pos++].trim();
        String chr = functionalAnnotationArray[pos++].trim();
        String clinvar_MedGen_id = functionalAnnotationArray[pos++].trim();
        addStrToJsonObject("clinvar_OMIM_id", functionalAnnotationArray[pos++], funcAnnotation, false);
        addStrToJsonObject("clinvar_Orphanet_id", functionalAnnotationArray[pos++], funcAnnotation, false);
        String clinvar_clnsig = functionalAnnotationArray[pos++].trim();
        String clinvar_hgvs = functionalAnnotationArray[pos++].trim();
        String clinvar_id = functionalAnnotationArray[pos++].trim();
        String clinvar_review = functionalAnnotationArray[pos++].trim();
        // 229
        String clinvar_trait = functionalAnnotationArray[pos++].trim();
        String clinvar_var_source = functionalAnnotationArray[pos++].trim();
        String codon_degeneracy = functionalAnnotationArray[pos++].trim();
        String codonpos = functionalAnnotationArray[pos++].trim();
        String fathmm_MKL_coding_group = functionalAnnotationArray[pos++].trim();
        String fathmm_MKL_coding_pred = functionalAnnotationArray[pos++].trim();
        String fathmm_MKL_coding_rankscore = functionalAnnotationArray[pos++].trim();
        String fathmm_MKL_coding_score = functionalAnnotationArray[pos++].trim();
        String fathmm_XF_coding_pred = functionalAnnotationArray[pos++].trim();
        String fathmm_XF_coding_rankscore = functionalAnnotationArray[pos++].trim();
        // 239
        String fathmm_XF_coding_score = functionalAnnotationArray[pos++].trim();
//        String genename = functionalAnnotationArray[pos++].trim();
        addStrToJsonObject("geneName", functionalAnnotationArray[pos++], funcAnnotation, false);

        addNumberToJsonObject("AC", functionalAnnotationArray[pos++] , frequencyGnomadEx, true);
        addNumberToJsonObject("AF", functionalAnnotationArray[pos++] , frequencyGnomadEx, true);
        addNumberToJsonObject("AFR_AC", functionalAnnotationArray[pos++] , frequencyGnomadEx, true);
        addNumberToJsonObject("AFR_AF", functionalAnnotationArray[pos++] , frequencyGnomadEx, true);

        String gnomAD_exomes_AFR_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_AFR_nhomalt = functionalAnnotationArray[pos++].trim();

        addNumberToJsonObject("AMR_AC", functionalAnnotationArray[pos++] , frequencyGnomadEx, true);
        addNumberToJsonObject("AMR_AF", functionalAnnotationArray[pos++] , frequencyGnomadEx, true);
        // 249
        String gnomAD_exomes_AMR_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_AMR_nhomalt = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_AN = functionalAnnotationArray[pos++].trim();

        addNumberToJsonObject("ASJ_AC", functionalAnnotationArray[pos++] , frequencyGnomadEx, true);
        addNumberToJsonObject("ASJ_AF", functionalAnnotationArray[pos++] , frequencyGnomadEx, true);
        String gnomAD_exomes_ASJ_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_ASJ_nhomalt = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_EAS_AC = functionalAnnotationArray[pos++].trim();

        addNumberToJsonObject("EAS_AF", functionalAnnotationArray[pos++] , frequencyGnomadEx, true);
        addNumberToJsonObject("EAS_AN", functionalAnnotationArray[pos++] , frequencyGnomadEx, true);
        // 259
        String gnomAD_exomes_EAS_nhomalt = functionalAnnotationArray[pos++].trim();

        addNumberToJsonObject("FIN_AC", functionalAnnotationArray[pos++] , frequencyGnomadEx, true);
        addNumberToJsonObject("FIN_AF", functionalAnnotationArray[pos++] , frequencyGnomadEx, true);
        String gnomAD_exomes_FIN_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_FIN_nhomalt = functionalAnnotationArray[pos++].trim();
        addNumberToJsonObject("NFE_AC", functionalAnnotationArray[pos++] , frequencyGnomadEx, true);
        addNumberToJsonObject("NFE_AF", functionalAnnotationArray[pos++] , frequencyGnomadEx, true);
        String gnomAD_exomes_NFE_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_NFE_nhomalt = functionalAnnotationArray[pos++].trim();
        // 269
        addNumberToJsonObject("POPMAX_AC", functionalAnnotationArray[pos++] , frequencyGnomadEx, true);
        // 269
        addNumberToJsonObject("POPMAX_AF", functionalAnnotationArray[pos++] , frequencyGnomadEx, true);
        String gnomAD_exomes_POPMAX_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_POPMAX_nhomalt = functionalAnnotationArray[pos++].trim();
        addNumberToJsonObject("SAS_AC", functionalAnnotationArray[pos++] , frequencyGnomadEx, true);
        addNumberToJsonObject("SAS_AF", functionalAnnotationArray[pos++] , frequencyGnomadEx, true);
        String gnomAD_exomes_SAS_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_SAS_nhomalt = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_AC = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_AF = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_AFR_AC = functionalAnnotationArray[pos++].trim();
        // 279
        String gnomAD_exomes_controls_AFR_AF = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_AFR_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_AFR_nhomalt = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_AMR_AC = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_AMR_AF = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_AMR_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_AMR_nhomalt = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_ASJ_AC = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_ASJ_AF = functionalAnnotationArray[pos++].trim();
        // 289
        String gnomAD_exomes_controls_ASJ_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_ASJ_nhomalt = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_EAS_AC = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_EAS_AF = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_EAS_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_EAS_nhomalt = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_FIN_AC = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_FIN_AF = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_FIN_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_FIN_nhomalt = functionalAnnotationArray[pos++].trim();
        // 299
        String gnomAD_exomes_controls_NFE_AC = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_NFE_AF = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_NFE_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_NFE_nhomalt = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_POPMAX_AC = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_POPMAX_AF = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_POPMAX_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_POPMAX_nhomalt = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_SAS_AC = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_SAS_AF = functionalAnnotationArray[pos++].trim();
        // 309
        String gnomAD_exomes_controls_SAS_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_SAS_nhomalt = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_controls_nhomalt = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_flag = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_nhomalt = functionalAnnotationArray[pos++].trim();
        addNumberToJsonObject("AC", functionalAnnotationArray[pos++] , frequencyGnomadGen, true);
        addNumberToJsonObject("AF", functionalAnnotationArray[pos++] , frequencyGnomadGen, true);
        addNumberToJsonObject("AFR_AC", functionalAnnotationArray[pos++] , frequencyGnomadGen, true);
        addNumberToJsonObject("AFR_AF", functionalAnnotationArray[pos++] , frequencyGnomadGen, true);
        String gnomAD_genomes_AFR_AN = functionalAnnotationArray[pos++].trim();
        // 319
        String gnomAD_genomes_AFR_nhomalt = functionalAnnotationArray[pos++].trim();
        addNumberToJsonObject("AMR_AC", functionalAnnotationArray[pos++] , frequencyGnomadGen, true);
        addNumberToJsonObject("AMR_AF", functionalAnnotationArray[pos++] , frequencyGnomadGen, true);
        String gnomAD_genomes_AMR_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_AMR_nhomalt = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_ASJ_AC = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_ASJ_AF = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_ASJ_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_ASJ_nhomalt = functionalAnnotationArray[pos++].trim();
        // 329
        addNumberToJsonObject("EAS_AC", functionalAnnotationArray[pos++] , frequencyGnomadGen, true);
        addNumberToJsonObject("EAS_AF", functionalAnnotationArray[pos++] , frequencyGnomadGen, true);
        String gnomAD_genomes_EAS_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_EAS_nhomalt = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_FIN_AC = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_FIN_AF = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_FIN_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_FIN_nhomalt = functionalAnnotationArray[pos++].trim();
        addNumberToJsonObject("NFE_AC", functionalAnnotationArray[pos++] , frequencyGnomadGen, true);
        addNumberToJsonObject("NFE_AF", functionalAnnotationArray[pos++] , frequencyGnomadGen, true);
        // 339
        String gnomAD_genomes_NFE_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_NFE_nhomalt = functionalAnnotationArray[pos++].trim();
        addNumberToJsonObject("POPMAX_AC", functionalAnnotationArray[pos++] , frequencyGnomadGen, true);
        addNumberToJsonObject("POPMAX_AF", functionalAnnotationArray[pos++] , frequencyGnomadGen, true);
        String gnomAD_genomes_POPMAX_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_POPMAX_nhomalt = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_controls_AC = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_controls_AF = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_controls_AFR_AC = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_controls_AFR_AF = functionalAnnotationArray[pos++].trim();
        // 349
        String gnomAD_genomes_controls_AFR_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_controls_AFR_nhomalt = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_controls_AMR_AC = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_controls_AMR_AF = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_controls_AMR_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_controls_AMR_nhomalt = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_controls_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_controls_ASJ_AC = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_controls_ASJ_AF = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_controls_ASJ_AN = functionalAnnotationArray[pos++].trim();
        // 359
        String gnomAD_genomes_controls_ASJ_nhomalt = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_controls_EAS_AC = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_controls_EAS_AF = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_controls_EAS_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_controls_EAS_nhomalt = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_controls_FIN_AC = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_controls_FIN_AF = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_controls_FIN_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_controls_FIN_nhomalt = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_controls_NFE_AC = functionalAnnotationArray[pos++].trim();
        // 369
        String gnomAD_genomes_controls_NFE_AF = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_controls_NFE_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_controls_NFE_nhomalt = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_controls_POPMAX_AC = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_controls_POPMAX_AF = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_controls_POPMAX_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_controls_POPMAX_nhomalt = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_controls_nhomalt = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_flag = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_nhomalt = functionalAnnotationArray[pos++].trim();
        // 379
        String hg18_chr = functionalAnnotationArray[pos++].trim();
        String hg18_pos_1based = functionalAnnotationArray[pos++].trim();
        String hg19_chr = functionalAnnotationArray[pos++].trim();
        String hg19_pos_1based = functionalAnnotationArray[pos++].trim();
        String integrated_confidence_value = functionalAnnotationArray[pos++].trim();
        String integrated_fitCons_rankscore = functionalAnnotationArray[pos++].trim();
        String integrated_fitCons_score = functionalAnnotationArray[pos++].trim();
        String phastCons100way_vertebrate = functionalAnnotationArray[pos++].trim();
        String phastCons100way_vertebrate_rankscore = functionalAnnotationArray[pos++].trim();
        String phastCons17way_primate = functionalAnnotationArray[pos++].trim();
        // 389
        String phastCons17way_primate_rankscore = functionalAnnotationArray[pos++].trim();
        String phastCons30way_mammalian = functionalAnnotationArray[pos++].trim();
        String phastCons30way_mammalian_rankscore = functionalAnnotationArray[pos++].trim();
        String phyloP100way_vertebrate = functionalAnnotationArray[pos++].trim();
        String phyloP100way_vertebrate_rankscore = functionalAnnotationArray[pos++].trim();
        String phyloP17way_primate = functionalAnnotationArray[pos++].trim();
        String phyloP17way_primate_rankscore = functionalAnnotationArray[pos++].trim();
        String phyloP30way_mammalian = functionalAnnotationArray[pos++].trim();
        String phyloP30way_mammalian_rankscore = functionalAnnotationArray[pos++].trim();
        String pos_1_based = functionalAnnotationArray[pos++].trim();
        // 399
//        String ref = functionalAnnotationArray[pos++].trim();
        addStrToJsonObject("ref", functionalAnnotationArray[pos++], funcAnnotation, false);
//        String refcodon = functionalAnnotationArray[pos++].trim();
        addStrToJsonObject("ref_codon", functionalAnnotationArray[pos++], funcAnnotation, false);
        String rs_dbSNP151 = functionalAnnotationArray[pos++].trim();
        // 402

        //System.out.println("csqArray.length="+csqArray.length);
        //String[] functionalAnnotationArray = csqArray[i].split("|");


        /* population frequencies Json Obj */



        frequencies.put("1000Gp3", frequency1000Gp3);
        frequencies.put("ExAc", frequencyExAc);
        frequencies.put("gnomAD_exomes", frequencyGnomadEx);
        frequencies.put("gnomAD_genomes", frequencyGnomadGen);
        frequencies.put("ESP6500", frequencyEsp6500);
        funcAnnotation.put("frequencies", frequencies);
        funcAnnotation.put("predictions", prediction);
        //funcAnnotation.put("properties", funcAnnoProperties);
        //funcAnnotation.(funcAnnoProperties);


        if ((boolean) frequencyGnomadEx.get("available")) System.out.println(funcAnnotation.toString(2));
        //if ((boolean) frequencyGnomadGen.get("available")) System.out.println(frequencyGnomadGen.toString(2));

//        System.out.print("\n");

        return funcAnnotation;

    }

    private static void addNumberToJsonObject(String popName, String var, JSONObject freqObject, boolean withAvailability) {
        if (var.length()>0) {
            try {
                freqObject.put(popName, NF.parse(var));
                if (withAvailability) freqObject.put("available", true);
            } catch (ParseException parseException) {
                parseException.printStackTrace();
                if (withAvailability) freqObject.put("available", false);
            }
        } else if (withAvailability) {
            freqObject.put("available", false);
        }
    }

    private static void addStrToJsonObject(String popName, String var, JSONObject freqObject, boolean withAvailability) {
        if (var.length()>0) {
                freqObject.put(popName, var);
                if (withAvailability) freqObject.put("available", true);
        } else if (withAvailability) {

            freqObject.put("available", false);
        }
    }



}
