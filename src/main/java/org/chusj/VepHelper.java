package org.chusj;


import org.json.JSONObject;

import java.io.BufferedReader;
import java.io.FileReader;
import java.text.NumberFormat;
import java.text.ParseException;

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


    private static void processVcfDataLine(String extractedLine) {

        //CHROM	POS	ID	REF	ALT	QUAL	FILTER	DP	MQ	MQRankSum
        // ReadPosRankSum	LOD	FractionInformativeReads	SNP	MNP	INS	DEL	MIXED	HOM	HET
        // GEN[*].GT	GEN[*].AD	GEN[*].AF	GEN[*].F1R2	GEN[*].F2R1	GEN[*].DP	GEN[*].SB	GEN[*].MB	CSQ
        //System.out.print(".");
        String[] lineValueArray = extractedLine.split("\t");
        // System.out.println("lineValueArray.length="+lineValueArray.length);
        // dynamic positioning system -- pos counter++
        int pos = 0;

        String chrom = lineValueArray[pos++];
        if ("CHROM".equalsIgnoreCase(chrom)) {
            return; // Meta data line
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

        JSONObject funcAnnotation = null;

        for (String s : csqArray) {
            funcAnnotation = processVepAnnotations(s, dnaChanges);

        }






//        String ac = lineValueArray[pos++].trim();
//        String an = lineValueArray[pos++].trim();
//        String db = lineValueArray[pos++].trim();
//        String gq = lineValueArray[pos++].trim();
//        String ac2 = lineValueArray[pos++].trim();

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
        JSONObject frequencies = new JSONObject();
        JSONObject frequencyExAc = new JSONObject();
        JSONObject frequency1000Gp3 = new JSONObject();
        JSONObject frequencyEsp6500 = new JSONObject();
        JSONObject frequencyGnomadEx = new JSONObject();
        JSONObject frequencyGnomadGen = new JSONObject();

        JSONObject predictions = new JSONObject();
        JSONObject prediction = new JSONObject();



        //0 -
        String Allele = functionalAnnotationArray[pos++].trim();
        String Consequence = functionalAnnotationArray[pos++].trim();
        String IMPACT = functionalAnnotationArray[pos++].trim();
        String SYMBOL = functionalAnnotationArray[pos++].trim();
        String Gene = functionalAnnotationArray[pos++].trim();
        String Feature_type = functionalAnnotationArray[pos++].trim();
        String Feature = functionalAnnotationArray[pos++].trim();
        String BIOTYPE = functionalAnnotationArray[pos++].trim();
        String EXON = functionalAnnotationArray[pos++].trim();
        String INTRON = functionalAnnotationArray[pos++].trim();
        // 9 -
        String HGVSc = functionalAnnotationArray[pos++].trim();
        String HGVSp = functionalAnnotationArray[pos++].trim();
        String cDNA_position = functionalAnnotationArray[pos++].trim();
        String CDS_position = functionalAnnotationArray[pos++].trim();
        String Protein_position = functionalAnnotationArray[pos++].trim();
        String Amino_acids = functionalAnnotationArray[pos++].trim();
        String Codons = functionalAnnotationArray[pos++].trim();
        String Existing_variation = functionalAnnotationArray[pos++].trim();
        String DISTANCE = functionalAnnotationArray[pos++].trim();
        String STRAND = functionalAnnotationArray[pos++].trim();
        //19 -
        String FLAGS = functionalAnnotationArray[pos++].trim();
        String SYMBOL_SOURCE = functionalAnnotationArray[pos++].trim();
        String HGNC_ID = functionalAnnotationArray[pos++].trim();
        String CANONICAL = functionalAnnotationArray[pos++].trim();
        String GIVEN_REF = functionalAnnotationArray[pos++].trim();
        String USED_REF = functionalAnnotationArray[pos++].trim();
        String BAM_EDIT = functionalAnnotationArray[pos++].trim();

        addNumberToJsonObject("AC", functionalAnnotationArray[pos++] , frequency1000Gp3);
        addNumberToJsonObject("AF", functionalAnnotationArray[pos++] , frequency1000Gp3);
        addNumberToJsonObject("AFR_AC", functionalAnnotationArray[pos++] , frequency1000Gp3);
        // 29
        addNumberToJsonObject("AFR_AF", functionalAnnotationArray[pos++] , frequency1000Gp3);
        addNumberToJsonObject("AMR_AC", functionalAnnotationArray[pos++] , frequency1000Gp3);
        addNumberToJsonObject("AMR_AF", functionalAnnotationArray[pos++] , frequency1000Gp3);
        addNumberToJsonObject("EAS_AC", functionalAnnotationArray[pos++] , frequency1000Gp3);
        addNumberToJsonObject("EAS_AF", functionalAnnotationArray[pos++] , frequency1000Gp3);
        addNumberToJsonObject("EUR_AC", functionalAnnotationArray[pos++] , frequency1000Gp3);
        addNumberToJsonObject("EUR_AF", functionalAnnotationArray[pos++] , frequency1000Gp3);
        addNumberToJsonObject("SAS_AC", functionalAnnotationArray[pos++] , frequency1000Gp3);
        addNumberToJsonObject("SAS_AF", functionalAnnotationArray[pos++] , frequency1000Gp3);
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
        addNumberToJsonObject("AA_AC", functionalAnnotationArray[pos++] , frequencyEsp6500);
        //59 -
        addNumberToJsonObject("AA_AF", functionalAnnotationArray[pos++] , frequencyEsp6500);
        addNumberToJsonObject("EA_AC", functionalAnnotationArray[pos++] , frequencyEsp6500);
        addNumberToJsonObject("EA_AF", functionalAnnotationArray[pos++] , frequencyEsp6500);
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
        addNumberToJsonObject("AC", functionalAnnotationArray[pos++] , frequencyExAc);
        addNumberToJsonObject("AF", functionalAnnotationArray[pos++] , frequencyExAc);
        addNumberToJsonObject("AFR_AC", functionalAnnotationArray[pos++] , frequencyExAc);
        addNumberToJsonObject("AFR_AF", functionalAnnotationArray[pos++] , frequencyExAc);
        addNumberToJsonObject("AMR_AC", functionalAnnotationArray[pos++] , frequencyExAc);
        addNumberToJsonObject("AMR_AF", functionalAnnotationArray[pos++] , frequencyExAc);
        addNumberToJsonObject("Adj_AC", functionalAnnotationArray[pos++] , frequencyExAc);
        addNumberToJsonObject("Adj_AF", functionalAnnotationArray[pos++] , frequencyExAc);
        // 79
        addNumberToJsonObject("EAS_AC", functionalAnnotationArray[pos++] , frequencyExAc);
        addNumberToJsonObject("EAS_AF", functionalAnnotationArray[pos++] , frequencyExAc);
        addNumberToJsonObject("FIN_AC", functionalAnnotationArray[pos++] , frequencyExAc);
        addNumberToJsonObject("FIN_AF", functionalAnnotationArray[pos++] , frequencyExAc);
        addNumberToJsonObject("NFE_AC", functionalAnnotationArray[pos++] , frequencyExAc);
        addNumberToJsonObject("NFE_AF", functionalAnnotationArray[pos++] , frequencyExAc);
        addNumberToJsonObject("SAS_AC", functionalAnnotationArray[pos++] , frequencyExAc);
        addNumberToJsonObject("SAS_AF", functionalAnnotationArray[pos++] , frequencyExAc);
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
        String FATHMM_pred = functionalAnnotationArray[pos++].trim();
        String FATHMM_score = functionalAnnotationArray[pos++].trim();
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
        String Polyphen2_HDIV_pred = functionalAnnotationArray[pos++].trim();
        String Polyphen2_HDIV_rankscore = functionalAnnotationArray[pos++].trim();
        String Polyphen2_HDIV_score = functionalAnnotationArray[pos++].trim();
        String Polyphen2_HVAR_pred = functionalAnnotationArray[pos++].trim();
        String Polyphen2_HVAR_rankscore = functionalAnnotationArray[pos++].trim();
        String Polyphen2_HVAR_score = functionalAnnotationArray[pos++].trim();
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
        String SIFT_pred = functionalAnnotationArray[pos++].trim();
        // 199
        String SIFT_score = functionalAnnotationArray[pos++].trim();
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
        String aaalt = functionalAnnotationArray[pos++].trim();
        String aapos = functionalAnnotationArray[pos++].trim();
        String aaref = functionalAnnotationArray[pos++].trim();
        String alt = functionalAnnotationArray[pos++].trim();
        String bStatistic = functionalAnnotationArray[pos++].trim();
        // 219
        String bStatistic_rankscore = functionalAnnotationArray[pos++].trim();
        String cds_strand = functionalAnnotationArray[pos++].trim();
        String chr = functionalAnnotationArray[pos++].trim();
        String clinvar_MedGen_id = functionalAnnotationArray[pos++].trim();
        String clinvar_OMIM_id = functionalAnnotationArray[pos++].trim();
        String clinvar_Orphanet_id = functionalAnnotationArray[pos++].trim();
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
        String genename = functionalAnnotationArray[pos++].trim();

        addNumberToJsonObject("AC", functionalAnnotationArray[pos++] , frequencyGnomadEx);
        addNumberToJsonObject("AF", functionalAnnotationArray[pos++] , frequencyGnomadEx);
        addNumberToJsonObject("AFR_AC", functionalAnnotationArray[pos++] , frequencyGnomadEx);
        addNumberToJsonObject("AFR_AF", functionalAnnotationArray[pos++] , frequencyGnomadEx);

        String gnomAD_exomes_AFR_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_AFR_nhomalt = functionalAnnotationArray[pos++].trim();

        addNumberToJsonObject("AMR_AC", functionalAnnotationArray[pos++] , frequencyGnomadEx);
        addNumberToJsonObject("AMR_AF", functionalAnnotationArray[pos++] , frequencyGnomadEx);
        // 249
        String gnomAD_exomes_AMR_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_AMR_nhomalt = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_AN = functionalAnnotationArray[pos++].trim();

        addNumberToJsonObject("ASJ_AC", functionalAnnotationArray[pos++] , frequencyGnomadEx);
        addNumberToJsonObject("ASJ_AF", functionalAnnotationArray[pos++] , frequencyGnomadEx);
        String gnomAD_exomes_ASJ_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_ASJ_nhomalt = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_EAS_AC = functionalAnnotationArray[pos++].trim();

        addNumberToJsonObject("EAS_AF", functionalAnnotationArray[pos++] , frequencyGnomadEx);
        addNumberToJsonObject("EAS_AN", functionalAnnotationArray[pos++] , frequencyGnomadEx);
        // 259
        String gnomAD_exomes_EAS_nhomalt = functionalAnnotationArray[pos++].trim();

        addNumberToJsonObject("FIN_AC", functionalAnnotationArray[pos++] , frequencyGnomadEx);
        addNumberToJsonObject("FIN_AF", functionalAnnotationArray[pos++] , frequencyGnomadEx);
        String gnomAD_exomes_FIN_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_FIN_nhomalt = functionalAnnotationArray[pos++].trim();
        addNumberToJsonObject("NFE_AC", functionalAnnotationArray[pos++] , frequencyGnomadEx);
        addNumberToJsonObject("NFE_AF", functionalAnnotationArray[pos++] , frequencyGnomadEx);
        String gnomAD_exomes_NFE_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_NFE_nhomalt = functionalAnnotationArray[pos++].trim();
        // 269
        addNumberToJsonObject("POPMAX_AC", functionalAnnotationArray[pos++] , frequencyGnomadEx);
        // 269
        addNumberToJsonObject("POPMAX_AF", functionalAnnotationArray[pos++] , frequencyGnomadEx);
        String gnomAD_exomes_POPMAX_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_exomes_POPMAX_nhomalt = functionalAnnotationArray[pos++].trim();
        addNumberToJsonObject("SAS_AC", functionalAnnotationArray[pos++] , frequencyGnomadEx);
        addNumberToJsonObject("SAS_AF", functionalAnnotationArray[pos++] , frequencyGnomadEx);
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
        addNumberToJsonObject("AC", functionalAnnotationArray[pos++] , frequencyGnomadGen);
        addNumberToJsonObject("AF", functionalAnnotationArray[pos++] , frequencyGnomadGen);
        addNumberToJsonObject("AFR_AC", functionalAnnotationArray[pos++] , frequencyGnomadGen);
        addNumberToJsonObject("AFR_AF", functionalAnnotationArray[pos++] , frequencyGnomadGen);
        String gnomAD_genomes_AFR_AN = functionalAnnotationArray[pos++].trim();
        // 319
        String gnomAD_genomes_AFR_nhomalt = functionalAnnotationArray[pos++].trim();
        addNumberToJsonObject("AMR_AC", functionalAnnotationArray[pos++] , frequencyGnomadGen);
        addNumberToJsonObject("AMR_AF", functionalAnnotationArray[pos++] , frequencyGnomadGen);
        String gnomAD_genomes_AMR_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_AMR_nhomalt = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_ASJ_AC = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_ASJ_AF = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_ASJ_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_ASJ_nhomalt = functionalAnnotationArray[pos++].trim();
        // 329
        addNumberToJsonObject("EAS_AC", functionalAnnotationArray[pos++] , frequencyGnomadGen);
        addNumberToJsonObject("EAS_AF", functionalAnnotationArray[pos++] , frequencyGnomadGen);
        String gnomAD_genomes_EAS_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_EAS_nhomalt = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_FIN_AC = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_FIN_AF = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_FIN_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_FIN_nhomalt = functionalAnnotationArray[pos++].trim();
        addNumberToJsonObject("NFE_AC", functionalAnnotationArray[pos++] , frequencyGnomadGen);
        addNumberToJsonObject("NFE_AF", functionalAnnotationArray[pos++] , frequencyGnomadGen);
        // 339
        String gnomAD_genomes_NFE_AN = functionalAnnotationArray[pos++].trim();
        String gnomAD_genomes_NFE_nhomalt = functionalAnnotationArray[pos++].trim();
        addNumberToJsonObject("POPMAX_AC", functionalAnnotationArray[pos++] , frequencyGnomadGen);
        addNumberToJsonObject("POPMAX_AF", functionalAnnotationArray[pos++] , frequencyGnomadGen);
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
        String ref = functionalAnnotationArray[pos++].trim();
        String refcodon = functionalAnnotationArray[pos++].trim();
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



        if ((boolean) frequencyGnomadEx.get("available")) System.out.println(funcAnnotation.toString(2));
        //if ((boolean) frequencyGnomadGen.get("available")) System.out.println(frequencyGnomadGen.toString(2));

//        System.out.print("\n");

        return; funcAnnotation;

    }

    private static void addNumberToJsonObject(String popName, String var, JSONObject freqObject) {
        if (var.length()>0) {
            try {
                freqObject.put(popName, NF.parse(var));
                freqObject.put("available", true);
            } catch (ParseException parseException) {
                parseException.printStackTrace();
                freqObject.put("available", false);
            }
        } else {
            freqObject.put("available", false);
        }
    }


}
