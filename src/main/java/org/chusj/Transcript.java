package org.chusj;

import org.json.JSONObject;

import static org.chusj.VepHelper.addStrToJsonObject;

public class Transcript {

    private String ensemblTranscriptId;
    private String aaPosition;
    private String FATHMM;
    private String SIFT;
    private String SIFT_score;
    private String ensemblProteinId;
    private String Polyphen2_HVAR_score;
    private String MutationAssessor_score;
    private String Polyphen2_HDIV;
    private String Polyphen2_HDIV_score;
    private String LRT_pred;
    private String FATHMM_score;
    private String MutationAssessor_pred;
    private String Polyphen2_HVAR_pred;
    private String LRT_score;

    public Transcript(String ensemblTranscriptId, String aaPosition, String FATHMM, String SIFT, String ensemblProteinId,
                      String polyphen2_HVAR_score, String mutationAssessor_score, String polyphen2_HDIV,
                      String polyphen2_HDIV_score, String LRT_pred, String FATHMM_score, String mutationAssessor_pred,
                      String polyphen2_HVAR_pred, String LRT_score, String SIFT_score) {
        this.ensemblTranscriptId = ensemblTranscriptId;
        this.aaPosition = aaPosition;
        this.FATHMM = FATHMM;
        this.SIFT = SIFT;
        this.ensemblProteinId = ensemblProteinId;
        Polyphen2_HVAR_score = polyphen2_HVAR_score;
        MutationAssessor_score = mutationAssessor_score;
        Polyphen2_HDIV = polyphen2_HDIV;
        Polyphen2_HDIV_score = polyphen2_HDIV_score;
        this.LRT_pred = LRT_pred;
        this.FATHMM_score = FATHMM_score;
        MutationAssessor_pred = mutationAssessor_pred;
        Polyphen2_HVAR_pred = polyphen2_HVAR_pred;
        this.LRT_score = LRT_score;
        this.SIFT_score = SIFT_score;
    }

    public Transcript() {
    }

    public Transcript(String ensemblTranscriptId) {
        this.ensemblTranscriptId = ensemblTranscriptId;
    }

    public String getEnsemblTranscriptId() {
        return ensemblTranscriptId;
    }

    public void setEnsemblTranscriptId(String ensemblTranscriptId) {
        this.ensemblTranscriptId = ensemblTranscriptId;
    }

    public String getAaPosition() {
        return aaPosition;
    }

    public void setAaPosition(String aaPosition) {
        this.aaPosition = aaPosition;
    }

    public String getFATHMM() {
        return FATHMM;
    }

    public void setFATHMM(String FATHMM) {
        this.FATHMM = FATHMM;
    }

    public String getSIFT() {
        return SIFT;
    }

    public void setSIFT(String SIFT) {
        this.SIFT = SIFT;
    }

    public String getEnsemblProteinId() {
        return ensemblProteinId;
    }

    public void setEnsemblProteinId(String ensemblProteinId) {
        this.ensemblProteinId = ensemblProteinId;
    }

    public String getPolyphen2_HVAR_score() {
        return Polyphen2_HVAR_score;
    }

    public void setPolyphen2_HVAR_score(String polyphen2_HVAR_score) {
        Polyphen2_HVAR_score = polyphen2_HVAR_score;
    }

    public String getMutationAssessor_score() {
        return MutationAssessor_score;
    }

    public void setMutationAssessor_score(String mutationAssessor_score) {
        MutationAssessor_score = mutationAssessor_score;
    }

    public String getPolyphen2_HDIV() {
        return Polyphen2_HDIV;
    }

    public void setPolyphen2_HDIV(String polyphen2_HDIV) {
        Polyphen2_HDIV = polyphen2_HDIV;
    }

    public String getPolyphen2_HDIV_score() {
        return Polyphen2_HDIV_score;
    }

    public void setPolyphen2_HDIV_score(String polyphen2_HDIV_score) {
        Polyphen2_HDIV_score = polyphen2_HDIV_score;
    }

    public String getLRT_pred() {
        return LRT_pred;
    }

    public void setLRT_pred(String LRT_pred) {
        this.LRT_pred = LRT_pred;
    }

    public String getFATHMM_score() {
        return FATHMM_score;
    }

    public void setFATHMM_score(String FATHMM_score) {
        this.FATHMM_score = FATHMM_score;
    }

    public String getMutationAssessor_pred() {
        return MutationAssessor_pred;
    }

    public void setMutationAssessor_pred(String mutationAssessor_pred) {
        MutationAssessor_pred = mutationAssessor_pred;
    }

    public String getPolyphen2_HVAR_pred() {
        return Polyphen2_HVAR_pred;
    }

    public void setPolyphen2_HVAR_pred(String polyphen2_HVAR_pred) {
        Polyphen2_HVAR_pred = polyphen2_HVAR_pred;
    }

    public String getLRT_score() {
        return LRT_score;
    }

    public void setLRT_score(String LRT_score) {
        this.LRT_score = LRT_score;
    }



    public JSONObject getPrediction() {
        JSONObject prediction = new JSONObject();

        addStrToJsonObject("FATHMM_score", FATHMM_score, prediction, true);
        addStrToJsonObject("FATHMM", FATHMM, prediction, true);
        addStrToJsonObject("SIFT", SIFT, prediction, true);
        addStrToJsonObject("SIFT_score", SIFT_score, prediction, true);
        addStrToJsonObject("MutationAssessor_pred", MutationAssessor_pred, prediction, true);
        addStrToJsonObject("MutationAssessor_score", MutationAssessor_score, prediction, true);
        addStrToJsonObject("Polyphen2_HDIV", Polyphen2_HDIV, prediction, true);
        addStrToJsonObject("Polyphen2_HDIV_score", Polyphen2_HDIV_score, prediction, true);
        addStrToJsonObject("Polyphen2_HVAR_pred", Polyphen2_HVAR_pred, prediction, true);
        addStrToJsonObject("Polyphen2_HVAR_score", Polyphen2_HVAR_score, prediction, true);
        addStrToJsonObject("LRT_pred", LRT_pred, prediction, true);
        addStrToJsonObject("LRT_score", LRT_score, prediction, true);


        return prediction;
    }


}
