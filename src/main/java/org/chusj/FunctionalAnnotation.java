package org.chusj;

import org.json.JSONArray;
import org.json.JSONObject;

import java.util.Objects;

public class FunctionalAnnotation {

    private String gene;
    private String aaChange;
    private String consequence;
    private String codingDNAChange;
    private Long strand;
    private JSONArray theRest;
    private String geneId;
    private String impact;
    private JSONObject scores;
    private JSONObject predictions;
    private String biotype;

    public String getBiotype() {
        return biotype;
    }

    public void setBiotype(String biotype) {
        this.biotype = biotype;
    }

    public JSONObject getScores() {
        return scores;
    }

    void setScores(JSONObject scores) {
        this.scores = scores;
    }

    JSONObject getPredictions() {
        return predictions;
    }

    void setPredictions(JSONObject predictions) {
        this.predictions = predictions;
    }

    String getImpact() {
        return impact;
    }

    void setImpact(String impact) {
        this.impact = impact;
    }

    public FunctionalAnnotation() {
    }

    FunctionalAnnotation(String gene, String aaChange, String consequence, String codingDNAChange, Long strand) {
        this.gene = gene;
        this.aaChange = aaChange;
        this.consequence = consequence;
        this.codingDNAChange = codingDNAChange;
        this.strand = strand;
    }

    public String getGene() {
        return gene;
    }

    public void setGene(String gene) {
        this.gene = gene;
    }

    private String getAaChange() {
        return aaChange;
    }

    public void setAaChange(String aaChange) {
        this.aaChange = aaChange;
    }

    private String getConsequence() {
        return consequence;
    }

    public void setConsequence(String consequence) {
        this.consequence = consequence;
    }

    private String getCodingDNAChange() {
        return codingDNAChange;
    }

    public void setCodingDNAChange(String codingDNAChange) {
        this.codingDNAChange = codingDNAChange;
    }

    public String getGeneId() {
        return geneId;
    }

    void setGeneId(String geneId) {
        this.geneId = geneId;
    }

    private Long getStrand() {
        return strand;
    }

    public void setStrand(Long strand) {
        this.strand = strand;
    }

    JSONArray getTheRest() {
        return theRest;
    }

    void setTheRest(JSONArray theRest) {
        this.theRest = theRest;
    }



    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof FunctionalAnnotation)) return false;
        FunctionalAnnotation that = (FunctionalAnnotation) o;
        return Objects.equals(getStrand(), that.getStrand()) &&
                Objects.equals(getGene(), that.getGene()) &&
                Objects.equals(getAaChange(), that.getAaChange()) &&
                Objects.equals(getConsequence(), that.getConsequence()) &&
                Objects.equals(getCodingDNAChange(), that.getCodingDNAChange());
    }

    @Override
    public int hashCode() {
        return Objects.hash(getGene(), getAaChange(), getConsequence(), getCodingDNAChange(), getStrand());
    }

    JSONObject getJsonObj() {

        JSONObject fa = new JSONObject();
        if (!gene.isEmpty()) fa.put("geneAffectedSymbol", gene);
        if (!geneId.isEmpty()) fa.put("geneAffectedId", geneId);
        if (!aaChange.isEmpty()) fa.put("aaChange", aaChange);
        fa.put("consequence", consequence);
        if (!codingDNAChange.isEmpty()) fa.put("codingDNAChange", codingDNAChange);
        fa.put("impact", impact);
        if (strand != null ) {
            fa.put("strand", strand);
        }
        fa.put("transcripts", theRest);
        if (scores != null && (boolean) scores.get("available")) fa.put("conservationsScores", scores);
        if (predictions != null && (boolean) predictions.get("available")) fa.put("predictions", predictions);
        fa.put("biotype", biotype);


        return fa;

    }



}
