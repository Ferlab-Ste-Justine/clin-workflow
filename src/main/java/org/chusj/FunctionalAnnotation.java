package org.chusj;

import org.json.JSONArray;
import org.json.JSONObject;

import java.util.*;

public class FunctionalAnnotation {

    private String gene;
    private String aaChange;
    private Set<String> consequence = new HashSet<>();
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

    void setBiotype(String biotype) {
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
        this.consequence = new HashSet<>(Arrays.asList(consequence.split("&")));
        this.codingDNAChange = codingDNAChange;
        this.strand = strand;
    }

    private String getGene() {
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

    private Set<String> getConsequence() {
        return consequence;
    }

    public void setConsequence(Set<String> consequence) {
        this.consequence = consequence;
    }

    public void addConsequence(String consequence) { this.consequence.add(consequence); }

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
        if (geneId != null && !geneId.isEmpty()) fa.put("geneAffectedId", geneId);
        if (!aaChange.isEmpty()) fa.put("aaChange", aaChange);
        fa.put("consequence", consequence.toArray());
        if (!codingDNAChange.isEmpty()) fa.put("cdnaChange", codingDNAChange);
        fa.put("impact", impact);
        if (strand != null ) {
            fa.put("strand", strand);
        }
        fa.put("transcripts", theRest);
        if (scores != null && (boolean) scores.get("available")) {
            scores.remove("available");
            fa.put("conservationsScores", scores);
        }
        if (predictions != null ) { //&& !predictions.isNull("available")
            fa.remove("available");
            fa.put("predictions", predictions);
        }
        fa.put("biotype", biotype);


        return fa;

    }



}
