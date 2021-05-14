package org.chusj;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

public class Variant implements Serializable {

    private List<Gene> genes = new ArrayList<>();
    private String mutationId;
    private String jsonObjInString;
    private String id;

    public List<Gene> getGenes() {
        return genes;
    }

    public void setGenes(List<Gene> genes) {
        this.genes = genes;
    }

    public String getMutationId() {
        return mutationId;
    }

    public void setMutationId(String mutationId) {
        this.mutationId = mutationId;
    }

    public String getJsonObjInString() {
        return jsonObjInString;
    }

    public void setJsonObjInString(String jsonObjInString) {
        this.jsonObjInString = jsonObjInString;
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public void addGene(Gene gene) {
        this.genes.add(gene);
    }

    @Override
    public String toString() {
        return "mutationId=" + mutationId;
    }
}
