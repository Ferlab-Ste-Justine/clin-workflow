package org.chusj;

import java.io.Serializable;

public class GeneVariants implements Serializable {

    private Gene gene;
    private Variant variant;


    public Gene getGene() {
        return gene;
    }

    public void setGene(Gene gene) {
        this.gene = gene;
    }

    public Variant getVariant() {
        return variant;
    }

    public void setVariant(Variant variant) {
        this.variant = variant;
    }


    @Override
    public String toString() {
        return "GeneVariants{" +
                "gene=" + gene +
                ", variant=" + variant +
                '}';
    }
}
