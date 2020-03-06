package org.chusj;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Objects;

public class Gene implements Serializable {

    private String ensemblId;
    private String geneSymbol;
    private String biotype;
    private String[] HPO;
    private String[] radboudumc;
    private String[] omim;
    private String[] orphanet;


    public Gene(String ensemblId, String geneSymbol, String biotype, String[] HPO, String[] radboudumc, String[] omim, String[] orphanet) {
        this.ensemblId = ensemblId;
        this.geneSymbol = geneSymbol;
        this.biotype = biotype;
        this.HPO = HPO;
        this.radboudumc = radboudumc;
        this.omim = omim;
        this.orphanet = orphanet;
    }

    public Gene(String ensemblId, String geneSymbol, String biotype) {
        this.ensemblId = ensemblId;
        this.geneSymbol = geneSymbol;
        this.biotype = biotype;
    }


    public String getEnsemblId() {
        return ensemblId;
    }

    public void setEnsemblId(String ensemblId) {
        this.ensemblId = ensemblId;
    }

    public String getGeneSymbol() {
        return geneSymbol;
    }

    public void setGeneSymbol(String geneSymbol) {
        this.geneSymbol = geneSymbol;
    }

    public String getBiotype() {
        return biotype;
    }

    public void setBiotype(String biotype) {
        this.biotype = biotype;
    }

    public String[] getHPO() {
        return HPO;
    }

    public void setHPO(String[] HPO) {
        this.HPO = HPO;
    }

    public String[] getRadboudumc() {
        return radboudumc;
    }

    public void setRadboudumc(String[] radboudumc) {
        this.radboudumc = radboudumc;
    }

    public String[] getOmim() {
        return omim;
    }

    public void setOmim(String[] omim) {
        this.omim = omim;
    }

    public String[] getOrphanet() {
        return orphanet;
    }

    public void setOrphanet(String[] orphanet) {
        this.orphanet = orphanet;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Gene gene = (Gene) o;
        return ensemblId.equals(gene.ensemblId) &&
                geneSymbol.equals(gene.geneSymbol);
    }

    @Override
    public int hashCode() {
        int result = Objects.hash(ensemblId, geneSymbol);
        result = 31 * result + Arrays.hashCode(HPO);
        result = 31 * result + Arrays.hashCode(radboudumc);
        result = 31 * result + Arrays.hashCode(omim);
        result = 31 * result + Arrays.hashCode(orphanet);
        return result;
    }

    @Override
    public String toString() {
        return "geneSymbol=" + geneSymbol ;
    }
}
