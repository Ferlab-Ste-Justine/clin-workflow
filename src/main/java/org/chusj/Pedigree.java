package org.chusj;

public class Pedigree {

    String familyId;
    String id;
    String paternalId;
    String maternalId;
    String sex;
    String phenotype;


    public String getFamilyId() {
        return familyId;
    }

    public void setFamilyId(String familyId) {
        this.familyId = familyId;
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getPaternalId() {
        return paternalId;
    }

    public void setPaternalId(String paternalId) {
        this.paternalId = paternalId;
    }

    public String getMaternalId() {
        return maternalId;
    }

    public void setMaternalId(String maternalId) {
        this.maternalId = maternalId;
    }

    public String getSex() {
        return sex;
    }

    public void setSex(String sex) {
        this.sex = sex;
    }

    public String getPhenotype() {
        return phenotype;
    }

    public void setPhenotype(String phenotype) {
        this.phenotype = phenotype;
    }

    @Override
    public String toString() {
        return "Pedigree{" +
                "familyId='" + familyId + '\'' +
                ", id='" + id + '\'' +
                ", paternalId='" + paternalId + '\'' +
                ", maternalId='" + maternalId + '\'' +
                ", sex='" + sex + '\'' +
                ", phenotype='" + phenotype + '\'' +
                '}';
    }
}
