package org.chusj;

import org.json.JSONObject;

import java.io.Serializable;
import java.util.Arrays;
import java.util.List;

public class Patient implements Serializable {

    private String specimenId;
    private String relation;
    private String patientId;
    private String familyId;
    private String studyId;
    private String orgId;
    private String labName;
    private String practitionerId;
    private boolean isProban;
    private boolean isAffected;
    private JSONObject patient;
    private List<String> hposTerms;
    private int qtyHposTermsFound=0;
    private String sequencingStrategy;
    private String[] linkToOthersPatientId;
    private String[] linkToOthersSpecimenId;
    private String gender;


    Patient(String specimenId) {
        this.specimenId = specimenId;
    }


    String getSpecimenId() {
        return specimenId;
    }

    void setSpecimenId(String specimenId) {
        this.specimenId = specimenId;
    }

    String getRelation() {
        return relation;
    }

    void setRelation(String relation) {
        this.relation = relation;
    }

    String getGender() {
        return gender;
    }

    void setGender(String gender) {
        this.gender = gender;
    }

    String[] getLinkToOthersSpecimenId() {
        return linkToOthersSpecimenId;
    }

    void setLinkToOthersSpecimenId(String[] linkToOthersSpecimenId) {
        this.linkToOthersSpecimenId = linkToOthersSpecimenId;
    }

    String getPatientId() {
        return patientId;
    }

    void setPatientId(String patientId) {
        this.patientId = patientId;
    }

    String getFamilyId() {
        return familyId;
    }

    String[] getLinkToOthersPatientId() {
        return linkToOthersPatientId;
    }

    void setLinkToOthersPatientId(String[] linkToOthersPatientId) {
        this.linkToOthersPatientId = linkToOthersPatientId;
    }


    void setFamilyId(String familyId) {
        this.familyId = familyId;
    }

    String getStudyId() {
        return studyId;
    }

    void setStudyId(String studyId) {
        this.studyId = studyId;
    }

    String getOrgId() {
        return orgId;
    }

    void setOrgId(String orgId) {
        this.orgId = orgId;
    }

    String getLabName() {
        return labName;
    }

    void setLabName(String labName) {
        this.labName = labName;
    }

    String getPractitionerId() {
        return practitionerId;
    }

    void setPractitionerId(String practitionerId) {
        this.practitionerId = practitionerId;
    }

    boolean isProban() {
        return isProban;
    }

    boolean isAffected() {
        return isAffected;
    }

    void setAffected(boolean affected) {
        isAffected = affected;
    }

    void setProban(boolean proban) {
        isProban = proban;
    }

    JSONObject getPatient() {
        return patient;
    }

    void setPatient(JSONObject patient) {
        this.patient = patient;
    }

    List<String> getHposTerms() {
        return hposTerms;
    }

    void setHposTerms(List<String> hposTerms) {
        this.hposTerms = hposTerms;
    }

    String getSequencingStrategy() {
        return sequencingStrategy;
    }

    void setSequencingStrategy(String sequencingStrategy) {
        this.sequencingStrategy = sequencingStrategy;
    }

    int getQtyHposTermsFound() {
        return qtyHposTermsFound;
    }

    void setQtyHposTermsFound(int qtyHposTermsFound) {
        this.qtyHposTermsFound = qtyHposTermsFound;
    }
    void addQtyOfHposTermsFound(int newQty) {
        this.qtyHposTermsFound += newQty;
    }

    @Override
    public String toString() {
        return "Patient{" +
                "specimenId='" + specimenId + '\'' +
                ", relation='" + relation + '\'' +
                ", patientId='" + patientId + '\'' +
                ", familyId='" + familyId + '\'' +
                ", studyId='" + studyId + '\'' +
                ", orgId='" + orgId + '\'' +
                ", labName='" + labName + '\'' +
                ", practitionerId='" + practitionerId + '\'' +
                ", isProban=" + isProban +
                ", isAffected=" + isAffected +
                ", hposTerms=" + hposTerms +
                ", qtyHposTermsFound=" + qtyHposTermsFound +
                ", sequencingStrategy='" + sequencingStrategy + '\'' +
                ", linkToOthersPatientId=" + Arrays.toString(linkToOthersPatientId) +
                ", linkToOthersSpecimenId=" + Arrays.toString(linkToOthersSpecimenId) +
                ", gender='" + gender + '\'' +
                '}';
    }
}
