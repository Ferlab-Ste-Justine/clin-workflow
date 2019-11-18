package org.chusj;

import org.json.JSONObject;

import java.io.Serializable;
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
    private JSONObject patient;
    private List<String> hposTerms;
    private String sequencingStrategy;

    Patient(String specimenId) {
        this.specimenId = specimenId;
    }


    public String getSpecimenId() {
        return specimenId;
    }

    public void setSpecimenId(String specimenId) {
        this.specimenId = specimenId;
    }

    public String getRelation() {
        return relation;
    }

    public void setRelation(String relation) {
        this.relation = relation;
    }

    public String getPatientId() {
        return patientId;
    }

    public void setPatientId(String patientId) {
        this.patientId = patientId;
    }

    public String getFamilyId() {
        return familyId;
    }

    public void setFamilyId(String familyId) {
        this.familyId = familyId;
    }

    public String getStudyId() {
        return studyId;
    }

    public void setStudyId(String studyId) {
        this.studyId = studyId;
    }

    public String getOrgId() {
        return orgId;
    }

    void setOrgId(String orgId) {
        this.orgId = orgId;
    }

    public String getLabName() {
        return labName;
    }

    public void setLabName(String labName) {
        this.labName = labName;
    }

    public String getPractitionerId() {
        return practitionerId;
    }

    public void setPractitionerId(String practitionerId) {
        this.practitionerId = practitionerId;
    }

    boolean isProban() {
        return isProban;
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

    public List<String> getHposTerms() {
        return hposTerms;
    }

    void setHposTerms(List<String> hposTerms) {
        this.hposTerms = hposTerms;
    }

    public String getSequencingStrategy() {
        return sequencingStrategy;
    }

    public void setSequencingStrategy(String sequencingStrategy) {
        this.sequencingStrategy = sequencingStrategy;
    }

    @Override
    public String toString() {
        return "Donor{" +
                "specimenId='" + specimenId + '\'' +
                ", relation='" + relation + '\'' +
                ", patientId='" + patientId + '\'' +
                ", familyId='" + familyId + '\'' +
                ", studyId='" + studyId + '\'' +
                ", orgId='" + orgId + '\'' +
                ", labName='" + labName + '\'' +
                ", practitionerId='" + practitionerId + '\'' +
                ", isProban=" + isProban +
                ", hposTerms=" + hposTerms +
                ", sequencingStrategy='" + sequencingStrategy + '\'' +
                '}';
    }
}
