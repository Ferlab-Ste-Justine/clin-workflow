package org.chusj;

import org.json.JSONArray;
import org.json.JSONObject;

import java.io.Serializable;
import java.util.*;

public class Clinvar implements Serializable {

    String id;
    Set<ClinvarInterpretation> interpretationSet = new HashSet<>();
    List<String> omimIds = new ArrayList<>();


    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }


    public List<String> getOmimIds() {
        return omimIds;
    }

    public void setOmimIds(List<String> omimIds) {
        this.omimIds = omimIds;
    }

    public void addOmim(String id) {
        this.omimIds.add(id);
    }

    public Set<ClinvarInterpretation> getInterpretationSet() {
        return interpretationSet;
    }

    public void setInterpretationSet(Set<ClinvarInterpretation> interpretationSet) {
        this.interpretationSet = interpretationSet;
    }

    public void addInterpretation(ClinvarInterpretation interpretation) {
        this.interpretationSet.add(interpretation);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Clinvar clinvar = (Clinvar) o;
        return Objects.equals(id, clinvar.id) &&
                Objects.equals(interpretationSet, clinvar.interpretationSet) &&
                Objects.equals(omimIds, clinvar.omimIds);
    }

    @Override
    public int hashCode() {
        return Objects.hash(id, interpretationSet, omimIds);
    }

    @Override
    public String toString() {
        return "Clinvar{" +
                "id='" + id + '\'' +
                ", interpretationSet=" + interpretationSet +
                ", omimIds=" + omimIds +
                '}';
    }

    JSONObject getJsonObj() {
        JSONObject clinvar = new JSONObject();
        JSONArray clinvarInt = new JSONArray();
        JSONArray omimIdArray = new JSONArray();

        clinvar.put("id", this.id);

        for (String omimId : this.omimIds) {
            omimIdArray.put(omimId);
        }
        if (omimIdArray.length() > 0) {
            clinvar.put("omim", omimIdArray);
        }
        for (ClinvarInterpretation ci : this.interpretationSet) {
            JSONObject ciObj = new JSONObject();
            ciObj.put("interpretation", ci.getInterpretation());
            //ciObj.put()
        }


        return clinvar;
    }
}
