package org.chusj;

import java.util.ArrayList;
import java.util.List;

public class Family {

    private String familyId;
    private List<Integer> familyComposition;

    public String getFamilyId() {
        return familyId;
    }

    public void setFamilyId(String familyId) {
        this.familyId = familyId;
    }

    Family(String familyId) {
        this.familyId = familyId;
        this.familyComposition = new ArrayList<>();
    }

    List<Integer> getFamilyComposition() {
        return familyComposition;
    }

    void setFamilyComposition(List<Integer> familyComposition) {
        this.familyComposition = familyComposition;
    }

    void addFamily(Integer donorPosition) {
        this.familyComposition.add(donorPosition);
    }
}
