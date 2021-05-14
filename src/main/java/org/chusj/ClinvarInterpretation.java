package org.chusj;

import java.io.Serializable;
import java.util.Objects;

public class ClinvarInterpretation implements Serializable {

    String interpretation;
    String condition;


    public String getInterpretation() {
        return interpretation;
    }

    public void setInterpretation(String interpretation) {
        this.interpretation = interpretation;
    }

    public String getCondition() {
        return condition;
    }

    public void setCondition(String condition) {
        this.condition = condition;
    }


    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        ClinvarInterpretation that = (ClinvarInterpretation) o;
        return Objects.equals(interpretation, that.interpretation) &&
                Objects.equals(condition, that.condition);
    }

    @Override
    public int hashCode() {
        return Objects.hash(interpretation, condition);
    }

    @Override
    public String toString() {
        return "ClinvarInterpretation{" +
                "interpretation='" + interpretation + '\'' +
                ", condition='" + condition + '\'' +
                '}';
    }
}
