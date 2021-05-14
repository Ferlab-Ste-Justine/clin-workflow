package org.chusj;

public class Frequencies {

    float ac =0;
    float pn =0;
    float an = 0;
    float hc = 0;
    float af =0f;

    public float getAc() {
        return ac;
    }

    public void setAc(float ac) {
        this.ac = ac;
    }

    public float getPn() {
        return pn;
    }

    public void setPn(float pn) {
        this.pn = pn;
    }

    public float getAn() {
        return an;
    }

    public void setAn(float an) {
        this.an = an;
    }

    public float getHc() {
        return hc;
    }

    public void setHc(float hc) {
        this.hc = hc;
    }

    public float getAf() {
        return af;
    }

    public void setAf(float af) {
        this.af = af;
    }

    public void updateAf() {
        if (this.an > 0) {
            this.af = this.ac / this.an;
        }
    }

    @Override
    public String toString() {
        return "Frequencies{" +
                "ac=" + ac +
                ", pn=" + pn +
                ", an=" + an +
                ", hc=" + hc +
                ", af=" + af +
                '}';
    }
}
