package edu.cornell.med.icb.goby.algorithmic.algorithm.dmr;

/**
 * @author Fabien Campagne
 *         Date: 5/6/12
 *         Time: 3:09 PM
 */
class WindowRange {
    /**
     * 0-based inclusive
     */
    int start;
    int end;
    int length;
    double pForRange;

    WindowRange(final int start, final int end) {
        this.start = start;
        this.end = end;
        this.length = end - start;
        pForRange = 1;
    }

    public String toString() {
        String res = "range: from ";
        res = res + start;
        res = res + " to ";
        res = res + end;
        res = res + " pvalue= ";
        res = res + pForRange;
        return res;
    }

    public void setPforRange(double pvalue) {
        pForRange = pvalue;
    }
}
