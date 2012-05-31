package edu.cornell.med.icb.goby.algorithmic.data;

/**
 * Stores methylation rate information.
 * User: nyasha
 * Date: 6/14/11
 * Time: 5:10 PM
 */
public class MethylRateInfo {

    int sampleIndex;
    int groupIndex;
    int position;
    float methylationRate;

    public MethylRateInfo(int sampleIndex, int groupIndex, int position, float methylationRate) {
        this.sampleIndex = sampleIndex;
        this.groupIndex = groupIndex;
        this.position = position;
        this.methylationRate = methylationRate;
    }
}
