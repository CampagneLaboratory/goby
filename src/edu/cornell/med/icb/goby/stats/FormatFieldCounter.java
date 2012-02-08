/*
 * Copyright (C) 2009-2012 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.cornell.med.icb.goby.stats;

/**
 * Keeps track of different statistics or counts for chosen format fields
 *
 * @author Nyasha Chambwe
 *         Date: 1/27/12
 *         Time: 10:17 AM
 */
public class FormatFieldCounter {

    private int annotationIndex;
    private int numSamples;

    private int[][] methylatedCCounterPerSample;
    private int[][] unmethylatedCCounterPerSample;

    private int[][] methylatedCCounterPerGroup;
    private int[][] unmethylatedCCounterPerGroup;

    private int[][] numberOfSites;

    private double[][] methylationRatePerSample;
    private double[][] methylationRatePerGroup;


    public FormatFieldCounter(int annotationIndex, int numSamples, int numGroups, int numContexts) {
        this(annotationIndex, numSamples, numContexts);
        methylatedCCounterPerGroup = new int[numContexts][numGroups];
        unmethylatedCCounterPerGroup = new int[numContexts][numGroups];
        methylationRatePerGroup = new double[numContexts][numGroups];
    }

    public FormatFieldCounter(int annotationIndex, int numSamples, int numContexts) {
        this.annotationIndex = annotationIndex;
        this.numSamples = numSamples;
        methylatedCCounterPerSample = new int[numContexts][numSamples];
        unmethylatedCCounterPerSample = new int[numContexts][numSamples];
        numberOfSites = new int[numContexts][numSamples];
        methylationRatePerSample = new double[numContexts][numSamples];
    }

    public void CalculateSampleMethylationRate(int sampleIndex, int contextIndex) {
        double denominator = methylatedCCounterPerSample[contextIndex][sampleIndex] +
                unmethylatedCCounterPerSample[contextIndex][sampleIndex];
        double MR;
        if (denominator == 0) {
            MR = Double.NaN;
        } else {
            MR = ((methylatedCCounterPerSample[contextIndex][sampleIndex] / denominator) * 100);
        }
        methylationRatePerSample[contextIndex][sampleIndex] = MR;
    }

    public void CalculateGroupMethylationRate(int groupIndex, int contextIndex) {
        double denominator = methylatedCCounterPerGroup[contextIndex][groupIndex] +
                unmethylatedCCounterPerGroup[contextIndex][groupIndex];
        double MR;
        if (denominator == 0) {
            MR = Double.NaN;
        } else {
            MR = ((methylatedCCounterPerGroup[contextIndex][groupIndex] / denominator) * 100);
        }
        methylationRatePerGroup[contextIndex][groupIndex] = MR;
    }

    public String toString(int contextIndex, int sampleIndex) {
        StringBuilder result = new StringBuilder();
        result.append("C: ");
        result.append(unmethylatedCCounterPerSample[contextIndex][sampleIndex]);
        result.append(" Cm: ");
        result.append(methylatedCCounterPerSample[contextIndex][sampleIndex]);
        return result.toString();
    }


    public void incrementCounts(int sampleIndex, int[] sampleIndexToGroupIndex, int c,
                                int cm, int contextIndex, boolean processGroups) {
        unmethylatedCCounterPerSample[contextIndex][sampleIndex] += c;
        methylatedCCounterPerSample[contextIndex][sampleIndex] += cm;
        if(processGroups){
        unmethylatedCCounterPerGroup[contextIndex][sampleIndexToGroupIndex[sampleIndex]] += c;
        methylatedCCounterPerGroup[contextIndex][sampleIndexToGroupIndex[sampleIndex]] += cm;
        }
        numberOfSites[contextIndex][sampleIndex] += 1;
    }

    public int[] getUnmethylatedCcounterPerGroup(int contextIndex) {
        return unmethylatedCCounterPerGroup[contextIndex];
    }

    public int[] getMethylatedCCounterPerGroup(int contextIndex) {
        return methylatedCCounterPerGroup[contextIndex];
    }


    public double[] getMethylationRatePerSample(int contextIndex) {
        return methylationRatePerSample[contextIndex];
    }

    public double[] getMethylationRatePerGroup(int contextIndex) {
        return methylationRatePerGroup[contextIndex];
    }

    public int getAnnotationIndex() {
        return annotationIndex;
    }

}
