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

    private int[][] methylatedCCountPerSample;
    private int[][] unmethylatedCCountPerSample;

    private int[][] methylatedCCountPerGroup;
    private int[][] unmethylatedCCountPerGroup;


    private int[][] numberOfSitesPerSample;
    private int[][] numberOfSitesPerGroup;

    private double[][] methylationRatePerSample;
    private double[][] methylationRatePerGroup;
    private String[] contexts;


    public FormatFieldCounter(int annotationIndex, int numSamples, int numGroups, String contexts[]) {

        this(annotationIndex, numSamples, contexts);
        int numContexts = contexts.length;
        methylatedCCountPerGroup = new int[numContexts][numGroups];
        unmethylatedCCountPerGroup = new int[numContexts][numGroups];
        methylationRatePerGroup = new double[numContexts][numGroups];
        numberOfSitesPerGroup = new int[numContexts][numGroups];

    }

    public FormatFieldCounter(int annotationIndex, int numSamples, String[] contexts) {
        this.annotationIndex = annotationIndex;
        this.numSamples = numSamples;
        this.contexts = contexts;
        int numContexts = contexts.length;

        methylatedCCountPerSample = new int[numContexts][numSamples];
        unmethylatedCCountPerSample = new int[numContexts][numSamples];
        numberOfSitesPerSample = new int[numContexts][numSamples];
        methylationRatePerSample = new double[numContexts][numSamples];
    }

    public void calculateSampleMethylationRate(int contextIndex, int sampleIndex) {
        final double denominator = methylatedCCountPerSample[contextIndex][sampleIndex] +
                unmethylatedCCountPerSample[contextIndex][sampleIndex];
        final double MR;
        if (denominator == 0) {
            MR = Double.NaN;
        } else {
            MR = (((double) methylatedCCountPerSample[contextIndex][sampleIndex]) / denominator) * 100.0;
        }
        methylationRatePerSample[contextIndex][sampleIndex] = MR;
    }

    public void calculateGroupMethylationRate(int contextIndex, int groupIndex) {
        double denominator = methylatedCCountPerGroup[contextIndex][groupIndex] +
                unmethylatedCCountPerGroup[contextIndex][groupIndex];
        double MR;
        if (denominator == 0) {
            MR = Double.NaN;
        } else {
            MR = ((((double) methylatedCCountPerGroup[contextIndex][groupIndex]) / denominator) * 100.0);
        }
        methylationRatePerGroup[contextIndex][groupIndex] = MR;
    }

    public String toString(int contextIndex, int sampleIndex) {
        StringBuilder result = new StringBuilder();
        result.append("C: ");
        result.append(unmethylatedCCountPerSample[contextIndex][sampleIndex]);
        result.append(" Cm: ");
        result.append(methylatedCCountPerSample[contextIndex][sampleIndex]);
        return result.toString();
    }


    public boolean countsAreNotNull(int c, int cm) {
        boolean flag = false;
        double doubleC = (double) c;
        double doubleCm = (double) cm;
        if (doubleC != doubleC && doubleCm != doubleCm) {
            // value is NaN
        } else {
            flag = true;
        }
        return flag;
    }


    public void incrementCounts(int sampleIndex, int[] sampleIndexToGroupIndex, int c,
                                int cm, int contextIndex) {
        // only increment counters if any cytosines are observed to pass the threshold
        if (countsAreNotNull(c, cm)) {
            unmethylatedCCountPerSample[contextIndex][sampleIndex] += c;
            methylatedCCountPerSample[contextIndex][sampleIndex] += cm;
            numberOfSitesPerSample[contextIndex][sampleIndex] += 1;
            if (sampleIndexToGroupIndex != null) {
                final int groupIndex = sampleIndexToGroupIndex[sampleIndex];
                // System.out.printf("increment: context=%s sampleIndex=%d groupIndex=%d c=%d cm=%d %n", contexts[contextIndex], sampleIndex, groupIndex, c, cm);
                unmethylatedCCountPerGroup[contextIndex][groupIndex] += c;
                methylatedCCountPerGroup[contextIndex][groupIndex] += cm;
                numberOfSitesPerGroup[contextIndex][groupIndex] += 1;
            }
        }
    }

    public int getNumberOfSitesPerSample(int contextIndex, int sampleIndex) {
        return numberOfSitesPerSample[contextIndex][sampleIndex];
    }

    public int getNumberOfSitesPerGroup(int contextIndex, int indexGroup) {
        return numberOfSitesPerGroup[contextIndex][indexGroup];
    }

    public int getUnmethylatedCCountPerSample(int contextIndex, int sampleIndex) {
        return unmethylatedCCountPerSample[contextIndex][sampleIndex];
    }

    public int getUnmethylatedCcountPerGroup(int contextIndex, int indexGroup) {
        return unmethylatedCCountPerGroup[contextIndex][indexGroup];
    }

    public int getMethylatedCCountPerGroup(int contextIndex, int indexGroup) {
        return methylatedCCountPerGroup[contextIndex][indexGroup];
    }

    public int getMethylatedCCountPerSample(int contextIndex, int sampleIndex) {
        return methylatedCCountPerSample[contextIndex][sampleIndex];
    }

    public double getMethylationRatePerSample(int contextIndex, int sampleIndex) {
        calculateSampleMethylationRate(contextIndex, sampleIndex);
        return methylationRatePerSample[contextIndex][sampleIndex];
    }

    public double getMethylationRatePerGroup(int contextIndex, int groupIndex) {
        calculateGroupMethylationRate(contextIndex, groupIndex);
        return methylationRatePerGroup[contextIndex][groupIndex];
    }

    public int getAnnotationIndex() {
        return annotationIndex;
    }

    /**
     * Return the rate in a sample, irrespective of contexts.
     *
     * @param sampleIndex index of the sample.
     * @return methylation rate estimate in the group.
     */
    public double getMethylationRatePerSample(final int sampleIndex) {
        double MR;
        double numerator = 0;
        double denominator = 0;
        for (int contextIndex = 0; contextIndex < methylatedCCountPerSample.length; contextIndex++) {
            denominator += methylatedCCountPerSample[contextIndex][sampleIndex] +
                    unmethylatedCCountPerSample[contextIndex][sampleIndex];
            numerator += methylatedCCountPerSample[contextIndex][sampleIndex];
        }
        if (denominator == 0) {
            MR = Double.NaN;
        } else {
            MR = numerator / denominator * 100.0;
        }
        return MR;
    }

    /**
     * Return the rate in a group, irrespective of contexts.
     *
     * @param groupIndex index of the group.
     * @return methylation rate estimate in the group.
     */
    public double getMethylationRatePerGroup(final int groupIndex) {
        double MR;
        double numerator = 0;
        double denominator = 0;
        for (int contextIndex = 0; contextIndex < methylatedCCountPerGroup.length; contextIndex++) {
            denominator += methylatedCCountPerGroup[contextIndex][groupIndex] +
                    unmethylatedCCountPerGroup[contextIndex][groupIndex];
            numerator += methylatedCCountPerGroup[contextIndex][groupIndex];
        }
        if (denominator == 0) {
            MR = Double.NaN;
        } else {
            MR = numerator / denominator * 100.0;
        }
        return MR;
    }


}
