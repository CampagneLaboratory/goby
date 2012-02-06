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
    int[] methylatedCCounterPerSample;
    int[] unmethylatedCCounterPerSample;
    int[] methylatedCCounterPerGroup;
    int[] unmethylatedCcounterPerGroup;
    int[] numberOfSites;

    public double[] getMethylationRatePerSample() {
        return methylationRatePerSample;
    }

    private double[] methylationRatePerSample;
    private double[] methylationRatePerGroup;

    public int getAnnotationIndex() {
        return annotationIndex;
    }

    public FormatFieldCounter(int annotationIndex, int numSamples, int numGroups) {
        this(annotationIndex, numSamples);
        methylatedCCounterPerGroup = new int[numGroups];
        unmethylatedCcounterPerGroup = new int[numGroups];
        methylationRatePerGroup = new double[numGroups];
    }

    public FormatFieldCounter(int annotationIndex, int numSamples) {
        this.annotationIndex = annotationIndex;
        this.numSamples = numSamples;
        methylatedCCounterPerSample = new int[numSamples];
        unmethylatedCCounterPerSample = new int[numSamples];
        numberOfSites = new int[numSamples];
        methylationRatePerSample = new double[numSamples];
    }

    public void CalculateSampleMethylationRate(int sampleIndex) {
        double denominator = methylatedCCounterPerSample[sampleIndex] + unmethylatedCCounterPerSample[sampleIndex];
        double MR;
        if(denominator==0){
            MR=Double.NaN;
        }else{
        MR = ((methylatedCCounterPerSample[sampleIndex] / denominator) * 100);
        }
        methylationRatePerSample[sampleIndex] = MR;
    }

    public void CalculateGroupMethylationRate(int groupIndex) {
        double denominator = methylatedCCounterPerGroup[groupIndex] + unmethylatedCcounterPerGroup[groupIndex];
        double MR;
        if (denominator == 0) {
            MR = Double.NaN;
        } else {
            MR = ((methylatedCCounterPerGroup[groupIndex] / denominator) * 100);
        }
        methylationRatePerGroup[groupIndex] = MR;
    }

    public String toString(int sampleIndex) {
        StringBuilder result = new StringBuilder();
        result.append("C: ");
        result.append(unmethylatedCCounterPerSample[sampleIndex]);
        result.append(" Cm: ");
        result.append(methylatedCCounterPerSample[sampleIndex]);
        return result.toString();
    }

    public double[] getMethylationRatePerGroup() {
        return methylationRatePerGroup;
    }
}
