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
    int[] methylatedCounter;
    int[] unmethylatedCounter;
    int[] numberOfSites;

    public double[] getMethylationRate() {
        return methylationRate;
    }

    private double [] methylationRate;

    public int getAnnotationIndex() {
        return annotationIndex;
    }


    public FormatFieldCounter(int annotationIndex, int numSamples) {
        this.annotationIndex = annotationIndex;
        this.numSamples = numSamples;
        methylatedCounter= new int[numSamples];
        unmethylatedCounter= new int[numSamples];
        numberOfSites=new int[numSamples];
        methylationRate=new double[numSamples];
    }

    public void CalculateSampleMethylationRate(int sampleIndex){
        double denominator= methylatedCounter[sampleIndex] + unmethylatedCounter[sampleIndex];
        double MR= ((methylatedCounter[sampleIndex]/denominator)*100);
        methylationRate[sampleIndex]=MR;
    }

    public String toString(int sampleIndex) {
        StringBuilder result = new StringBuilder();
        result.append("C: ");
        result.append(unmethylatedCounter[sampleIndex]);
        result.append(" Cm: ");
        result.append(methylatedCounter[sampleIndex]);
        return result.toString();
    }
}
