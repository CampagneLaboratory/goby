/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.xml;

import javax.xml.bind.annotation.XmlRootElement;
import java.util.Arrays;

/**
 * Stores statistics about CpG sites in a genome and matching RRBS or MethylSeq sample.
 *
 * @author Fabien Campagne
 *         Date: 9/18/11
 *         Time: 1:14 PM
 */
@XmlRootElement
public class MethylStats {

    /**
     * Number of cytosine in non-CpG context that were found converted.
     */
    public long numConvertedNotCpGContext;

    /**
     * Number of cytosine in non-CpG context that were observed in any state.
     */
    public long numNotCpGContext;

    public MethylStats(final int[] depths, final int[] fragmentLengthBins) {
        this.depths = depths;
        this.fragmentLengthBins = fragmentLengthBins;
        numberCpGsPerDepth = new long[depths.length];
        numberCpGsPerFragmentBinGenome = new long[fragmentLengthBins.length];
        numberCpGsPerFragmentBinObserved = new long[fragmentLengthBins.length];
    }

    /**
     * The sample for which these stats were collected.
     */
    public String sampleId;
    /**
     * The number of CpGs in the genome.
     */
    long numberCpGsInGenome;

    /**
     * No arg constructor required by JAXB.
     */
    public MethylStats() {
    }

    public long getNumberCpGsInGenome() {
        return numberCpGsInGenome;
    }

    public long getNumberCpGsObserved() {
        return numberCpGsObserved;
    }

    public int[] getDepths() {
        return depths;
    }

    public long[] getNumberCpGsPerDepth() {
        return numberCpGsPerDepth;
    }

    public int[] getFragmentLengthBins() {
        return fragmentLengthBins;
    }

    public long[] getNumberCpGsPerFragmentBinGenome() {
        return numberCpGsPerFragmentBinGenome;
    }

    public long[] getNumberCpGsPerFragmentBinObserved() {
        return numberCpGsPerFragmentBinObserved;
    }

    /**
     * The number of CpGs observed (i.e., with depth >=1) in the sample.
     */
    long numberCpGsObserved;
    /**
     * The depths for which statistics are collected.
     */
    int[] depths;
    /**
     * The number of CpGs observed with depth >=d
     */
    long[] numberCpGsPerDepth;
    /**
     * The length that separates a CpG from the next CpG. This encodes a range,
     * with fragmentLengthBins[i] representing fragments of length [ fragmentLengthBins[i] - fragmentLengthBins[i+1]  [
     * Indices larger than fragmentLengthBins.length should be considered to represent infinite length.
     */
    int[] fragmentLengthBins;
    /**
     * The number of CpGs (numberCpGsPerFragmentBin[i]) that start a fragment of length [ fragmentLengthBins[i] - fragmentLengthBins[i+1]  [
     */
    long[] numberCpGsPerFragmentBinGenome;
    /**
     * The number of CpGs (numberCpGsPerFragmentBin[i]) that start a fragment of length [ fragmentLengthBins[i] - fragmentLengthBins[i+1]  [
     */
    long[] numberCpGsPerFragmentBinObserved;

    public void genomeHasCpG(final int fragmentLength) {
        numberCpGsInGenome += 1;
        numberCpGsPerFragmentBinGenome[indexToIncrement(fragmentLengthBins, fragmentLength)] += 1;
    }

    private int indexToIncrement(final int[] array, final int depth) {
        final int r = Arrays.binarySearch(array, depth);
        final int indexToIncrement = r >= 0 ? r : -(r + 1);

        if (indexToIncrement >= array.length) {
            return array.length - 1;
        } else {
            assert indexToIncrement >= 0 && indexToIncrement < array.length :
                    String.format("index must be within array bounds: indexToIncrement: %d r: %d %n", indexToIncrement, r);
            return indexToIncrement;

        }
    }

    public void setNumberCpGsInGenome(long numberCpGsInGenome) {
        this.numberCpGsInGenome = numberCpGsInGenome;
    }

    public void setNumberCpGsObserved(long numberCpGsObserved) {
        this.numberCpGsObserved = numberCpGsObserved;
    }

    public void setDepths(int[] depths) {
        this.depths = depths;
    }

    public void setNumberCpGsPerDepth(long[] numberCpGsPerDepth) {
        this.numberCpGsPerDepth = numberCpGsPerDepth;
    }

    public void setFragmentLengthBins(int[] fragmentLengthBins) {
        this.fragmentLengthBins = fragmentLengthBins;
    }

    public void setNumberCpGsPerFragmentBinGenome(long[] numberCpGsPerFragmentBinGenome) {
        this.numberCpGsPerFragmentBinGenome = numberCpGsPerFragmentBinGenome;
    }

    public void setNumberCpGsPerFragmentBinObserved(long[] numberCpGsPerFragmentBinObserved) {
        this.numberCpGsPerFragmentBinObserved = numberCpGsPerFragmentBinObserved;
    }

    /**
     * This method should be called when a CpG has been observed in a sample.
     *
     * @param depth          Number of passing filter bases that covered the CpG site.
     * @param fragmentLength Length of the fragment between the observed CpG and the next CpG in the genome.
     */
    public void observedInSample(final int depth, final int fragmentLength) {
        if (depth >= 1) {
            numberCpGsPerFragmentBinObserved[indexToIncrement(fragmentLengthBins, fragmentLength)] += 1;
            numberCpGsPerDepth[indexToIncrement(depths, depth)] += 1;
            numberCpGsObserved++;
        }
    }

    /**
     * Generate a copy of this object with the same genome stats as the original.
     *
     * @return a copy.
     */
    public MethylStats copy() {
        MethylStats copy = new MethylStats(depths, fragmentLengthBins);
        System.arraycopy(numberCpGsPerFragmentBinGenome, 0, copy.numberCpGsPerFragmentBinGenome, 0, numberCpGsPerFragmentBinGenome.length);
        copy.numberCpGsInGenome = numberCpGsInGenome;
        return copy;
    }

    /**
     * Frequency of methylated C followed by G A C or T (in array order).
     */
    final long[] mcpXFrequencies = new long[4];
        public static final int CPMIN = 0;
        public static final int CPG = 0;
        public static final int CPA = 1;
        public static final int CPC = 2;
        public static final int CPT = 3;
        public static final int CPMAX = 4;

    public long[] getMethylCpXFreqs() {
        return mcpXFrequencies;
    }

    /**
     * The number of cytosines in CpG or TpG in reference CpG context seen in the sample.
     */
    public long numCTpG;

}
