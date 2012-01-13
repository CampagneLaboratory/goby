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

package edu.cornell.med.icb.goby.alignments;

import edu.cornell.med.icb.goby.algorithmic.data.EquivalentIndelRegion;
import it.unimi.dsi.fastutil.objects.ObjectSet;

import java.util.Arrays;

/**
 * Provide a strategy for filtering genotype observations (bases or indels) and reduce the impact of sequencing errors
 * on downstream statistics.
 *
 * @author Fabien Campagne
 *         Date: Mar 23, 2011
 *         Time: 11:07:11 AM
 */
public abstract class GenotypeFilter {


    /**
     * Adjust genotypes and sampleCounts to remove/reduce the effect of sequencing errors.
     * Observations that are filtered by this strategy are added to filteredSet. When a previous
     * filter removed genotypes, they can be found in filteredSet.
     *
     * @param list         Variation or reference genotype observations at position
     * @param sampleCounts Counts for alleles at position each each sample under study.
     * @param filteredSet  Set of genotype observations that have been filtered
     */
    public abstract void filterGenotypes(DiscoverVariantPositionData list,
                                         SampleCountInfo[] sampleCounts,
                                         ObjectSet<PositionBaseInfo> filteredSet);

    /**
     * Returns a short description of the filtering criteria.
     *
     * @return a short description of the filtering criteria.
     */
    public String describe() {
        return this.getClass().getSimpleName();
    }

    int[] varCountRemovedPerSample;
    int[] refCountRemovedPerSample;
    int numFiltered = 0;
    int numScreened = 0;

    public double getPercentFilteredOut() {
        double rate = numFiltered;
        rate /= numScreened;
        return rate * 100d;
    }

    public String getName() {
        return "filter (" + describe() + ")";
    }

    void initStorage(int numSamples) {

        if (varCountRemovedPerSample == null) {
            varCountRemovedPerSample = new int[numSamples];
            refCountRemovedPerSample = new int[numSamples];
        } else {
            Arrays.fill(varCountRemovedPerSample, 0);
            Arrays.fill(refCountRemovedPerSample, 0);
        }
    }

    void resetCounters() {

        numScreened = 0;
        numFiltered = 0;
    }

    protected void adjustRefVarCounts(SampleCountInfo[] sampleCounts) {
        for (SampleCountInfo sci : sampleCounts) {
            sci.varCount -= varCountRemovedPerSample[sci.sampleIndex];
            sci.refCount -= refCountRemovedPerSample[sci.sampleIndex];
            // TODO remove the max statements, the code should work without
            sci.varCount = Math.max(0, sci.varCount);
            sci.refCount = Math.max(0, sci.refCount);
            assert sci.refCount >= 0 : "refCount negative: " + sci.refCount;
            assert sci.varCount >= 0 : "varCount negative: " + sci.varCount;
        }
    }

    protected void filterIndels(final DiscoverVariantPositionData list) {
        if (list.hasCandidateIndels()) {
            // remove candidate indels if they don't make the frequency threshold (threshold determined by bases observed
            // at that position):
            for (final EquivalentIndelRegion indel : list.getIndels()) {
                if (indel != null && indel.getFrequency() < getThresholdForSample(indel.sampleIndex)) {

                    list.failIndel(indel);
                }
            }
        }
    }

    public abstract int getThresholdForSample(final int sampleIndex);
}
