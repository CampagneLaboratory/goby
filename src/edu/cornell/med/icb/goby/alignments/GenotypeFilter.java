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
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

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
     * Used to log debug and informational messages.
     */
    protected static final Log LOG = LogFactory.getLog(DiscoverVariantIterateSortedAlignments.class);

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
    /**
     * List of sampleIndex,genotypeIndex tuples. Items in this list indicate genotypes suggested for removal by some
     * filter.
     */
    protected LongArrayList removalSuggestions = new LongArrayList();

    public double getPercentFilteredOut() {
        double rate = numFiltered;
        rate /= numScreened;
        return rate * 100d;
    }

    public String getName() {
        return "filter (" + describe() + ")";
    }

    public void initStorage(int numSamples) {

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

    private ObjectArraySet<EquivalentIndelRegion> toBeRemoved = new ObjectArraySet<EquivalentIndelRegion>();
    private IntArraySet doNotRemoveSampleIndices = new IntArraySet();

    protected void filterIndels(final DiscoverVariantPositionData list, SampleCountInfo[] sampleCounts) {
        if (list.hasCandidateIndels()) {
            // remove candidate indels if they don't make the frequency threshold (threshold determined by bases observed
            // at that position):
            for (final EquivalentIndelRegion indel : list.getIndels()) {
                if (indel != null && indel.getFrequency() < getThresholdForSample(indel.sampleIndex)) {
                    indel.markFiltered();
                    list.failIndel(indel);
                    sampleCounts[indel.sampleIndex].removeIndel(indel);
                }
            }
        }
    }

    /**
     * Adjust the genotypes considering genotype removals and thresholding effects.
     *
     * @param list
     * @param filteredList
     * @param sampleCounts
     */
    protected void adjustGenotypes(DiscoverVariantPositionData list, ObjectSet<PositionBaseInfo> filteredList,
                                   SampleCountInfo[] sampleCounts) {


    }

    /**
     * Use this method to remove a genotype in a sub-class filter.
     *
     * @param info         the base to remove from consideration, according to the filter logic.
     * @param filteredList the list of filtered bases to add to.
     */
    protected void removeGenotype(PositionBaseInfo info, ObjectSet<PositionBaseInfo> filteredList) {
        filteredList.add(info);

        final int sampleIndex = info.readerIndex;
        if (info.matchesReference) {
            refCountRemovedPerSample[sampleIndex]++;
        } else {
            varCountRemovedPerSample[sampleIndex]++;
        }
        numFiltered++;
    }


    public abstract int getThresholdForSample(final int sampleIndex);
}
