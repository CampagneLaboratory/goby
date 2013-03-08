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
 * This filter considers whether the remaining base calls of each allele (left-over)
 * have a count larger than the number of call corrections done by previous filters.
 * If this is not the case, the allele is considered mis-called and its observations
 * filtered out. This is a simple and efficient base call correction strategy that
 * improves agreement between genotypes in technical replicates.
 *
 * @author Fabien Campagne
 *         Date: Mar 24, 2011
 *         Time: 11:42:42 AM
 */
public class LeftOverFilter extends GenotypeFilter {
    private static final int MULTIPLIER = 2;
    private int minVariationSupport = 0;

    public LeftOverFilter(int minVariationSupport) {
        this.minVariationSupport = minVariationSupport;
    }

    @Override
    void initStorage(int numSamples) {
        super.initStorage(numSamples);
        if (thresholdsPerSample == null) {
            thresholdsPerSample = new int[numSamples];
        }
        Arrays.fill(thresholdsPerSample, minVariationSupport);
    }

    private int[] thresholdsPerSample;


    @Override
    public void filterGenotypes(DiscoverVariantPositionData list,
                                SampleCountInfo[] sampleCounts,
                                ObjectSet<PositionBaseInfo> filteredList) {
        resetCounters();
        initStorage(sampleCounts.length);
        for (SampleCountInfo sci : sampleCounts) {
            sci.clearFiltered();
        }
        for (PositionBaseInfo positionBaseInfo : filteredList) {
            thresholdsPerSample[positionBaseInfo.readerIndex] += 1;
        }
        if (list.hasCandidateIndels()) {
            for (EquivalentIndelRegion indel : list.getIndels()) {
                if (indel.isFiltered()) {
                    thresholdsPerSample[indel.sampleIndex] += 1;
                }
            }
        }
        for (int sampleIndex = 0; sampleIndex < sampleCounts.length; sampleIndex++) {
            thresholdsPerSample[sampleIndex] *= MULTIPLIER;
        }

        for (PositionBaseInfo positionBaseInfo : list) {

            numScreened++;
            final int sampleIndex = positionBaseInfo.readerIndex;
            char base = positionBaseInfo.matchesReference ? positionBaseInfo.from : positionBaseInfo.to;

            final SampleCountInfo sampleCountInfo = sampleCounts[sampleIndex];
            // how many of this base have we seen in this sample?
            final int baseIndex = sampleCountInfo.baseIndex(base);
            final int count = sampleCountInfo.counts[baseIndex];
            if (count == 0) continue;
            if (count < thresholdsPerSample[sampleIndex]) {

                // less counts remaining for this allele than were removed on average by previous filters, still likely
                // an error.
                // We remove this call

                sampleCountInfo.suggestRemovingGenotype(baseIndex);
                removeGenotype(positionBaseInfo, filteredList);
            }
        }

        filterIndels(list, sampleCounts);
        adjustGenotypes(list, filteredList, sampleCounts);
        adjustRefVarCounts(sampleCounts);
    }


    @Override
    public String describe() {
        return String.format("#count(allele) < (%d *#filtered)", MULTIPLIER);
    }

    @Override
    public int getThresholdForSample(int sampleIndex) {
        return thresholdsPerSample[sampleIndex];
    }

}
