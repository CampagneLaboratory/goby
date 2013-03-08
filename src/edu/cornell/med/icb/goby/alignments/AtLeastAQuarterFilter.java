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

import it.unimi.dsi.fastutil.objects.ObjectSet;

import java.util.Arrays;

/**
 * Filter that rejects alleles if their count is not at least a quarter of the count of the allele with the most counts.
 * This filter is useful when the samples are from unpooled diploid genomes. It will remove the occasional sequencing
 * error that is not flagged by low base quality score.
 *
 * @author Fabien Campagne
 *         Date: Apr 23, 2011
 *         Time: 1:57:30 PM
 */
public class AtLeastAQuarterFilter extends GenotypeFilter {
    private int[] maxAlleleCountsPerSample;

    @Override
    public String describe() {
        return "count(allele in sample) < 1/4 * max_count over alleles in sample";
    }

    void initStorage(int numSamples) {
        super.initStorage(numSamples);
        if (maxAlleleCountsPerSample == null) {
            maxAlleleCountsPerSample = new int[numSamples];
        } else {
            Arrays.fill(maxAlleleCountsPerSample, 0);
        }
    }

    @Override
    public int getThresholdForSample(final int sampleIndex) {
        return maxAlleleCountsPerSample[sampleIndex] / 4;
    }

    public void filterGenotypes(DiscoverVariantPositionData list,
                                SampleCountInfo[] sampleCounts,
                                ObjectSet<PositionBaseInfo> filteredList) {

        resetCounters();
        initStorage(sampleCounts.length);

        for (SampleCountInfo sci : sampleCounts) {
            sci.clearFiltered();
            for (int genotypeIndex = 0; genotypeIndex < sci.getGenotypeMaxIndex(); ++genotypeIndex) {
                final int count = sci.getGenotypeCount(genotypeIndex);
                maxAlleleCountsPerSample[sci.sampleIndex] = Math.max(maxAlleleCountsPerSample[sci.sampleIndex], count);
            }
        }

        for (PositionBaseInfo positionBaseInfo : list) {

            final int sampleIndex = positionBaseInfo.readerIndex;

            final int removedBaseCountThreshold = getThresholdForSample(sampleIndex);

            numScreened++;
            char base = positionBaseInfo.matchesReference ? positionBaseInfo.from : positionBaseInfo.to;

            final SampleCountInfo sampleCountInfo = sampleCounts[sampleIndex];
            // how many of this base have we seen in this sample?
            final int baseIndex = sampleCountInfo.baseIndex(base);
            final int count = sampleCountInfo.counts[baseIndex];
            if (count == 0) continue;
            if (count < removedBaseCountThreshold) {

                // this allele has less than 1/4 of the counts of the allele with the most counts in this sample.
                // remove.
                sampleCountInfo.suggestRemovingGenotype(baseIndex);
                removeGenotype(positionBaseInfo, filteredList);
            }
        }
        filterIndels(list, sampleCounts);
        adjustGenotypes(list, filteredList, sampleCounts);
        adjustRefVarCounts(sampleCounts);
    }

}
