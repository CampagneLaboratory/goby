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

package edu.cornell.med.icb.goby.alignments;

import it.unimi.dsi.fastutil.objects.ObjectSet;

import java.util.Arrays;

/**
 * Filters called genotypes to keep the two highest frequency genotypes. This filter is appropriate for diploid genomes
 * as it will force the calls to have at most two alleles.
 *
 * @author Fabien Campagne
 *         Date: 2/7/12
 *         Time: 1:35 PM
 */
public class DiploidFilter extends GenotypeFilter {
    // indexed per sample:
    private int firstMaxFrequency[];
    private int secondMaxFrequency[];

    @Override
    void initStorage(final int numSamples) {
        super.initStorage(numSamples);
        if (firstMaxFrequency == null) {
            firstMaxFrequency = new int[numSamples];
            secondMaxFrequency = new int[numSamples];
        } else {
            Arrays.fill(firstMaxFrequency, 0);
            Arrays.fill(secondMaxFrequency, 0);
        }
    }

    @Override
    public void filterGenotypes(final DiscoverVariantPositionData list,
                                final SampleCountInfo[] sampleCounts,
                                final ObjectSet<PositionBaseInfo> filteredSet) {
        resetCounters();
        initStorage(sampleCounts.length);
        for (final SampleCountInfo sci : sampleCounts) {
            final int sampleIndex = sci.sampleIndex;
            int maxIndex = sci.getGenotypeMaxIndex();

            for (int genotypeIndex = 0; genotypeIndex < maxIndex; genotypeIndex++) {
                final int count = sci.getGenotypeCount(genotypeIndex);
                firstMaxFrequency[sampleIndex] = Math.max(firstMaxFrequency[sampleIndex], count);
            }
            for (int genotypeIndex = 0; genotypeIndex < maxIndex; genotypeIndex++) {
                final int count = sci.getGenotypeCount(genotypeIndex);
                if (count < firstMaxFrequency[sampleIndex]) {
                    secondMaxFrequency[sampleIndex] = Math.max(secondMaxFrequency[sampleIndex], count);
                }
            }
            int genotypesWithFirstOrSecondMax = 0;
            for (int genotypeIndex = 0; genotypeIndex < maxIndex; genotypeIndex++) {
                final int count = sci.getGenotypeCount(genotypeIndex);
                if (count == firstMaxFrequency[sampleIndex] || count == secondMaxFrequency[sampleIndex]) {
                    genotypesWithFirstOrSecondMax++;
                }
            }
            if (genotypesWithFirstOrSecondMax > 2) {
                // more than two genotypes with top or second best frequency. Keep only the best frequency:
                secondMaxFrequency[sampleIndex] = -1;  // -1 will disable matching the genotpe to secondMaxFrequency
            }
        }

        for (final PositionBaseInfo positionBaseInfo : list) {

            numScreened++;
            final int sampleIndex = positionBaseInfo.readerIndex;
            char base = positionBaseInfo.matchesReference ? positionBaseInfo.from : positionBaseInfo.to;

            final SampleCountInfo sampleCountInfo = sampleCounts[sampleIndex];
            // how many of this base have we seen in this sample?
            final int baseIndex = sampleCountInfo.baseIndex(base);
            final int count = sampleCountInfo.counts[baseIndex];
            if (count == 0) {
                continue;
            }
            if (count != firstMaxFrequency[sampleIndex] && count != secondMaxFrequency[sampleIndex]) {
                // remove this genotype.
                if (!filteredSet.contains(positionBaseInfo)) {
                    sampleCountInfo.counts[baseIndex]--;
                    if (positionBaseInfo.matchesReference) {
                        refCountRemovedPerSample[sampleIndex]++;
                    } else {

                        varCountRemovedPerSample[sampleIndex]++;
                    }
                    filteredSet.add(positionBaseInfo);
                    numFiltered++;
                }

            }
        }

        filterIndels(list, sampleCounts);
        adjustRefVarCounts(sampleCounts);
    }

    @Override
    public int getThresholdForSample(final int sampleIndex) {
        if (secondMaxFrequency[sampleIndex] == -1) {
            return firstMaxFrequency[sampleIndex];
        } else {
            return secondMaxFrequency[sampleIndex];
        }
    }
}
