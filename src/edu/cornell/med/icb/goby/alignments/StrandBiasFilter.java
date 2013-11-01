package edu.cornell.med.icb.goby.alignments;

import it.unimi.dsi.fastutil.objects.ObjectSet;

/**
 * Removes genotypes when only one strand supports the genotype in a sample.
 * @author Fabien Campagne
 *         Date: 11/1/13
 *         Time: 2:24 PM
 */
public class StrandBiasFilter extends GenotypeFilter {

    /**
     * Returns a short description of the filtering criteria.
     *
     * @return a short description of the filtering criteria.
     */
    public String describe() {
        return "remove variations seen with support from a single strand.";
    }

    @Override
    public void filterGenotypes(DiscoverVariantPositionData list, SampleCountInfo[] sampleCounts, ObjectSet<PositionBaseInfo> filteredSet) {
        resetCounters();
        initStorage(sampleCounts.length);

        for (PositionBaseInfo positionBaseInfo : list) {

            numScreened++;
            final int sampleIndex = positionBaseInfo.readerIndex;
            char base = positionBaseInfo.matchesReference ? positionBaseInfo.from : positionBaseInfo.to;

            final SampleCountInfo sampleCountInfo = sampleCounts[sampleIndex];

            // how many of this base have we seen in this sample?
            final int baseIndex = sampleCountInfo.baseIndex(base);
            final int countForward = sampleCountInfo.getGenotypeCount(baseIndex,true);
            final int countReverse = sampleCountInfo.getGenotypeCount(baseIndex,false);

            if ((countForward+countReverse)>=2 && (countForward==0 || countReverse==0)) {

                // the variation is not represented on one of the strands, filter
                sampleCountInfo.suggestRemovingGenotype(baseIndex, positionBaseInfo.matchesForwardStrand);
                removeGenotype(positionBaseInfo, filteredSet);
            }
        }
        filterIndels(list, sampleCounts);
        adjustGenotypes(list, filteredSet, sampleCounts);
        adjustRefVarCounts(sampleCounts);
    }

    @Override
    public int getThresholdForSample(int sampleIndex) {
        return 0;
    }
}
