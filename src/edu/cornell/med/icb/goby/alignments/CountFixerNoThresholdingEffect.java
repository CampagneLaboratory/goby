package edu.cornell.med.icb.goby.alignments;

import edu.cornell.med.icb.goby.algorithmic.data.EquivalentIndelRegion;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.ObjectSet;

/**
 * @author Fabien Campagne
 *         Date: 3/15/13
 *         Time: 11:21 AM
 */
public class CountFixerNoThresholdingEffect extends CountFixer {
    /**
     * This implementation of CountFixer prevents thresholding effect that can happen when a genotype is filtered
     * in a sample, but not in others because the frequency is just below some filter threshold. To this end, the
     * fixer only removes genotypes when the filters removed the genotype from all samples under consideration.
     *
     * @param list
     * @param sampleCounts
     * @param likelyErrors
     */
    @Override
    public void fix(final DiscoverVariantPositionData list,
                    final SampleCountInfo[] sampleCounts,
                    final ObjectSet<PositionBaseInfo> likelyErrors) {
        assert beforeFilterCounts != null : "The preserveCounts method must be called before any filter was applied. beforeFilterCounts cannot be null.";
        int maxGenotypeIndex = 0;
        for (SampleCountInfo sci : sampleCounts) {
            maxGenotypeIndex = Math.max(maxGenotypeIndex, sci.getGenotypeMaxIndex());
        }
        for (int genotypeIndex = 0; genotypeIndex < maxGenotypeIndex; genotypeIndex++) {
            // determine if the genotype was present before filtering:
            boolean genotypePresentInSomeSample = false;
            boolean genotypeRemainsInSomeSample = false;
            for (SampleCountInfo sci : sampleCounts) {

                genotypePresentInSomeSample |= beforeFilterCounts[sci.sampleIndex].getInt(genotypeIndex) > 0;
                genotypeRemainsInSomeSample |= sci.getGenotypeCount(genotypeIndex) > 0;
            }
            // when the genotype was present, and it remains after filtering in some sample, we revert the effect of the
            // filters in all samples:
            if (genotypePresentInSomeSample && genotypeRemainsInSomeSample) {
                // we need to revert filtering for this genotype to prevent artifactual thresholding effects.
                for (SampleCountInfo sci : sampleCounts) {
                    final int countBeforeFiltering = beforeFilterCounts[sci.sampleIndex].getInt(genotypeIndex);
                    sci.setGenotypeCount(genotypeIndex, countBeforeFiltering);
                    if (sci.isIndel(genotypeIndex)) {
                   //     System.out.println("reverting indel genotype.");
                        EquivalentIndelRegion indel = sci.getIndelGenotype(genotypeIndex);
                        indel.removeFiltered();
                        list.getFailedIndels().remove(indel);
                    }
                }
            }
        }
        // do not decrement counts again. The Filters have done this already..

        for (SampleCountInfo sci : sampleCounts) {
            for (int i = 0; i < sci.counts.length; i++) {
                assert sci.counts[i] >= 0 : ("Counts must never be negative. This would happen if a GenotypeFilter removed counts directly. Value was " + sci.counts[i]);
            }
        }
        // calculate failed Count in each sample:
        for (PositionBaseInfo failed : likelyErrors) {
            ++(sampleCounts[failed.readerIndex].failedCount);
        }
        for (final EquivalentIndelRegion failedIndel : list.getFailedIndels()) {
            sampleCounts[failedIndel.sampleIndex].removeIndel(failedIndel);
            ++(sampleCounts[failedIndel.sampleIndex].failedCount);
        }
        list.removeAll(likelyErrors);

    }

    /**
     * Stores counts before any filter was applied. This array is indexed by sci.sampleIndex. Elements
     * of each sample beforeCount list are indexed by genotype index in the sample.
     */
    private IntArrayList[] beforeFilterCounts;

    @Override
    public void preserveCounts(SampleCountInfo[] sampleCounts) {
        int sampleIndex = 0;

        if (beforeFilterCounts == null) {
            beforeFilterCounts=new IntArrayList[sampleCounts.length];
            for (SampleCountInfo sci : sampleCounts) {
                beforeFilterCounts[sci.sampleIndex] = new IntArrayList();
            }
        }
        for (SampleCountInfo sci : sampleCounts) {

            beforeFilterCounts[sci.sampleIndex].clear();

        }
        for (SampleCountInfo sci : sampleCounts) {
            for (int genotypeIndex = 0; genotypeIndex < sci.getGenotypeMaxIndex(); genotypeIndex++) {
                beforeFilterCounts[sci.sampleIndex].add(sci.getGenotypeCount(genotypeIndex));
            }
        }

    }
    // cache counts before filters in the sample counts:

}

