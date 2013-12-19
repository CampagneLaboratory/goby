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
    private static final int POSITIVE_STRAND = 1;
    private static final int NEGATIVE_STRAND = 0;

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
            boolean genotypeFilteredInAllSamples = true;
            for (SampleCountInfo sci : sampleCounts) {

                genotypePresentInSomeSample |= getBeforeFilterCount(genotypeIndex, sci) > 0;
                genotypeRemainsInSomeSample |= sci.getGenotypeCount(genotypeIndex) > 0;
                genotypeFilteredInAllSamples &= sci.isFiltered(genotypeIndex);
            }
            // when the genotype was present, and it remains after filtering in some sample, we revert the effect of the
            // filters in all samples:
            if (genotypePresentInSomeSample && genotypeRemainsInSomeSample) {
                // we need to revert filtering for this genotype to prevent artifactual thresholding effects.
                for (SampleCountInfo sci : sampleCounts) {
                    resetCount(genotypeIndex, sci);

                    if (sci.isIndel(genotypeIndex)) {
                        //     System.out.println("reverting indel genotype.");
                        EquivalentIndelRegion indel = sci.getIndelGenotype(genotypeIndex);
                        indel.removeFiltered();
                        list.getFailedIndels().remove(indel);
                    }
                }
            } else {
                if (genotypeFilteredInAllSamples) {

                    for (SampleCountInfo sci : sampleCounts) {
                        if (sci.isIndel(genotypeIndex)) {
                            EquivalentIndelRegion indel = sci.getIndelGenotype(genotypeIndex);
                            list.failIndel(indel);
                            sci.removeIndel(indel);
                            sci.setGenotypeCount(genotypeIndex, 0);
                        }
                    }
                }
            }
        }
        // do not decrement counts again. The Filters have done this already..

        for (SampleCountInfo sci : sampleCounts) {
            for (int i = 0; i < sci.getCountsSize(); i++) {
                assert sci.getGenotypeCount(i) >= 0 : ("Counts must never be negative. This would happen if a GenotypeFilter removed counts directly. Value was " + sci.getGenotypeCount(i));
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

    private void resetCount(int genotypeIndex, SampleCountInfo sci) {

        sci.setGenotypeCount(genotypeIndex, beforeFilterCounts[POSITIVE_STRAND][sci.sampleIndex].getInt(genotypeIndex), true);
        sci.setGenotypeCount(genotypeIndex, beforeFilterCounts[NEGATIVE_STRAND][sci.sampleIndex].getInt(genotypeIndex), false);
    }

    private int getBeforeFilterCount(int genotypeIndex, SampleCountInfo sci) {
        return beforeFilterCounts[0][sci.sampleIndex].getInt(genotypeIndex)+beforeFilterCounts[1][sci.sampleIndex].getInt(genotypeIndex);
    }

    /**
     * Stores counts before any filter was applied. This array is indexed by strand (0 or 1, where 1 is positive strand)
     * and sci.sampleIndex . Elements of each sample beforeCount list are indexed by genotype index in the sample.
     */
    private IntArrayList[][] beforeFilterCounts;

    @Override
    public void preserveCounts(SampleCountInfo[] sampleCounts) {
        int sampleIndex = 0;

        if (beforeFilterCounts == null) {
            beforeFilterCounts = new IntArrayList[2][];
            beforeFilterCounts[NEGATIVE_STRAND] = new IntArrayList[sampleCounts.length];
            beforeFilterCounts[POSITIVE_STRAND] = new IntArrayList[sampleCounts.length];
            for (SampleCountInfo sci : sampleCounts) {
                beforeFilterCounts[NEGATIVE_STRAND][sci.sampleIndex] = new IntArrayList();
                beforeFilterCounts[POSITIVE_STRAND][sci.sampleIndex] = new IntArrayList();
            }
        }
        for (SampleCountInfo sci : sampleCounts) {

            beforeFilterCounts[NEGATIVE_STRAND][sci.sampleIndex].clear();
            beforeFilterCounts[POSITIVE_STRAND][sci.sampleIndex].clear();

        }
        for (SampleCountInfo sci : sampleCounts) {
            for (int genotypeIndex = 0; genotypeIndex < sci.getGenotypeMaxIndex(); genotypeIndex++) {
                beforeFilterCounts[POSITIVE_STRAND][sci.sampleIndex].add(sci.getGenotypeCount(genotypeIndex,true));
                beforeFilterCounts[NEGATIVE_STRAND][sci.sampleIndex].add(sci.getGenotypeCount(genotypeIndex,false));
            }
        }

    }
    // cache counts before filters in the sample counts:

}

