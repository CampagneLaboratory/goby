package edu.cornell.med.icb.goby.alignments;

import edu.cornell.med.icb.goby.util.dynoptions.DynamicOptionClient;
import edu.cornell.med.icb.goby.util.dynoptions.RegisterThis;
import it.unimi.dsi.fastutil.objects.ObjectSet;

/**
 * Removes genotypes when only one strand supports the genotype in a sample.
 *
 * @author Fabien Campagne
 *         Date: 11/1/13
 *         Time: 2:24 PM
 */
public class StrandBiasFilter extends GenotypeFilter {
    private final int jPlusOne;
    private int jParameter = 9;

    @RegisterThis
    public static final DynamicOptionClient doc = new DynamicOptionClient(QualityScoreFilter.class, "j:The maximum " +
            "number of bases below after which we require base support from both strands. When j=9, variations with 10 bases" +
            "must have representation from both strands. Please note that values of j lower than 1 are not valid.:9");

    public StrandBiasFilter(int jParameter) {
        this.jParameter = jParameter;
        jPlusOne = jParameter + 1;
    }

    public StrandBiasFilter() {
        this.jParameter = doc.getInteger("j");
        assert jParameter >= 1;
        jPlusOne = jParameter + 1;
    }

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
            final int countForward = sampleCountInfo.getGenotypeCount(baseIndex, true);
            final int countReverse = sampleCountInfo.getGenotypeCount(baseIndex, false);

            if ((countForward + countReverse) >= (jPlusOne) && (countForward == 0 || countReverse == 0)) {

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
