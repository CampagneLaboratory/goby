package edu.cornell.med.icb.goby.alignments;

import edu.cornell.med.icb.goby.algorithmic.data.EquivalentIndelRegion;
import edu.cornell.med.icb.goby.util.dynoptions.DynamicOptionClient;
import edu.cornell.med.icb.goby.util.dynoptions.RegisterThis;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.PoissonDistributionImpl;

import java.util.Arrays;

/**
 * Filter to remove indels at a site where a sample shows lots of distinct possible indels. Indels at these
 * sites are very likely to be artefactual. We count the number of samples where three distinct indel genotypes
 * are seen. If more than 1/4 of the samples have likely indel artifacts, we remove all indel candidates at the site.
 *
 * @author Fabien Campagne
 *         Date: 3/15/13
 *         Time: 3:20 PM
 */
public class EntropicIndelArtifactFilter extends GenotypeFilter {
    @RegisterThis
    public static final DynamicOptionClient doc = new DynamicOptionClient(EntropicIndelArtifactFilter.class,
            "maxIndelPerSite:Maximum number of distinct indels at a given genomic site.:1",
            "fractionOfSamples:Maximum fraction of samples that can have an indel candidate for the indel to be considered (indel candidates that occur in many samples are more likely to be spurious).:0.25");

    private int maxIndelPerSite=1;
    private double fractionOfSamples;

    public static DynamicOptionClient doc() {
        return doc;
    }

    public EntropicIndelArtifactFilter() {
        maxIndelPerSite = doc.getInteger("maxIndelPerSite");
        fractionOfSamples = doc.getDouble("fractionOfSamples");
    }
    private int[] distinctIndelsWithCount;

    public void initStorage(int numSamples) {
        super.initStorage(numSamples);
        if (distinctIndelsWithCount == null) {
            distinctIndelsWithCount = new int[numSamples];
            candidateIndels = new ObjectArraySet<EquivalentIndelRegion>();

        } else {
            Arrays.fill(distinctIndelsWithCount, 0);
            candidateIndels.clear();
        }
    }

    ObjectArraySet<EquivalentIndelRegion> candidateIndels;

    @Override
    public void filterGenotypes(DiscoverVariantPositionData list,
                                SampleCountInfo[] sampleCounts,
                                ObjectSet<PositionBaseInfo> filteredSet) {
        resetCounters();
        initStorage(sampleCounts.length);
        int likelyIndelArtifact = 0;
        for (SampleCountInfo sci : sampleCounts) {
            if (sci.hasIndels()) {
                for (EquivalentIndelRegion indel : sci.getEquivalentIndelRegions()) {
                    if (indel.getFrequency() > 0) {
                        distinctIndelsWithCount[sci.sampleIndex]++;
                        candidateIndels.add(indel);
                        numScreened++;
                    }
                }
                if (distinctIndelsWithCount[sci.sampleIndex] > 1) {
                    likelyIndelArtifact++;
                }
            }
        }


        if (candidateIndels.size() > maxIndelPerSite || likelyIndelArtifact >= (sampleCounts.length *fractionOfSamples)) {
            // more than a quarter of samples seem to have indel artifacts at this genomic site. Filter everything.

            //   System.out.printf("Filtering indels at site with %d samples with likely indel artifact %n",                    likelyIndelArtifact);
            for (SampleCountInfo sci : sampleCounts) {
                final int disctingIndelsInSample = distinctIndelsWithCount[sci.sampleIndex];
                for (EquivalentIndelRegion indel : sci.getEquivalentIndelRegions()) {

                    numFiltered++;
                    indel.markFiltered();
                }
            }
        }
    }


    @Override
    public int getThresholdForSample(int sampleIndex) {
        return 0;
    }

}
