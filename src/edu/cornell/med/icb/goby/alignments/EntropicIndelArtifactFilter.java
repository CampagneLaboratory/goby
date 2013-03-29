package edu.cornell.med.icb.goby.alignments;

import edu.cornell.med.icb.goby.algorithmic.data.EquivalentIndelRegion;
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
    private int[] distinctIndelsWithCount;

    public void initStorage(int numSamples) {
        super.initStorage(numSamples);
        if (distinctIndelsWithCount == null) {
            distinctIndelsWithCount = new int[numSamples];
        } else {
            Arrays.fill(distinctIndelsWithCount, 0);
        }
    }

    @Override
    public void filterGenotypes(DiscoverVariantPositionData list,
                                SampleCountInfo[] sampleCounts,
                                ObjectSet<PositionBaseInfo> filteredSet) {
        resetCounters();
        initStorage(sampleCounts.length);

        int likelyIndelArtifact = 0;
        for (SampleCountInfo sci : sampleCounts) {

            for (EquivalentIndelRegion indel : sci.getEquivalentIndelRegions()) {
                if (indel.getFrequency() > 0) {
                    distinctIndelsWithCount[sci.sampleIndex]++;
                }
            }
            if (distinctIndelsWithCount[sci.sampleIndex] >= 2) {
                likelyIndelArtifact++;
            }

        }


        if (likelyIndelArtifact >= sampleCounts.length / 4) {
            // more than a quarter of samples seem to have indel artifacts at this genomic site. Filter everything.

            //   System.out.printf("Filtering indels at site with %d samples with likely indel artifact %n",                    likelyIndelArtifact);
            for (SampleCountInfo sci : sampleCounts) {
                final int disctingIndelsInSample = distinctIndelsWithCount[sci.sampleIndex];
                for (EquivalentIndelRegion indel : sci.getEquivalentIndelRegions()) {


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
