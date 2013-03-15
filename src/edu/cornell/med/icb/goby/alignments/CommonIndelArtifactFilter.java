package edu.cornell.med.icb.goby.alignments;

import edu.cornell.med.icb.goby.algorithmic.data.EquivalentIndelRegion;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.PoissonDistributionImpl;

/**
 * Filter to remove common sequencing artifacts that look like indels. Such artifacts occur in a region of the genome
 * where a base is repeated in a stretch of residue. Many sequencing platforms loose phase in these regions, resulting
 * in apparent sequence with one or more bases added or removed from the true sequence.
 * Filtering these artifacts consists in counting how many repeat bases have been seen, and the number of indels.
 * This counts provide the expected number of indel artifacts. An indel is called only when the observed frequency is
 * larger than expected (P<0.05, with Poisson cumulative distribution) from prior observations (most of them are assumed to be errors).
 *
 * @author Fabien Campagne
 *         Date: 3/15/13
 *         Time: 3:20 PM
 */
public class CommonIndelArtifactFilter extends GenotypeFilter {
    /**
     * The number of bases in indel repeats seen so far.
     */
    private int cumulativeIndelBaseRepeatLength = 1;
    private int cumulativeIndelFrequency = 1;
    private double proportionFreqOverLength = 1;

    @Override
    public void filterGenotypes(DiscoverVariantPositionData list, SampleCountInfo[] sampleCounts, ObjectSet<PositionBaseInfo> filteredSet) {
        if (list.hasCandidateIndels()) {
            for (EquivalentIndelRegion indel : list.getIndels()) {
                if (indel != null) {
                    int lengthRepeatBases = countRepetitiveBases(indel) * indel.getFrequency();
                    int expectedFrequency = (int) (lengthRepeatBases * proportionFreqOverLength);
                    if (expectedFrequency > 0) {
                        if (indel.from.contains("AAAAAAAAAAAA") || indel.to.equals("AAAAAAAAAAAA")) {
                            System.out.println("STOP");
                        }
                        try {

                            cumulativeIndelBaseRepeatLength += countRepetitiveBases(indel) * indel.getFrequency();
                            cumulativeIndelFrequency += indel.getFrequency();

                            int actualObservedBaseIndels = lengthRepeatBases;
                            PoissonDistributionImpl poisson = new PoissonDistributionImpl(actualObservedBaseIndels);

                            if (actualObservedBaseIndels <= expectedFrequency || poisson.cumulativeProbability(expectedFrequency) > 0.05) {

                                // fail this index, since it is expected given the number of repeat bases.
                                System.out.printf("Failing INDEL=%s Expected prop=%g freq=%d observed freq=%d%n",
                                        indel,
                                        proportionFreqOverLength, expectedFrequency, indel.getFrequency());

                                indel.markFiltered();
                                list.failIndel(indel);
                                sampleCounts[indel.sampleIndex].removeIndel(indel);
                            } else {
                                System.out.println("Not filtering " + indel);
                            }
                        } catch (MathException e) {
                            System.err.println(e);
                            LOG.error(e);
                        }
                    }
                }

            }

            // update cumulative counts and proportions:
            for (EquivalentIndelRegion indel : list.getIndels()) {
                if (indel != null) {

                }
            }

            estimateProportions();
        }


    }

    private void estimateProportions() {
        //  System.out.printf("cumulative %d %d %n", cumulativeIndelFrequency, cumulativeIndelBaseRepeatLength);

        proportionFreqOverLength = ((double) cumulativeIndelFrequency) / (double) cumulativeIndelBaseRepeatLength;
    }

    protected int countRepetitiveBases(EquivalentIndelRegion indel) {
        if (indel == null || indel.from == null) return 0;
        if (indel.from.contains("-")) {
            return countRepetitiveBases(indel.to);
        } else {
            return countRepetitiveBases(indel.from);
        }
    }

    protected int countRepetitiveBases(String bases) {
        if (bases.length() == 0) return 0;
        int count = 0;
        char previous = bases.charAt(0);
        char currentBase = '\0';
        int length = bases.length();
        for (int i = 1; i < length; i++) {
            currentBase = bases.charAt(i);
            if (currentBase == previous && currentBase != '-') {
                count++;

            }
            previous = currentBase;
        }
        return count;
    }

    @Override
    public int getThresholdForSample(int sampleIndex) {
        return 0;
    }
}
