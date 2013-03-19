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
                    int depthAtPosition = 0;
                    for (SampleCountInfo sci : sampleCounts) {
                        depthAtPosition += sci.getSumCounts();
                    }
                    final int repeatLength = countRepetitiveBases(indel);
                    final int actualIndelFrequency = indel.getFrequency();
                    int gapLength = calculateGapLength(indel);
                    double indelFractionOfDepth = actualIndelFrequency / depthAtPosition;
                    int depthLessThan10 = depthAtPosition < 10 ? 1 : 0;
                    int expectedFrequency = getPredictedIndelFrequency3(repeatLength, gapLength, depthAtPosition, repeatPatternLength(indel), depthLessThan10);
                    if (expectedFrequency > 0) {

                        try {
                            cumulativeIndelBaseRepeatLength += repeatLength * actualIndelFrequency;
                            cumulativeIndelFrequency += actualIndelFrequency;

                            PoissonDistributionImpl poisson = new PoissonDistributionImpl(actualIndelFrequency);
                            double pValue = poisson.cumulativeProbability(expectedFrequency);
                         /*   System.out.printf("@INDEL_DATA@%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%g\t%s%n",
                                    indel.from, indel.to,  depthAtPosition, repeatLength, gapLength,
                                    repeatPatternLength(indel), depthLessThan10, expectedFrequency, actualIndelFrequency, pValue,
                                    pValue > 0.05 ? "FILTER" : "KEPT");
                           */
                            if (actualIndelFrequency <= expectedFrequency || pValue > 0.05) {

                                // fail this index, since it is expected given the number of repeat bases.
                             /*   System.out.printf("Failing INDEL=%s Expected prop=%g freq=%d observed freq=%d%n",
                                        indel,
                                        proportionFreqOverLength, expectedFrequency, actualIndelFrequency);
                               */
                                indel.markFiltered();

                                //sampleCounts[indel.sampleIndex].removeIndel(indel);
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

    protected int repeatPatternLength(EquivalentIndelRegion indel) {
        return Math.max(repeatPatternLength(indel.from, indel.to),
                repeatPatternLength(indel.to, indel.from));
    }

    private int calculateGapLength(EquivalentIndelRegion indel) {
        return Math.max(gapLength(indel.from), gapLength(indel.to));
    }

    protected int gapLength(String to) {
        int length = 0;
        for (int i = 0; i < to.length(); i++) {
            if (to.charAt(i) == '-') {
                length++;
            }
        }
        return length;
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

    // use a trained neural net to determine the expected indel frequency, from features about the indel:
    public int getPredictedIndelFrequency3(int repeatLength, int gapLength, int depthAtPosition, int repeatPatternLengthIndel, int depthLessThan10) {
      /*%PRODUCER: JMP - Neural */
      /*%TARGET: Log10_actual_indel_freq */
      /*%INPUT: repeatLength */
      /*%INPUT: gapLength */
      /*%INPUT: repeatPatternLengthIndel */
      /*%INPUT: log10_Depth_at_position */
      /*%OUTPUT: Log10_actual_indel_freq_Predicted */
        double log10_Depth_at_position = Math.log10(depthAtPosition);
      /* Transformation Code */

      /* Hidden Layer Code */
        double H1 = Math.tanh(.5 * (-0.217511440307164 * repeatLength + -0.153912432260887 * gapLength + -0.626691444965826 * repeatPatternLengthIndel + 1.07492523248614 * log10_Depth_at_position + -1.69089686523462));
        double H2 = 0.068941619222792 * repeatLength + -0.0863283200708544 * gapLength + 0.00229027639005201 * repeatPatternLengthIndel + -3.84733744230361 * log10_Depth_at_position + 11.7784564182034;
        double HH1 = Math.tanh(.5 * (7.08109026103453 * H1 + 1.27426071706873 * H2 + 1.10315736842476));
        double HH2 = Math.tanh(.5 * (-0.0150959912141104 * H1 + -1.46380340536445 * H2 + -1.2308961426117));
        double HH3 = Math.tanh(.5 * (1.53535193119845 * H1 + -0.562710832434124 * H2 + 1.51259843151645));
        double HH4 = -0.480384063647369 * H1 + -0.0935191076999823 * H2 + -0.373846155675067;
        double HH5 = -0.323686814060226 * H1 + -0.0216726416714448 * H2 + -0.154578788760765;
        double HH6 = -0.0364869192793619 * H1 + 0.0235850981885036 * H2 + -0.768639170264522;

      /* Final Layer Code */
        double THETA1 = -0.200181906310874 * HH1 + 0.194417036932698 * HH2 + -0.394347363656623 * HH3 + 0.897881984682189 * HH4 + -1.14470827545174 * HH5 + -1.68508250216669 * HH6 + 0.360199064456832;

      /* Response Mapping Code */
        double Log10_actual_indel_freq_Predicted = THETA1;
        int actualIndelFrequency_Predicted = (int) Math.max(0, Math.pow(10, Log10_actual_indel_freq_Predicted));
        return actualIndelFrequency_Predicted;
    }

    public int repeatPatternLength(String seq1, String seq2) {
        int gapLength = gapLength(seq1);
        String seq1NoGap = seq1.replaceAll("-", "");
        String seq2NoGap = seq2.replaceAll("-", "");
        int i;
        int maxIndex = Math.min(Math.min(gapLength, seq1NoGap.length()), seq1NoGap.length());

        for (i = 0; i < maxIndex; i++) {

            if (seq1NoGap.charAt(i) != seq2NoGap.charAt(i)) {
                break;
            }

        }
        return i;
    }
}
