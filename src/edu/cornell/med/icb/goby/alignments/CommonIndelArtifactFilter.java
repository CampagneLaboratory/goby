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
                    int lengthRepeatBases = repeatLength * actualIndelFrequency;
                    int gapLength = calculateGapLength(indel);
                    double indelFractionOfDepth = actualIndelFrequency / depthAtPosition;
                    //getPredictedIndelFrequency3(int repeatLength, int gapLength, double indelFractionOfDepth, int repeatPatternLength, int depthLessThan10)
                    int expectedFrequency = getPredictedIndelFrequency3(repeatLength, gapLength, indelFractionOfDepth, repeatPatternLength(indel),depthAtPosition<10?1:0);
                    if (expectedFrequency > 0) {

                        try {
                            cumulativeIndelBaseRepeatLength += repeatLength * actualIndelFrequency;
                            cumulativeIndelFrequency += actualIndelFrequency;


                            PoissonDistributionImpl poisson = new PoissonDistributionImpl(actualIndelFrequency);
                            double pValue = poisson.cumulativeProbability(expectedFrequency);
                            System.out.printf("@INDEL_DATA@%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%g\t%s%n",
                                    indel.from, indel.to, actualIndelFrequency, depthAtPosition, repeatLength, gapLength,
                                    lengthRepeatBases, repeatPatternLength(indel),
                                    expectedFrequency, actualIndelFrequency, pValue, pValue > 0.05 ? "FILTER" : "KEPT");

                            if (actualIndelFrequency <= expectedFrequency || poisson.cumulativeProbability(expectedFrequency) > 0.05) {

                                // fail this index, since it is expected given the number of repeat bases.
                             /*   System.out.printf("Failing INDEL=%s Expected prop=%g freq=%d observed freq=%d%n",
                                        indel,
                                        proportionFreqOverLength, expectedFrequency, actualIndelFrequency);
                               */
                                indel.markFiltered();
                                list.failIndel(indel);
                                sampleCounts[indel.sampleIndex].removeIndel(indel);
                            } else {
                                //  System.out.println("Not filtering " + indel);
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

    // use a trained neural net to determine the expected indel frequency, from the fraction of indel at depth, repeat length and gapLength
    public int getPredictedIndelFrequency(int repeatLength, int gapLength, double indelFractionOfDepth, int depthLessThan10) {
        /*%PRODUCER: JMP - Neural */
        /*%TARGET: indel_getFrequency_ */
        /*%INPUT: repeatLength */
        /*%INPUT: gapLength */
        /*%INPUT: indelFractionOfDepth */
        /*%INPUT: depthLessThan10 */
        /*%OUTPUT: indel_getFrequency__Predicted */

        /* Transformation Code */

        /* Hidden Layer Code */
        double H1 = Math.tanh(.5 * (0.27356506270027 * repeatLength + 0.166529009810623 * gapLength + 27.3279278140353 * indelFractionOfDepth + -3.2012277396908 * depthLessThan10 + -2.36610247334202));
        double H2 = Math.tanh(.5 * (-0.198507001359369 * repeatLength + -0.0762741349243423 * gapLength + -4.30342582957599 * indelFractionOfDepth + -8.93908165427444 * depthLessThan10 + 2.32056339497897));
        double H3 = Math.tanh(.5 * (-0.193335839754507 * repeatLength + -0.0705430261002909 * gapLength + -23.8750875209862 * indelFractionOfDepth + -13.3575848303556 * depthLessThan10 + 4.85244737726472));

        /* Final Layer Code */
        double THETA1 = 18.3130689781713 * H1 + 46.7611929833778 * H2 + -26.3989806176432 * H3 + 8.92550399206593;

        /* Response Mapping Code */
        int indel_getFrequency__Predicted = (int) THETA1;
        return Math.max(0, indel_getFrequency__Predicted);
    }

    // use a trained neural net to determine the expected indel frequency, from the fraction of indel at depth, repeat length and gapLength
    public int getPredictedIndelFrequency2(int repeatLength, int gapLength, double depthAtPosition, int lengthRepeatBases) {
       /*%PRODUCER: JMP - Neural */
       /*%TARGET: actualIndelFrequency */
       /*%INPUT: depthAtPosition */
       /*%INPUT: repeatLength */
       /*%INPUT: gapLength */
       /*%INPUT: lengthRepeatBases */
       /*%OUTPUT: actualIndelFrequency_Predicted */
       /* Transformation Code */

       /* Hidden Layer Code */
        double H1 = Math.tanh(.5 * (-0.000086554144926645 * depthAtPosition + 3.16907492832398 * repeatLength + -0.0250392599204096 * gapLength + 0.00860053382931063 * lengthRepeatBases + -1.8028891458221));
        double H2 = Math.tanh(.5 * (0.0000599121636027239 * depthAtPosition + -1.79882893444531 * repeatLength + 0.0221281375539389 * gapLength + 0.00637589716134292 * lengthRepeatBases + -0.348968572114815));
        double H3 = Math.tanh(.5 * (-0.000001671794178104 * depthAtPosition + -0.0015960570685604 * repeatLength + -0.00126693572993562 * gapLength + 0.00465859744456237 * lengthRepeatBases + 0.274169573724267));

       /* Final Layer Code */
        double THETA1 = 109.812429539664 * H1 + 233.811671806761 * H2 + 144.232437536586 * H3 + 100.624807326541;

       /* Response Mapping Code */
        int actualIndelFrequency_Predicted = (int) Math.max(0, THETA1);
        return actualIndelFrequency_Predicted;
    }  // use a trained neural net to determine the expected indel frequency, from the fraction of indel at depth, repeat length and gapLength
    public int getPredictedIndelFrequency3(int repeatLength, int gapLength, double indelFractionOfDepth, int repeatPatternLength, int depthLessThan10) {
      /*%PRODUCER: JMP - Neural */
      /*%TARGET: actualIndelFrequency */
      /*%INPUT: depthLessThan10 */
      /*%INPUT: indelFractionOfDepth */
      /*%INPUT: repeatLength */
      /*%INPUT: gapLength */
      /*%INPUT: repeatPatternLength */
      /*%OUTPUT: actualIndelFrequency_Predicted */

      /* Transformation Code */

      /* Hidden Layer Code */
    double  H1 = Math.tanh(.5 * (30.7569539696584 * depthLessThan10 + -28.6613224227517 * indelFractionOfDepth + -0.141307074837274 * repeatLength + 0.0127209881460977 * gapLength + 0.0815326097409295 * repeatPatternLength + 4.39285626019316));
        double  H2 = 147.329633240459*depthLessThan10 + -80.1455280830311*indelFractionOfDepth + -4.14627675030666*repeatLength + 0.437086814147535*gapLength + 2.69260298511829*repeatPatternLength + 22.5654455372969;
        double  HH1 = Math.tanh(.5 * (2.46731575079385 * H1 + -0.0861261422718692 * H2 + -0.14432842314448));
        double  HH2 = Math.tanh(.5 * (-0.976232159013453 * H1 + 1.51964305777285 * H2 + -18.4356764270507));
        double  HH3 = Math.tanh(.5 * (-74.5927953058719 * H1 + -0.0722113892561223 * H2 + 75.1015067951385));
        double  HH4 = -0.0955052634966312*H1 + -0.0935448245948543*H2 + -3.14711649420219;;
        double  HH5 = 0.168513127901274*H1 + 0.28526253608911*H2 + 6.47779522886649;;

      /* Final Layer Code */
        double  THETA1=-25.5467286027933*HH1 + -4.66183256485626*HH2 + 19.1258767129895*HH3 + 60.8855184647536*HH4 + 19.8522989361628*HH5 + 69.8630756526183;

      /* Response Mapping Code */

        int actualIndelFrequency_Predicted = (int) Math.max(0, THETA1);
        return actualIndelFrequency_Predicted;
    }

    public int repeatPatternLength(String seq1, String seq2) {
        int gapLength = gapLength(seq1);
        String seq1NoGap = seq1.replaceAll("-", "");
        String seq2NoGap = seq2.replaceAll("-", "");
        int i;
        int maxIndex=Math.min(Math.min(gapLength, seq1NoGap.length()),seq1NoGap.length());

        for (i = 0; i < maxIndex; i++) {

            if (seq1NoGap.charAt(i) != seq2NoGap.charAt(i)) {
                break;
            }

        }
        return i;
    }
}
