package edu.cornell.med.icb.goby.alignments;

import edu.cornell.med.icb.goby.R.GobyRengine;
import edu.cornell.med.icb.goby.algorithmic.data.CovariateInfo;
import edu.cornell.med.icb.goby.modes.DiscoverSequenceVariantsMode;
import edu.cornell.med.icb.goby.modes.SequenceVariationOutputFormat;
import edu.cornell.med.icb.goby.readers.vcf.ColumnType;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import edu.cornell.med.icb.goby.stats.FisherExactRCalculator;
import edu.cornell.med.icb.goby.stats.VCFWriter;
import edu.cornell.med.icb.goby.util.OutputInfo;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.fastutil.objects.ObjectSortedSet;
import org.apache.commons.math.MathException;
import org.apache.commons.math.MaxIterationsExceededException;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;
import org.apache.commons.math.distribution.PoissonDistributionImpl;
import org.apache.commons.math.stat.inference.ChiSquareTest;
import org.apache.commons.math.stat.inference.ChiSquareTestImpl;
import org.apache.log4j.Logger;
import org.rosuda.JRI.Rengine;

import java.util.Arrays;

/**
 * File format to output genotypes for somatic variations. The format must be used together with the covariates option
 * and a covariate file with the following columns: sample-id, patient-id, kind-of-sample=Germline|Somatic,
 * Parents=P1|P2 (pipe separated list
 * of parents of a somatic sample).
 * <p/>
 * Assuming the following table, genotype frequencies will be compared across the following pairs of samples:
 * S1|S2 vs S3 [determines if S3 has variations not explained by either parent S1 or S2]
 * S3 vs S4    [determines if S3 has variations not found in germline DNA]
 * <p/>
 * These comparisons are determined from the covariates because: P1 and P2 are parents of P3
 * <p/>
 * <table border="1" style="border-collapse:collapse">
 * <tr><td>sample-id</td><td>patient-id</td><td>gender</td><td>type</td><td>kind-of-sample</td><td>tissue</td><td>parents</td><td>offspring</td></tr>
 * <tr><td>S1</td><td>P1</td><td>Male</td><td>Father</td><td>Germline</td><td>Blood</td><td>N/A</td><td>S3|S4</td></tr>
 * <tr><td>S2</td><td>P2</td><td>Female</td><td>Mother</td><td>Germline</td><td>Blood</td><td>N/A</td><td>S3|S4</td></tr>
 * <tr><td>S3</td><td>P3</td><td>Male</td><td>Patient</td><td>Somatic</td><td>Blood</td><td>S1|S2</td><td>N/A</td></tr>
 * <tr><td>S4</td><td>P3</td><td>Male</td><td>Patient</td><td>Germline</td><td>Skin</td><td>S1|S2</td><td>N/A</td></tr>
 * </table>
 *
 * @author Fabien Campagne
 *         Date: 3/9/13
 *         Time: 10:13 AM
 */
public class SomaticVariationOutputFormat implements SequenceVariationOutputFormat {


    protected void setSomaticPValueIndex(int[] somaticPValueIndex) {
        this.somaticPValueIndex = somaticPValueIndex;
    }

    private int somaticPValueIndex[];
    private CovariateInfo covInfo;
    private ObjectArraySet<String> somaticSampleIds;
    GenotypesOutputFormat genotypeFormatter = new GenotypesOutputFormat();
    private static final Logger LOG = Logger.getLogger(SomaticVariationOutputFormat.class);

    /**
     * Given the index of a somatic sample, provides the index of the patient's father sample in sampleCounts.
     * Indices that are not defined have value -1.
     */

    private int[] sample2FatherSampleIndex;
    /**
     * Given the index of a somatic sample, provides the index of the patient's mother sample in sampleCounts.
     * Indices that are not defined have value -1.
     */
    private int[] sample2MotherSampleIndex;
    /**
     * Given the index of a somatic sample, provides the indices of the patient's other germline samples.
     * empty arrays represent samples that have no associated germline samples.
     */
    private int[][] sample2GermlineSampleIndices;
    private int numSamples;
    private boolean fisherRInstalled;
    private int pos;
    private CharSequence currentReferenceId;
    private int referenceIndex;
    private int[] sampleIndex2SomaticSampleIndex;
    private boolean[] isSomatic;
    /**
     * Proportion of total bases observed in a given sample.
     */
    private double[] proportionCountsIn;
    /**
     * Cumulative count of total bases observed in a given sample.
     */
    private double[] countsInSample;

    /**
     * Hook to install mock statsWriter.
     */
    protected void setStatsWriter(VCFWriter statsWriter) {
        this.statsWriter = statsWriter;
    }

    private VCFWriter statsWriter;
    String[] samples;
    private int igvFieldIndex;

    /**
     * Hook to install the somatic sample indices for testing.
     *
     * @param somaticSampleIndices
     */
    protected void setSomaticSampleIndices(IntArrayList somaticSampleIndices) {
        this.somaticSampleIndices = somaticSampleIndices;
    }

    private IntArrayList somaticSampleIndices;

    public void defineColumns(OutputInfo outputInfo, DiscoverSequenceVariantsMode mode) {
        // define columns for genotype format
        samples = mode.getSamples();
        statsWriter = new VCFWriter(outputInfo.getPrintWriter());

        igvFieldIndex = statsWriter.defineField("INFO", "BIOMART_COORDS", 1, ColumnType.String, "Coordinates formatted for use with IGV.");
        genotypeFormatter.defineInfoFields(statsWriter);
        genotypeFormatter.defineGenotypeField(statsWriter);

        covInfo = covInfo != null ? covInfo : mode.getCovariateInfo();
        numSamples = samples.length;
        ObjectSet<String> allCovariates = covInfo.getCovariateKeys();
        if (!allCovariates.contains("patient-id") ||
                !allCovariates.contains("kind-of-sample") ||
                !allCovariates.contains("gender")) {

            System.err.println("SomaticVariationOutputFormat requires the following covariate columns: patient-id, kind-of-sample={Germline|Somatic}, gender={Male|Female}. Please fix the covariate information file. Aborting.");
            System.out.println("The following columns were found in the provided covariate file: " + allCovariates.toString());
            System.exit(1);

        }
        isSomatic = new boolean[numSamples];
        somaticSampleIds = covInfo.samplesWithExactCovariate("kind-of-sample", "Somatic");
        boolean error = false;
        for (String somaticSampleId : somaticSampleIds) {
            int sampleIndex = locateSampleIndex(somaticSampleId);
            if (sampleIndex == -1) {
                System.err.println("Sample id must match between covariate file and alignment basesnames. Mismatch detected for " + somaticSampleId);
                error = true;
            }
        }
        if (error) {
            System.exit(1);
        }
        for (String somaticSampleId : somaticSampleIds) {
            int sampleIndex = locateSampleIndex(somaticSampleId);

            isSomatic[sampleIndex] = true;
        }
        // add column(s) for p-values of somatic variation:
        somaticPValueIndex = new int[numSamples];
        Arrays.fill(somaticPValueIndex, -1);
        setupR();

        for (String sample : somaticSampleIds) {
            int sampleIndex = locateSampleIndex(sample);
            assert sampleIndex != -1 : "sample-id must match between covariate file and alignment basenames.";
            somaticPValueIndex[sampleIndex] = statsWriter.defineField("INFO",
                    String.format("Somatic-P-value(%s)[%s]", fisherRInstalled ? "Fisher" : "Poisson", sample),
                    1, ColumnType.Float,
                    "P-value that a variation is somatic in this particular sample, compared to other germline samples (e.g., germline skin, or mother/father).", "p-value", "statistic", "indexed");
        }

        sample2FatherSampleIndex = new int[numSamples];
        sample2MotherSampleIndex = new int[numSamples];


        Arrays.fill(sample2FatherSampleIndex, -1);
        Arrays.fill(sample2MotherSampleIndex, -1);

        somaticSampleIndices = new IntArrayList();
        for (String somaticSampleId : somaticSampleIds) {
            int sampleIndex = locateSampleIndex(somaticSampleId);
            somaticSampleIndices.add(sampleIndex);
            String parentString = covInfo.getCovariateValue(somaticSampleId, "parents");
            if (parentString != null) {
                // parents column was defined:
                String[] parentIds = parentString.split("[|]");
                for (String parentId : parentIds) {

                    ObjectArraySet<String> parentSamples = covInfo.samplesWithExactCovariate("patient-id", parentId);
                    if (parentSamples.size() != 0) {

                        String parentSampleId = parentSamples.iterator().next();
                        String genderOfParent = covInfo.getCovariateValue(parentSampleId, "gender");
                        int parentSampleIndex = locateSampleIndex(parentSampleId);
                        if (genderOfParent.equals("Male")) {

                            sample2FatherSampleIndex[sampleIndex] = parentSampleIndex;
                        }
                        if (genderOfParent.equals("Female")) {

                            sample2MotherSampleIndex[sampleIndex] = parentSampleIndex;
                        }
                    } else {
                        LOG.warn("Parent could not be found for id:"+parentId);
                    }
                }
            }
        }

        sample2GermlineSampleIndices = new int[numSamples][];
        for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {
            sample2GermlineSampleIndices[sampleIndex] = new int[0];

        }
        for (String somaticSampleId : somaticSampleIds) {
            int sampleIndex = locateSampleIndex(somaticSampleId);
            String patientId = covInfo.getCovariateValue(somaticSampleId, "patient-id");
            ObjectArraySet<String> allSamplesForPatient = covInfo.samplesWithExactCovariate("patient-id", patientId);
            int count = 0;
            for (String sampleId : allSamplesForPatient) {
                if (covInfo.hasCovariateValue(sampleId, "kind-of-sample", "Germline")) {
                    count++;
                }
            }
            int index = 0;
            sample2GermlineSampleIndices[sampleIndex] = new int[count];
            //   System.out.printf("SampleIndex=%d count=%d%n", sampleIndex, count);
            for (String sampleId : allSamplesForPatient) {
                int germlineSampleIndex = locateSampleIndex(sampleId);
                if (covInfo.hasCovariateValue(sampleId, "kind-of-sample", "Germline")) {
                    sample2GermlineSampleIndices[sampleIndex][index++] = germlineSampleIndex;
                }
            }
            assert index == count : "all germline indices must be filled";
        }
        statsWriter.defineSamples(samples);
        statsWriter.setWriteFieldGroupAssociations(true);
        statsWriter.writeHeader();
        countsInSample = new double[numSamples];
        proportionCountsIn = new double[numSamples];
        // initialize sample counts to equal counts across all samples:
        Arrays.fill(countsInSample, 1);

    }

    protected void setupR() {
        final Rengine rEngine = GobyRengine.getInstance().getRengine();
        fisherRInstalled = rEngine != null && rEngine.isAlive();
        //assert fisherRInstalled : "Somatic format requires a working R connection.";
        if (fisherRInstalled) {
            System.err.println("Using FISHER statistics to estimate somatic variation p-values.");
        } else {
            System.err.println("Using Poisson statistics to estimate somatic variation p-values.");
        }
    }

    private int locateSampleIndex(String somaticSampleId) {
        int index = 0;
        for (String sample : samples) {
            if (sample.equals(somaticSampleId)) return index;
            index++;
        }
        return -1;
    }

    public void allocateStorage(int numberOfSamples, int numberOfGroups) {
        this.numSamples = numberOfSamples;
        countsInSample = new double[numberOfSamples];
        proportionCountsIn = new double[numberOfSamples];
        // initialize sample counts to equal counts across all samples:
        Arrays.fill(countsInSample, 1);
        this.numSamples = numberOfSamples;
        genotypeFormatter.allocateStorage(numberOfSamples, numberOfGroups);
    }


    public void writeRecord(final DiscoverVariantIterateSortedAlignments iterator, final SampleCountInfo[] sampleCounts,
                            final int referenceIndex, int position, final DiscoverVariantPositionData list,
                            final int groupIndexA, final int groupIndexB) {
        updateSampleProportions();
        this.pos = position;
        position = position + 1; // report  1-based position
        genotypeFormatter.fillVariantCountArrays(sampleCounts);

        currentReferenceId = iterator.getReferenceId(referenceIndex);

        statsWriter.setId(".");
        statsWriter.setInfo(igvFieldIndex,
                String.format("%s:%d-%d", currentReferenceId, position,
                        position));
        statsWriter.setChromosome(currentReferenceId);

        statsWriter.setPosition(position);

        genotypeFormatter.writeGenotypes(statsWriter, sampleCounts, position);
        if (!statsWriter.hasAlternateAllele()) {
            // do not write a record if the position does not have an alternate allele.
            return;
        }

        // Do not write record if alleleSet is empty, IGV VCF track cannot handle that.

        estimateSomaticPValue(sampleCounts);
        statsWriter.writeRecord();
        updateSampleCumulativeCounts(sampleCounts);
    }

    private void updateSampleCumulativeCounts(SampleCountInfo[] sampleCounts) {
        for (SampleCountInfo info : sampleCounts) {
            // estimate sample proportion with number of reference bases that matched.
            countsInSample[info.sampleIndex] += info.refCount;
        }
    }

    private void updateSampleProportions() {
        long sumAllCounts = 0;

        for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {
            sumAllCounts += countsInSample[sampleIndex];
        }
        for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {

            proportionCountsIn[sampleIndex] = ((double) countsInSample[sampleIndex]) / (double) sumAllCounts;
        }
    }

    /**
     * Calculate the proportion of bases that we would expect across a pair of sample.
     *
     * @param sampleIndexA first sample index of the pair under consideration
     * @param sampleIndexB second sample index of the pair under consideration
     * @param request      index of the sample for which the proportion should be returned.
     * @return proportion of bases for request sample.
     */
    private double getSpecificSampleProportion(int sampleIndexA, int sampleIndexB, int request) {
        long sumAllCounts = 0;

        for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {
            if (sampleIndex == sampleIndexA || sampleIndex == sampleIndexB) {
                sumAllCounts += countsInSample[sampleIndex];
            }
        }
        if (request == sampleIndexA) return ((double) countsInSample[sampleIndexA]) / (double) sumAllCounts;
        else if (request == sampleIndexB) return ((double) countsInSample[sampleIndexB]) / (double) sumAllCounts;
        else throw new IllegalArgumentException("request must be one of sampleIndexA or sampleIndexB");
    }


    public void close() {

        statsWriter.close();
    }

    public void setGenome(RandomAccessSequenceInterface genome) {
        genotypeFormatter.setGenome(genome);
    }

    public void setGenomeReferenceIndex(int index) {
        genotypeFormatter.setGenomeReferenceIndex(index);
    }

    protected DoubleArrayList pValues = new DoubleArrayList();

    protected void setSample2FatherSampleIndex(int[] sample2FatherSampleIndex) {
        this.sample2FatherSampleIndex = sample2FatherSampleIndex;
    }

    protected void setSample2GermlineSampleIndices(int[][] sample2GermlineSampleIndices) {
        this.sample2GermlineSampleIndices = sample2GermlineSampleIndices;
    }

    protected void setSample2MotherSampleIndex(int[] sample2MotherSampleIndex) {
        this.sample2MotherSampleIndex = sample2MotherSampleIndex;
    }

    protected void estimateSomaticPValue(SampleCountInfo[] sampleCounts) {
       // compare somatic to father:

        for (int sampleIndex : somaticSampleIndices) {
            pValues.clear();

            SampleCountInfo somaticCounts = sampleCounts[sampleIndex];
            int fatherSampleIndex = sample2FatherSampleIndex[sampleIndex];
            if (fatherSampleIndex != -1) {

                SampleCountInfo fatherCounts = sampleCounts[fatherSampleIndex];
                double fatherP = estimateP(somaticCounts, fatherCounts);
                pValues.add(fatherP);
            }
            int motherSampleIndex = sample2MotherSampleIndex[sampleIndex];
            if (fatherSampleIndex != -1) {

                SampleCountInfo motherCounts = sampleCounts[motherSampleIndex];
                double motherP = estimateP(somaticCounts, motherCounts);
                pValues.add(motherP);
            }
            int germlineSampleIndices[] = sample2GermlineSampleIndices[sampleIndex];
            for (int germlineSampleIndex : germlineSampleIndices) {

                SampleCountInfo germlineCounts = sampleCounts[germlineSampleIndex];
                double germlineP = estimateP(somaticCounts, germlineCounts);
                pValues.add(germlineP);
            }
            // use the max of the above p-values:
            double pValue = max(pValues);
            statsWriter.setInfo(somaticPValueIndex[sampleIndex], pValue);
        }
    }

    private double max(DoubleArrayList pValues) {
        if (pValues.size() == 0) return Double.NaN;
        double max = Double.MIN_VALUE;
        for (double v : pValues) {
            max = Math.max(v, max);
        }
        return max;
    }

    protected double estimateP(SampleCountInfo somaticCounts, SampleCountInfo germlineCounts) {

        double fisherP = 1;


        for (int genotypeIndex = 0; genotypeIndex < somaticCounts.getGenotypeMaxIndex(); genotypeIndex++) {


            boolean ok = checkCounts(somaticCounts, germlineCounts, genotypeIndex);
            final int a = germlineCounts.getGenotypeCount(genotypeIndex);
            final int b = somaticCounts.getGenotypeCount(genotypeIndex);
            if (a >= 1) {
                // at least 1 base in germline, we don't need to score this:
                fisherP = Math.min(fisherP, 1);
                continue;
            }
            //    final int c = (int) ((a + b) * proportionCountsIn[germlineCounts.sampleIndex]);
            //    final int d = (int) ((a + b) * proportionCountsIn[somaticCounts.sampleIndex]);
            double germlineSampleExpectedProportion = getSpecificSampleProportion(germlineCounts.sampleIndex,
                    somaticCounts.sampleIndex, germlineCounts.sampleIndex);
            double somaticSampleExpectedProportion = getSpecificSampleProportion(germlineCounts.sampleIndex,
                    somaticCounts.sampleIndex, somaticCounts.sampleIndex);
            assert germlineSampleExpectedProportion + somaticSampleExpectedProportion == 1 : "proportions must sum to 1.0";
            final int c = (int) ((a + b) * germlineSampleExpectedProportion);
            final int d = (int) ((a + b) * somaticSampleExpectedProportion);
            if (ok) {
                if (fisherRInstalled) {
                    //    System.out.printf("%n fisher=? %n%d|%d%n%d|%d ", a, b - a, c, d - c);
                    fisherP = Math.min(fisherP, FisherExactRCalculator.getFisherOneTailedLesserPValue(
                            a, b,
                            c, d));
                    //System.out.printf("%n fisher=%g %n%d|%d%n%d|%d ", fisherP, a, b, c, d);
                    //System.out.flush();
                } else {

                    int mean = Math.max(1, (a + c) / 2);

                    final PoissonDistributionImpl poissonDistribution = new PoissonDistributionImpl(mean);
                    try {
                        fisherP = Math.min(fisherP, poissonDistribution.cumulativeProbability(c));
                        // System.out.printf("%nmean=%d%n%d|%d%n%d|%d  x=%d  P(X<=x) %g",mean, a,b,c,d,c, fisherP);
                    } catch (MathException e) {
                        fisherP = Double.NaN;
                    }

                }

            } else {
                System.err.printf("An exception was caught evaluating the Fisher Exact test P-value. Details are provided below%n" +
                        "referenceId=%s referenceIndex=%d position=%d %n" +
                        "a=%d b=%d%n" +
                        "c=%d, d=%d",
                        currentReferenceId, referenceIndex,
                        pos + 1,
                        a, b, c, d
                );
            }

        }
        return fisherP != -1 ? fisherP : Double.NaN;
    }

    private boolean checkCounts(SampleCountInfo aCounts, SampleCountInfo bCounts, int genotypeIndex) {

        boolean ok = true;
        // detect if any count is negative (that's a bug)
        int count = aCounts.getGenotypeCount(genotypeIndex);

        if (count < 0) ok = false;
        count = bCounts.getGenotypeCount(genotypeIndex);

        if (count < 0) ok = false;


        return ok;

    }

    public void setCovariateInfo(CovariateInfo covInfo) {
        this.covInfo = covInfo;
    }
}
