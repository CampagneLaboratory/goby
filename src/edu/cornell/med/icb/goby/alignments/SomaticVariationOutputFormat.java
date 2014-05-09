package edu.cornell.med.icb.goby.alignments;

import edu.cornell.med.icb.goby.R.GobyRengine;
import edu.cornell.med.icb.goby.algorithmic.data.CovariateInfo;
import edu.cornell.med.icb.goby.modes.DiscoverSequenceVariantsMode;
import edu.cornell.med.icb.goby.modes.SequenceVariationOutputFormat;
import edu.cornell.med.icb.goby.readers.vcf.ColumnType;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import edu.cornell.med.icb.goby.stats.VCFWriter;
import edu.cornell.med.icb.goby.util.OutputInfo;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.PoissonDistributionImpl;
import org.apache.log4j.Logger;
import org.rosuda.JRI.Rengine;

import java.io.IOException;
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

    /**
     * We will store the largest candidate somatic frequency here.
     */
    private int[] candidateFrequencyIndex;

    /**
     * A priority score for the somatic site. Larger integer values indicate more support for the site being a
     * somatic variation. Indexed by sampleIndex.
     */
    private int[] maxGenotypeSomaticPriority;

    private int strictSomaticCandidateFieldIndex;
    /**
     * If a somatic candidate has more bases in a parent that this threshold, the candidate is not marked as
     * STRICT_SOMATIC. A reasonable default is 1.
     */
    private int strictThresholdParents = 0;
    /**
     * If a somatic candidate has more bases in a matched germline sample that this threshold, the candidate is
     * not marked as STRICT_SOMATIC. A reasonable default is 10.
     */
    private int strictThresholdGermline = 10;

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
        recordNumAlignedReads(mode.getInputFilenames());

        igvFieldIndex = statsWriter.defineField("INFO", "BIOMART_COORDS", 1, ColumnType.String, "Coordinates formatted for use with IGV.");
        strictSomaticCandidateFieldIndex = statsWriter.defineField("FILTER", "STRICT_SOMATIC", 1, ColumnType.Flag, "Indicates that the site is not a strict somatic candidate. Strict somatic candidates are not detected in the parents and only poorly in the matched germline. False otherwise.");
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
        candidateFrequencyIndex = new int[numSamples];
        maxGenotypeSomaticPriority = new int[numSamples];
        Arrays.fill(somaticPValueIndex, -1);

        for (String sample : somaticSampleIds) {
            int sampleIndex = locateSampleIndex(sample);
            assert sampleIndex != -1 : "sample-id must match between covariate file and alignment basenames.";
            candidateFrequencyIndex[sampleIndex] = statsWriter.defineField("INFO",
                    String.format("somatic-frequency[%s]", sample),
                    1, ColumnType.Float,
                    "Frequency of a somatic variation (%), valid only when the p-value is significant.", "statistic", "indexed");
            maxGenotypeSomaticPriority[sampleIndex] = statsWriter.defineField("INFO",
                    String.format("priority[%s]", sample),
                    1, ColumnType.Float,
                    "Somatic priority, larger numbers indicate more support for somatic variation in sample (%)", "statistic", "indexed");

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
                        LOG.warn("Parent could not be found for id:" + parentId);
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
                    assert germlineSampleIndex != -1 : "A sampleId in the covariate file (" + sampleId + ") does not match the input alignment basenames.";
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

    int[] numMatchedReads;

    private void recordNumAlignedReads(String[] inputFilenames) {
        numMatchedReads = new int[samples.length];
        String[] sampleIds = AlignmentReaderImpl.getBasenames(inputFilenames);
        // also remove the path to the file to keep only filenames:
        for (int i = 0; i < sampleIds.length; i++) {
            sampleIds[i] = FilenameUtils.getName(sampleIds[i]
            );
            if (sampleIds[i].equals(samples[i])) {
                // put a minimum of one read to prevent divisions by zero:
                numMatchedReads[i] = Math.max(1, getNumMatchedReads(inputFilenames[i]));
            }
        }
    }

    private int getNumMatchedReads(String inputFilename) {
        AlignmentReaderImpl reader = null;
        try {
            reader = new AlignmentReaderImpl(inputFilename);
            reader.readHeader();
            return reader.getNumberOfAlignedReads();
        } catch (IOException e) {

            LOG.error("Cannot read header file for alignment " + inputFilename, e);
            return 0;
        } finally {
            if (reader != null) {
                reader.close();
            }
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
                        position)
        );
        statsWriter.setChromosome(currentReferenceId);

        statsWriter.setPosition(position);

        genotypeFormatter.writeGenotypes(statsWriter, sampleCounts, position);
        if (!statsWriter.hasAlternateAllele()) {
            // do not write a record if the position does not have an alternate allele.
            return;
        }
        allocateIsSomaticCandidate(sampleCounts);

        // Do not write record if alleleSet is empty, IGV VCF track cannot handle that.
        if (isPossibleSomaticVariation(sampleCounts)) {
            estimateSomaticFrequencies(sampleCounts);
            estimatePriority(sampleCounts);
            if (isSomaticCandidate()) {
                statsWriter.writeRecord();
            }
        }

        updateSampleCumulativeCounts(sampleCounts);
    }

    void allocateIsSomaticCandidate(SampleCountInfo[] sampleCounts) {
        int maxGenotypeIndex = 0;
        for (SampleCountInfo sci : sampleCounts) {
            maxGenotypeIndex = Math.max(maxGenotypeIndex, sci.getGenotypeMaxIndex());
        }
        isSomaticCandidate = new boolean[sampleCounts.length][maxGenotypeIndex];
        isStrictSomaticCandidate = new boolean[sampleCounts.length][maxGenotypeIndex];
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

    public void estimatePriority(SampleCountInfo[] sampleCounts) {

        for (int sampleIndex : somaticSampleIndices) {
            SampleCountInfo somaticCounts = sampleCounts[sampleIndex];
            double maxPriority = -10;

            for (int genotypeIndex = 0; genotypeIndex < somaticCounts.getGenotypeMaxIndex(); ++genotypeIndex) {
                if (isSomaticCandidate[sampleIndex][genotypeIndex]) {
                    double parentGenotypePriorityContribution = Double.MAX_VALUE;
                    double germlineGenotypePriorityContribution = Double.MAX_VALUE;
                    int fatherSampleIndex = sample2FatherSampleIndex[sampleIndex];
                    int numParents = 0;
                    if (fatherSampleIndex != -1) {

                        SampleCountInfo fatherCounts = sampleCounts[fatherSampleIndex];
                        double fatherPriority = estimatePriorityComponent(genotypeIndex, somaticCounts, fatherCounts, sampleIndex, fatherSampleIndex);
                        parentGenotypePriorityContribution = Math.min(parentGenotypePriorityContribution, fatherPriority);
                        numParents += 1;
                    }

                    int motherSampleIndex = sample2MotherSampleIndex[sampleIndex];
                    if (motherSampleIndex != -1) {

                        SampleCountInfo motherCounts = sampleCounts[motherSampleIndex];
                        double motherPriority = estimatePriorityComponent(genotypeIndex, somaticCounts, motherCounts, sampleIndex, motherSampleIndex);
                        parentGenotypePriorityContribution = Math.min(parentGenotypePriorityContribution, motherPriority);
                        numParents += 1;
                    }
                    if (numParents == 0) {
                        parentGenotypePriorityContribution = 0;
                    }
                    //     parentGenotypePriorityContribution /= Math.max(1, numParents);
                    int germlineSampleIndices[] = sample2GermlineSampleIndices[sampleIndex];
                    int numGermlineSamples = 0;
                    for (int germlineSampleIndex : germlineSampleIndices) {
                        if (germlineSampleIndex != -1) {
                            SampleCountInfo germlineCounts = sampleCounts[germlineSampleIndex];
                            double germlinePriority = estimatePriorityComponent(genotypeIndex, somaticCounts, germlineCounts, sampleIndex, germlineSampleIndex);
                            germlineGenotypePriorityContribution = Math.min(germlineGenotypePriorityContribution, germlinePriority);
                            numGermlineSamples += 1;
                        }
                    }
                    if (numGermlineSamples == 0) {
                        germlineGenotypePriorityContribution = 0;
                    }
                    maxPriority = Math.max(parentGenotypePriorityContribution + germlineGenotypePriorityContribution, maxPriority);
                }
            }
            statsWriter.setInfo(maxGenotypeSomaticPriority[sampleIndex], maxPriority);
            statsWriter.setFilter(isStrictSomaticCandidate() ? "PASS" : "STRICT_SOMATIC");
        }
    }

    /**
     * Calculate a priority score normalized by the total number of reads mapped in a sample. Mulitply by 100 million
     * to get scores in a more convenient range.
     *
     * @param priority    Priority score contribution calculated for one sample
     * @param sampleIndex sample index of the corresponding sample.
     * @return
     */
    private double normalizePriority(int priority, int sampleIndex) {
        return 100000000d * ((double) priority) / (Math.max(1, numMatchedReads[sampleIndex]));
    }

    private double estimatePriorityComponent(int genotypeIndex, SampleCountInfo somaticCounts,
                                             SampleCountInfo parentOrGermlineCounts,
                                             int somaticSampleIndex, int germlineSampleIndex) {

        double somaticNormalized = normalizePriority(somaticCounts.getGenotypeCount(genotypeIndex), somaticSampleIndex);
        double germlineNormalized = normalizePriority(parentOrGermlineCounts.getGenotypeCount(genotypeIndex), germlineSampleIndex);

        return somaticNormalized - germlineNormalized;
    }


    boolean isPossibleSomaticVariation(SampleCountInfo[] sampleCounts) {

        // In cases where both parents are homozygous and the patient can be heterozygous, which creates low fisher p-values
        // in the contingency table of base counts.

        // We estimate the frequency of each base detected in the somatic patient sample. We calculate the frequency
        // of this same base in the father or mother and take the maximum parent frequency.

        // if the frequency of any base in the somatic sample is larger than the parent frequency, we output the
        // record. The p-value will inform about the strength of the somatic observation.
        // otherwise, we do not output the variation in the somatic report.

        for (int sampleIndex : somaticSampleIndices) {
            SampleCountInfo somaticCounts = sampleCounts[sampleIndex];
            for (int genotypeIndex = 0; genotypeIndex < somaticCounts.getGenotypeMaxIndex(); genotypeIndex++) {
                boolean parentHasGenotype = false;
                boolean strict = true;
                float maxGermlineOrParentsFrequency = 0;
                int fatherSampleIndex = sample2FatherSampleIndex[sampleIndex];
                int minGermlineCoverage = Integer.MAX_VALUE;

                if (fatherSampleIndex != -1) {

                    SampleCountInfo fatherCounts = sampleCounts[fatherSampleIndex];
                    minGermlineCoverage = Math.min(fatherCounts.coverage(), minGermlineCoverage);
                    int fatherCount = fatherCounts.getGenotypeCount(genotypeIndex);
                    parentHasGenotype = fatherCount > fatherCounts.failedCount || fatherCount > 5;
                    strict = fatherCount <= strictThresholdParents;
                    maxGermlineOrParentsFrequency = Math.max(maxGermlineOrParentsFrequency, fatherCounts.frequency(genotypeIndex));
                }
                int motherSampleIndex = sample2MotherSampleIndex[sampleIndex];
                if (motherSampleIndex != -1) {

                    SampleCountInfo motherCounts = sampleCounts[motherSampleIndex];
                    minGermlineCoverage = Math.min(motherCounts.coverage(), minGermlineCoverage);
                    int motherCount = motherCounts.getGenotypeCount(genotypeIndex);
                    parentHasGenotype |= motherCount > motherCounts.failedCount || motherCount > 5;
                    strict &= motherCount <= strictThresholdParents;
                    maxGermlineOrParentsFrequency = Math.max(maxGermlineOrParentsFrequency, motherCounts.frequency(genotypeIndex));

                }
                boolean germlineHasPhenotype = false;
                int germlineSampleIndices[] = sample2GermlineSampleIndices[sampleIndex];
                for (int germlineSampleIndex : germlineSampleIndices) {
                    if (germlineSampleIndex != -1) {
                        SampleCountInfo germlineCounts = sampleCounts[germlineSampleIndex];
                        minGermlineCoverage = Math.min(germlineCounts.coverage(), minGermlineCoverage);
                        int germlineCount = germlineCounts.getGenotypeCount(genotypeIndex);
                        germlineHasPhenotype |= germlineCount >= 10 && germlineCount >= 1.5 * somaticCounts.failedCount;
                        strict &= germlineCount <= strictThresholdGermline;
                        maxGermlineOrParentsFrequency = Math.max(maxGermlineOrParentsFrequency, germlineCounts.frequency(genotypeIndex));
                    }
                }
                boolean somaticHasGenotype = somaticCounts.getGenotypeCount(genotypeIndex) > 0;
                if (parentHasGenotype || germlineHasPhenotype || !somaticHasGenotype) {
                    isSomaticCandidate[sampleIndex][genotypeIndex] = false;
                    isStrictSomaticCandidate[sampleIndex][genotypeIndex] = false;
                    //   System.out.println(explainSomaticCandidateChoice(sampleCounts));
                } else {
                    int somaticCoverage = sampleCounts[sampleIndex].coverage();
                    if (minGermlineCoverage < somaticCoverage / 2) {
                        // not enough coverage in germline samples to call this site confidently
                        isSomaticCandidate[sampleIndex][genotypeIndex] = false;
                    } else {
                        if (somaticCounts.frequency(genotypeIndex) > 3 * maxGermlineOrParentsFrequency) {

                            isSomaticCandidate[sampleIndex][genotypeIndex] = true;
                        }
                        isStrictSomaticCandidate[sampleIndex][genotypeIndex] = strict && isSomaticCandidate[sampleIndex][genotypeIndex];
                    }
                }
            }
        }
        //System.out.println(explainSomaticCandidateChoice(sampleCounts));
        return isSomaticCandidate();
    }

    private String explainSomaticCandidateChoice(SampleCountInfo[] sampleCounts) {
        int sampleIndex = 0;
        int genotypeIndex = 0;
        String output = "";

        for (boolean[] someGenotypeIsSomatic : isSomaticCandidate) {
            genotypeIndex = 0;
            for (boolean candidate : someGenotypeIsSomatic) {
                if (candidate) {
                    output += "genotype " + sampleCounts[sampleIndex].baseString(genotypeIndex) + " is candidate somatic.\n";
                }
                genotypeIndex++;
            }
            sampleIndex++;

        }
        return output;
    }

    public void estimateSomaticFrequencies(SampleCountInfo[] sampleCounts) {
        // force recalculation of the isSomaticCandidate arrays:
        isPossibleSomaticVariation(sampleCounts);

        for (int sampleIndex : somaticSampleIndices) {

            SampleCountInfo somaticCounts = sampleCounts[sampleIndex];

            float somaticFrequency = 0;
            for (int genotypeIndex = 0; genotypeIndex < somaticCounts.getGenotypeMaxIndex(); ++genotypeIndex) {
                if (isSomaticCandidate[sampleIndex][genotypeIndex]) {
                    somaticFrequency = Math.max(somaticCounts.frequency(genotypeIndex), somaticFrequency);
                }
            }
            if (!isSomaticCandidate()) {
                somaticFrequency = 0;
            }
            statsWriter.setInfo(candidateFrequencyIndex[sampleIndex], somaticFrequency * 100);
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

    boolean isSomaticCandidate[][];
    boolean isStrictSomaticCandidate[][];

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

    public boolean isSomaticCandidate() {
        for (boolean[] someGenotypeIsSomatic : isSomaticCandidate) {
            for (boolean candidate : someGenotypeIsSomatic) {
                if (candidate) return true;
            }
        }
        return false;
    }

    public boolean isStrictSomaticCandidate() {
        for (boolean[] someGenotypeIsSomatic : isStrictSomaticCandidate) {
            for (boolean candidate : someGenotypeIsSomatic) {
                if (candidate) return true;
            }
        }
        return false;
    }

    public void setCandidateFrequencyIndex(int[] candidateFrequencyIndex) {
        this.candidateFrequencyIndex = candidateFrequencyIndex;
    }

    public int[] getCandidateFrequencyIndex() {
        return candidateFrequencyIndex;
    }


    public void estimateSomaticPValue(SampleCountInfo[] sampleCounts) {
        // this method does nothing. Kept for compatibility with JUnit tests.
    }

    public void setupR() {
        // this method does nothing. Kept for compatibility with JUnit tests.
    }
}
