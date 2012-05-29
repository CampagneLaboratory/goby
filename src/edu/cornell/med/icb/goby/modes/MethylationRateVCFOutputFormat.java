/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.cornell.med.icb.goby.modes;

import edu.cornell.med.icb.goby.R.GobyRengine;
import edu.cornell.med.icb.goby.algorithmic.algorithm.dmr.EstimatedDistribution;
import edu.cornell.med.icb.goby.algorithmic.algorithm.dmr.SitesInFixedWindow;
import edu.cornell.med.icb.goby.algorithmic.data.GroupComparison;
import edu.cornell.med.icb.goby.algorithmic.data.MethylCountInfo;
import edu.cornell.med.icb.goby.algorithmic.data.SamplePairEnumerator;
import edu.cornell.med.icb.goby.alignments.*;
import edu.cornell.med.icb.goby.readers.vcf.ColumnType;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import edu.cornell.med.icb.goby.stats.*;
import edu.cornell.med.icb.goby.util.dynoptions.DynamicOptionClient;
import edu.cornell.med.icb.goby.util.OutputInfo;
import edu.cornell.med.icb.goby.util.dynoptions.RegisterThis;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.apache.commons.io.output.NullWriter;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.rosuda.JRI.Rengine;

import java.io.IOException;
import java.util.ArrayList;

/**
 * A Variant Call Format output to estimate methylation rates for a set of samples and find methylation rate differences
 * between group.
 *
 * @author Fabien Campagne
 *         Date: April 4 2011
 *         Time: 2:38:13 AM
 */
public class MethylationRateVCFOutputFormat extends AbstractOutputFormat implements SequenceVariationOutputFormat, MethylationFormat {

    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(MethylationRateVCFOutputFormat.class);
    @RegisterThis
    public static final DynamicOptionClient doc = new DynamicOptionClient(MethylationRateVCFOutputFormat.class,
            EmpiricalPValueEstimator.LOCAL_DYNAMIC_OPTIONS,
            "window-length: int, length of the fixed genomic window used to detect DMR. Default 1000 bp:1000",
            "significance-threshold: double, significance threshold to consider a site significant for DMR detection purposes. Default 0.01:0.01"
    );
    private boolean writeFieldGroupAssociations=true;

    public static DynamicOptionClient doc() {
        return doc;
    }

    VCFWriter statWriter;
    private int refIdColumnIndex;
    private int positionColumnIndex;

    private boolean fisherRInstalled;
    private String[] groups;
    private String[] samples;
    private ObjectArrayList<ReadIndexStats> readIndexStats;
    private int log2OddsRatioColumnIndex[];
    private int fisherExactPValueColumnIndex[];
    private int numberOfGroups;

    private DifferentialExpressionAnalysis deAnalyzer;
    private DifferentialExpressionCalculator deCalculator;

    private int log2OddsRatioStandardErrorColumnIndex[];
    private int log2OddsRatioZColumnIndex[];
    private int deltaMRColumnIndex[];
    int[] readerIndexToGroupIndex;
    private SitesInFixedWindow[] fixedWindow;
    private IntSet[] distinctReadIndicesCountPerGroup;

    private MethylCountInfo mci;
    private int numberOfSamples;
    private int biomartFieldIndex;
    private GenotypesOutputFormat genotypeFormatter;
    private int depthFieldIndex;
    private int[] methylatedCCountsIndex;
    private int[] notMethylatedCCountsIndex;

    private int methylationRateFieldIndex;
    private int[] methylationRatePerGroupIndex;

    private int strandFieldIndex;
    private ArrayList<GroupComparison> groupComparisons = new ArrayList<GroupComparison>();
    private int convertedCytosineFieldIndex;
    private int unconvertedCytosineFieldIndex;
    private int convertedCytosinePlusFieldIndex;
    private int unconvertedCytosinePlusFieldIndex;
    private int convertedCytosineMinusFieldIndex;
    private int unconvertedCytosineMinusFieldIndex;
    private RandomAccessSequenceInterface genome;
    private int genomicContextIndex;
    private CharSequence chromosome;
    private int genomeReferenceIndex;
    private int[] empiricalPValueColumnIndex;
    private boolean estimateIntraGroupDifferences;
    private boolean estimateEmpiricalP;
    private OutputInfo outputInfo;
    private int[] lastOfDMRIndex;
    private double fixedWindowEmpiricalPSignificanceThreshold;

    @Override
    public void setMinimumEventThreshold(final int minimumEventThreshold) {
        this.minimumEventThreshold = minimumEventThreshold;
    }

    private int minimumEventThreshold = 1;


    public void defineColumns(final OutputInfo outputInfo, final DiscoverSequenceVariantsMode mode) {
        deAnalyzer = mode.getDiffExpAnalyzer();
        deCalculator = mode.getDiffExpCalculator();
        groups = mode.getGroups();
        samples = mode.getSamples();

        this.outputInfo = outputInfo;
        estimateIntraGroupDifferences = doc.getBoolean("estimate-intra-group-differences");
        estimateEmpiricalP = doc.getBoolean("estimate-empirical-P");
        fixedWindowEmpiricalPSignificanceThreshold = doc.getDouble("significance-threshold");
        final int windowLength = doc.getInteger("window-length");


        readerIndexToGroupIndex = mode.getReaderIndexToGroupIndex();
        final ObjectArrayList<ReadIndexStats> readIndexStats = mode.getReadIndexStats();
        final VCFWriter vcfWriter = new VCFWriter(estimateIntraGroupDifferences ? new NullWriter() : outputInfo.getPrintWriter());
        this.statWriter = vcfWriter;
        groupComparisons = mode.getGroupComparisons();
        int maxComparisonIndex = -1;
        for (final GroupComparison comparison : groupComparisons) {
            maxComparisonIndex = Math.max(comparison.index, maxComparisonIndex);
        }
        fixedWindow = new SitesInFixedWindow[maxComparisonIndex + 1];
        for (final GroupComparison comparison : groupComparisons) {
            fixedWindow[comparison.index] = new SitesInFixedWindow(windowLength);
        }

        try {
            //activate R only if we need it:
            final Rengine rEngine = GobyRengine.getInstance().getRengine();
            fisherRInstalled = rEngine != null && rEngine.isAlive();
        } catch (java.lang.UnsatisfiedLinkError e) {
            System.out.println("Cannot initialize R");
            e.printStackTrace();
            throw e;
        }
        if (groups.length < 1) {
            System.err.println("Methylation format requires at least one group.");
            System.exit(1);
        }

        this.readIndexStats = readIndexStats;

        log2OddsRatioColumnIndex = new int[groupComparisons.size()];
        log2OddsRatioStandardErrorColumnIndex = new int[groupComparisons.size()];
        log2OddsRatioZColumnIndex = new int[groupComparisons.size()];
        fisherExactPValueColumnIndex = new int[groupComparisons.size()];
        empiricalPValueColumnIndex = new int[groupComparisons.size()];
        deltaMRColumnIndex = new int[groupComparisons.size()];
        lastOfDMRIndex = new int[groupComparisons.size()];

        numberOfGroups = groups.length;
        statWriter.setWriteFieldGroupAssociations(writeFieldGroupAssociations);
        biomartFieldIndex = statWriter.defineField("INFO", "BIOMART_COORDS", 1, ColumnType.String, "Coordinates for use with Biomart.","biomart");
        strandFieldIndex = statWriter.defineField("INFO", "Strand", 1, ColumnType.String, "Strand of the cytosine site on the reference sequence.", "genomic-coordinate");
        genomicContextIndex = statWriter.defineField("INFO", "Context", 1, ColumnType.String, "Site genomic context");

        for (GroupComparison comparison : groupComparisons) {
            log2OddsRatioColumnIndex[comparison.index] = statWriter.defineField("INFO", String.format("LOD[%s/%s]", comparison.nameGroup1, comparison.nameGroup2),
                    1, ColumnType.Float, String.format("Log2 of the odds-ratio of observing methylation in  group %s versus group %s", comparison.nameGroup1, comparison.nameGroup2),
                    "statistic");

            log2OddsRatioStandardErrorColumnIndex[comparison.index] = statWriter.defineField("INFO", String.format("LOD_SE[%s/%s]", comparison.nameGroup1, comparison.nameGroup2),
                    1, ColumnType.Float, String.format("Standard Error of the log2 of the odds-ratio between group %s and group %s", comparison.nameGroup1, comparison.nameGroup2)
                    ,"LOD","statistic");

            log2OddsRatioZColumnIndex[comparison.index] = statWriter.defineField("INFO", String.format("LOD_Z[%s/%s]", comparison.nameGroup1, comparison.nameGroup2),
                    1, ColumnType.Float, String.format("Z value of the odds-ratio between group %s and group %s", comparison.nameGroup1, comparison.nameGroup2),
                    "LOD","statistic");

            fisherExactPValueColumnIndex[comparison.index] = statWriter.defineField("INFO", String.format("FisherP[%s/%s]", comparison.nameGroup1, comparison.nameGroup2),
                    1, ColumnType.Float,
                    String.format("Fisher exact P-value of observing as large a difference by chance between group %s and group %s.", comparison.nameGroup1, comparison.nameGroup2),"p-value");

            deltaMRColumnIndex[comparison.index] = statWriter.defineField("INFO", String.format("Delta_MR[%s/%s]", comparison.nameGroup1, comparison.nameGroup2),
                    1, ColumnType.Integer, String.format("Absolute Difference in methylation between group %s and group %s", comparison.nameGroup1, comparison.nameGroup2),
                    "statistic");
            if (estimateEmpiricalP) {
                empiricalPValueColumnIndex[comparison.index] = statWriter.defineField("INFO", String.format("empiricalP[%s/%s]", comparison.nameGroup1, comparison.nameGroup2),
                        1, ColumnType.Float,
                        String.format("Empirical P-value of observing as large a difference by chance between group %s and group %s.", comparison.nameGroup1, comparison.nameGroup2),"p-value","statistic");

                lastOfDMRIndex[comparison.index] = statWriter.defineField("INFO", String.format("DMR[%s/%s]", comparison.nameGroup1, comparison.nameGroup2),
                        1, ColumnType.Integer, String.format("Number of significant sites observed in a fixed window of length %d before this significant site (inclusive). Comparing group %s and group %s.",
                        windowLength, comparison.nameGroup1, comparison.nameGroup2),"DMR");
            }

        }
        methylationRatePerGroupIndex = new int[groups.length];
        for (int groupIndex = 0; groupIndex < groups.length; groupIndex++) {
            methylationRatePerGroupIndex[groupIndex] = statWriter.defineField("INFO", String.format("MR[%s]", groups[groupIndex]),
                    1, ColumnType.Float, String.format("Methylation rate in group %s.", groups[groupIndex]),
                    "MethylationRate");

        }
        methylatedCCountsIndex = new int[groups.length];
        notMethylatedCCountsIndex = new int[groups.length];

        for (int groupIndex = 0; groupIndex < groups.length; groupIndex++) {
            notMethylatedCCountsIndex[groupIndex] = statWriter.defineField("INFO", String.format("#C_Group[%s]", groups[groupIndex]),
                    1, ColumnType.Integer, String.format("Number of unmethylated Cytosines at site in group %s.", groups[groupIndex]),"#C");
            methylatedCCountsIndex[groupIndex] = statWriter.defineField("INFO", String.format("#Cm_Group[%s]", groups[groupIndex]),
                    1, ColumnType.Integer, String.format("Number of methylated Cytosines at site in group %s.", groups[groupIndex]),"#Cm");
        }
        depthFieldIndex = statWriter.defineField("INFO", "DP",
                1, ColumnType.Integer, "Total depth of sequencing across groups at this site","Depth");
        {
            // define one VCF 'sample' for each strand of each input sample:
            String[] sampleTwice = new String[samples.length * 2];
            int sampleIndex = 0;
            for (String sample : samples) {
                sampleTwice[sampleIndex++] = "+|" + sample;
            }
            for (String sample : samples) {
                sampleTwice[sampleIndex++] = "-|" + sample;
            }

            statWriter.defineSamples(sampleTwice);
        }
        // define Genotype as first field (required by VCF specification):
        genotypeFormatter.defineGenotypeField(statWriter);
        methylationRateFieldIndex = statWriter.defineField("FORMAT", "MR", 1, ColumnType.Integer, "Methylation rate. 0-100%, 100% indicate fully methylated.",
                "MethylationRate");
        convertedCytosineFieldIndex = statWriter.defineField("FORMAT", "C", 1, ColumnType.Integer, "Number of converted cytosines at site.","#C");
        unconvertedCytosineFieldIndex = statWriter.defineField("FORMAT", "Cm", 1, ColumnType.Integer, "Number of unconverted cytosines at site", "#Cm");

        // setup empirical p-value estimator:
        String statName = doc.getString("statistic");
        if (!statName.endsWith("_mci")) {

            // do not append several times
            statName = statName + "_mci";
        }
        // override the stat name with the mci provider variant:
        doc.setValue("statistic", statName);
        empiricalPValueEstimator = new EmpiricalPValueEstimator();
        empiricalPValueEstimator.configure(1, doc);
        empiricalPValueEstimator.setGroupEnumerator(new SamplePairEnumerator(readerIndexToGroupIndex, samples.length, groups.length, groupComparisons.size()));

        if (estimateIntraGroupDifferences || estimateEmpiricalP) {
            empiricalPValueEstimator.recordWithinGroupSamplePairs(groups);
            for (final GroupComparison comparison : groupComparisons) {

                empiricalPValueEstimator.recordBetweenGroupsSamplePairs(comparison);
            }
        }
        statWriter.writeHeader();
    }

    @Override
    public void allocateStorage(final int numberOfSamples, final int numberOfGroups) {
        this.numberOfGroups = numberOfGroups;
        this.numberOfSamples = numberOfSamples;
        mci = new MethylCountInfo(numberOfSamples, numberOfGroups);
        distinctReadIndicesCountPerGroup = new IntSet[numberOfGroups];
        for (int i = 0; i < numberOfGroups; i++) {
            distinctReadIndicesCountPerGroup[i] = new IntArraySet();
        }
        genotypeFormatter = new GenotypesOutputFormat();
        genotypeFormatter.allocateStorage(numberOfSamples, numberOfGroups);
    }

    @Override
    public void writeRecord(final DiscoverVariantIterateSortedAlignments iterator,
                            final SampleCountInfo[] sampleCounts,
                            final int referenceIndex,
                            int position,
                            final DiscoverVariantPositionData list,
                            final int groupIndexA,
                            final int groupIndexB) {

        final int oneBasedPosition = position + 1;


        final char refBase = sampleCounts[0].referenceBase;
        if (refBase != 'C' && refBase != 'G') {
            mci.reset();
            return;
        }
        fillMethylationCountArrays(sampleCounts, list, oneBasedPosition, refBase, mci, readerIndexToGroupIndex);
        doEmpiricalP(mci, referenceIndex, position);
        if (estimateIntraGroupDifferences) {
            return;
        }
        if (mci.eventCountAtSite < minimumEventThreshold) {
            return;
        }
        statWriter.setInfo(depthFieldIndex, list.size());
        final CharSequence currentReferenceId = iterator.getReferenceId(referenceIndex);

        statWriter.setChromosome(currentReferenceId);
        statWriter.setPosition(oneBasedPosition);
        statWriter.setReferenceAllele(Character.toString(sampleCounts[0].referenceBase));

        // construct a biomart region span in the format chr:pos1:chr:pos
        final String biomartRegionSpan = String.format("%s:%s:%s", currentReferenceId, oneBasedPosition,
                oneBasedPosition);

        statWriter.setInfo(biomartFieldIndex, biomartRegionSpan);
        statWriter.setInfo(strandFieldIndex, Character.toString(mci.strandAtSite));

        final String genomicContext = findGenomicContext(genomeReferenceIndex, oneBasedPosition);
        statWriter.setInfo(genomicContextIndex, genomicContext);

        for (int groupIndex = 0; groupIndex < numberOfGroups; groupIndex++) {
            statWriter.setInfo(notMethylatedCCountsIndex[groupIndex], mci.unmethylatedCCountPerGroup[groupIndex]);
            statWriter.setInfo(methylatedCCountsIndex[groupIndex], mci.methylatedCCountPerGroup[groupIndex]);
            statWriter.setInfo(methylationRatePerGroupIndex[groupIndex], getMR(groupIndex));
        }

        for (int sampleIndex = 0; sampleIndex < numberOfSamples; sampleIndex++) {
            final float numerator = mci.methylatedCCountPerSample[sampleIndex];
            final float denominator = numerator + mci.unmethylatedCCountPerSample[sampleIndex];

            final float methylationRate = numerator * 100 / denominator;
            statWriter.setSampleValue(methylationRateFieldIndex, sampleIndex, Math.round(methylationRate));
            statWriter.setSampleValue(unconvertedCytosineFieldIndex, sampleIndex, mci.methylatedCCountPerSample[sampleIndex]);
            statWriter.setSampleValue(convertedCytosineFieldIndex, sampleIndex, mci.unmethylatedCCountPerSample[sampleIndex]);
        }
        for (final GroupComparison comparison : groupComparisons) {
            final int indexGroup1 = comparison.indexGroup1;
            final int indexGroup2 = comparison.indexGroup2;
            final double denominator = (double) (mci.unmethylatedCCountPerGroup[indexGroup1]) * (double) (mci.methylatedCCountPerGroup[indexGroup2]);
            final double oddsRatio = denominator == 0 ? Double.NaN :
                    ((double) (mci.unmethylatedCCountPerGroup[indexGroup2]) * (double) (mci.methylatedCCountPerGroup[indexGroup1])) /
                            denominator;
            final double logOddsRatioSE;

            if (mci.methylatedCCountPerGroup[indexGroup1] < 10 ||
                    mci.methylatedCCountPerGroup[indexGroup2] < 10 ||
                    mci.unmethylatedCCountPerGroup[indexGroup1] < 10 ||
                    mci.unmethylatedCCountPerGroup[indexGroup2] < 10) {
                // standard error estimation is unreliable when any of the counts are less than 10.
                logOddsRatioSE = Double.NaN;
            } else {
                logOddsRatioSE = Math.sqrt(1d / mci.unmethylatedCCountPerGroup[indexGroup2] +
                        1d / mci.methylatedCCountPerGroup[indexGroup1] +
                        1d / mci.methylatedCCountPerGroup[indexGroup2] +
                        1d / mci.unmethylatedCCountPerGroup[indexGroup1]);
            }
            final double log2OddsRatio = Math.log(oddsRatio) / Math.log(2);
            final double log2OddsRatioZValue = log2OddsRatio / logOddsRatioSE;
            double fisherP = Double.NaN;
            if (mci.eventCountAtSite >= minimumEventThreshold) {
                // estimate Fisher only if we have seen at least 10 events.

                final boolean ok = checkCounts();
                if (ok) {
                    fisherP = fisherRInstalled ? FisherExactRCalculator.getFisherPValue(
                            mci.unmethylatedCCountPerGroup[indexGroup2], mci.methylatedCCountPerGroup[indexGroup2],
                            mci.unmethylatedCCountPerGroup[indexGroup1], mci.methylatedCCountPerGroup[indexGroup1]) : Double.NaN;
                } else {
                    System.err.printf("An exception was caught evaluating the Fisher Exact test P-value. Details are provided below%n" +
                            "referenceId=%s referenceIndex=%d position=%d %n" +
                            "unmethylatedCCountsPerGroup[1]=%d methylatedCCountPerGroup[1]=%d%n" +
                            "unmethylatedCCountsPerGroup[0]=%d, methylatedCCountPerGroup[0]=%d",
                            currentReferenceId, referenceIndex,
                            position,
                            mci.unmethylatedCCountPerGroup[indexGroup2], mci.methylatedCCountPerGroup[indexGroup2],
                            mci.unmethylatedCCountPerGroup[indexGroup1], mci.methylatedCCountPerGroup[indexGroup1]
                    );
                }
            }

            final double group1MR = getMR(indexGroup1);
            final double group2MR = getMR(indexGroup2);
            final int deltaMR = (int) Math.round(Math.abs(group1MR - group2MR));

            statWriter.setInfo(log2OddsRatioColumnIndex[comparison.index], log2OddsRatio);
            statWriter.setInfo(log2OddsRatioStandardErrorColumnIndex[comparison.index], logOddsRatioSE);
            statWriter.setInfo(log2OddsRatioZColumnIndex[comparison.index], log2OddsRatioZValue);
            statWriter.setInfo(fisherExactPValueColumnIndex[comparison.index], fisherP);
            statWriter.setInfo(deltaMRColumnIndex[comparison.index], deltaMR);

        }
        genotypeFormatter.writeGenotypes(statWriter, sampleCounts, oneBasedPosition);
        for (int sampleIndex = 0; sampleIndex < numberOfSamples; sampleIndex++) {
            final int firstIndex = sampleIndex;
            final int secondIndex = convertIndex(sampleIndex, mci.strandAtSite);
            if (mci.strandAtSite == '+') {
                final int otherIndex = convertIndex(sampleIndex, '-');
                statWriter.setSampleValue(methylationRateFieldIndex, otherIndex, "100");
                statWriter.setSampleValue(convertedCytosineFieldIndex, otherIndex, "0");
                statWriter.setSampleValue(unconvertedCytosineFieldIndex, otherIndex, "0");
                statWriter.setSampleValue(genotypeFormatter.getGenotypeFieldIndex(), otherIndex, "0/0");
                statWriter.setSampleValue(genotypeFormatter.getBaseCountFieldIndex(), otherIndex, "ignore");
                statWriter.setSampleValue(genotypeFormatter.getFailBaseCountFieldIndex(), otherIndex, "ignore");
                statWriter.setSampleValue(genotypeFormatter.getGoodBaseCountFieldIndex(), otherIndex, "0");
            } else {
                statWriter.switchSampleValue(methylationRateFieldIndex, firstIndex, secondIndex, "100");
                statWriter.switchSampleValue(convertedCytosineFieldIndex, firstIndex, secondIndex, "0");
                statWriter.switchSampleValue(unconvertedCytosineFieldIndex, firstIndex, secondIndex, "0");
                statWriter.switchSampleValue(genotypeFormatter.getGenotypeFieldIndex(), firstIndex, secondIndex, "0/0");
                statWriter.switchSampleValue(genotypeFormatter.getBaseCountFieldIndex(), firstIndex, secondIndex, "ignore");
                statWriter.switchSampleValue(genotypeFormatter.getFailBaseCountFieldIndex(), firstIndex, secondIndex, "ignore");
                statWriter.switchSampleValue(genotypeFormatter.getGoodBaseCountFieldIndex(), firstIndex, secondIndex, "0");
            }
        }
        statWriter.writeRecord();

    }

    // calculate methylation rate in a group:
    private float getMR(int groupIndex) {
        float groupMR;
        final float totalCsObservedgroup1 = mci.methylatedCCountPerGroup[groupIndex] + mci.unmethylatedCCountPerGroup[groupIndex];
        if (totalCsObservedgroup1 == 0) {
            groupMR = Float.NaN;
        } else {
            groupMR = (float) mci.methylatedCCountPerGroup[groupIndex] * 100f / totalCsObservedgroup1;
        }
        return groupMR;
    }

    EmpiricalPValueEstimator empiricalPValueEstimator;

    private void doEmpiricalP(final MethylCountInfo mci, int referenceIndex, int position) {
        if (estimateIntraGroupDifferences) {
            int groupIndex = 0;
            // use all groups available groups to estimate empirical p:
            for (String group : groups) {
                empiricalPValueEstimator.estimateNullDensity(0, groupIndex, mci);
                empiricalPValueEstimator.estimateNullDensity(0, groupIndex, mci);
                groupIndex++;
            }
        }
        if (estimateEmpiricalP) {
            for (final GroupComparison comparison : groupComparisons) {

                final double p = empiricalPValueEstimator.estimateEmpiricalPValue(0, comparison, mci);
                statWriter.setInfo(empiricalPValueColumnIndex[comparison.index], p);
                if (p <= fixedWindowEmpiricalPSignificanceThreshold) {
                    fixedWindow[comparison.index].add(referenceIndex, position);
                } else {
                    fixedWindow[comparison.index].prune(referenceIndex, position);
                }
                // record how many significant sites in the window before this site
                statWriter.setInfo(lastOfDMRIndex[comparison.index], fixedWindow[comparison.index].n());
            }
        }
    }

    private double getDeltaMR(GroupComparison comparison, MethylCountInfo mci) {
        int Cma = mci.methylatedCCountPerGroup[comparison.indexGroup1];
        int Cmb = mci.methylatedCCountPerGroup[comparison.indexGroup2];
        int Ca = mci.unmethylatedCCountPerGroup[comparison.indexGroup1];
        int Cb = mci.unmethylatedCCountPerGroup[comparison.indexGroup2];
        double mra = 100d * ((double) Ca) / ((double) (Cma + Ca));
        double mrb = 100d * ((double) Ca) / ((double) (Cma + Ca));
        return Math.abs(mra - mrb);
    }

    private String findGenomicContext(int referenceIndex, int position) {
        int zeroBasedPos = position - 1;
        char currentBase = genome.get(referenceIndex, zeroBasedPos);
        int referenceLength = genome.getLength(referenceIndex);
        char nextBase = '?';
        String tempContext = new StringBuilder().append('C').append('p').toString();

        char concatBase = '?';

        if (currentBase == 'C') {
            if (referenceLength == position) {
                return Character.toString(currentBase);
            }
            nextBase = genome.get(referenceIndex, (zeroBasedPos + 1));
            concatBase = nextBase;
        } else {
            if (currentBase == 'G') {
                if (zeroBasedPos == 0) {
                    return Character.toString(currentBase);
                }
                nextBase = genome.get(referenceIndex, (zeroBasedPos - 1));
                switch (nextBase) {
                    case 'C':
                        concatBase = 'G';
                        break;
                    case 'A':
                        concatBase = 'T';
                        break;
                    case 'T':
                        concatBase = 'A';
                        break;
                    case 'G':
                        concatBase = 'C';
                        break;
                }
            }
        }
        tempContext = tempContext.concat(Character.toString(concatBase));
        return tempContext;
    }


    private char flipStrand(final char strandAtSite) {
        if (strandAtSite == '+') {
            return '-';
        } else {
            return '+';
        }
    }

    private int convertIndex(final int sampleIndex, final char strandAtSite) {
        final int numSamples = samples.length;
        if (strandAtSite == '+') {
            return sampleIndex;
        } else {
            return sampleIndex + numSamples;
        }
    }


    @Override
    public void close() {
        statWriter.close();
        if (estimateIntraGroupDifferences) {
            // when estimating intra-group differences, we  serialize the estimator to the output.
            try {
                LOG.debug("Storing density to filename: " + outputInfo.getFilename());
                EstimatedDistribution.store(empiricalPValueEstimator.getNullDistribution(), outputInfo.getFilename());
            } catch (IOException e) {
                LOG.error("Unable to write estimator to file", e);
            }
        }
    }

    @Override
    public void setGenome(RandomAccessSequenceInterface genome) {
        this.genome = genome;
    }

    @Override
    public void setGenomeReferenceIndex(int index) {
        genomeReferenceIndex = index;
    }

    public static void fillMethylationCountArrays(final SampleCountInfo[] sampleCounts, final ObjectArrayList<PositionBaseInfo> list,
                                                  final int position, final char refBase, final MethylCountInfo mci,
                                                  final int[] readerIndexToGroupIndex) {

        mci.reset();


        if (refBase == 'C') {
            mci.strandAtSite = '+';
        } else if (refBase == 'G') {
            mci.strandAtSite = '-';
        }

        for (final PositionBaseInfo info : list) {
            //@@@! only look at the strand that the read matches:
            if (refBase == 'G' && info.matchesForwardStrand) {
                continue;
            }
            if (refBase == 'C' && !info.matchesForwardStrand) {
                continue;
            }

            if (refBase == 'C' || refBase == 'G') {
                // readBase is always given in the forward strand..
                final char readBase = info.matchesReference ? refBase : info.to;
                final int sampleIndex = info.readerIndex;
                final int groupIndex = readerIndexToGroupIndex[sampleIndex];

                if (readBase == refBase) {
                    // C staying C on forward strand stayed so because they were either (1) methylated or (2) not converted.
                    ++mci.methylatedCCountPerSample[sampleIndex];
                    ++mci.methylatedCCountPerGroup[groupIndex];
                    ++mci.eventCountAtSite;
                } else {
                    // C became T on forward strand (G->A on reverse) indicates that the Cytosine was not methylated and was converted.
                    ++mci.unmethylatedCCountPerSample[sampleIndex];
                    ++mci.unmethylatedCCountPerGroup[groupIndex];
                    ++mci.eventCountAtSite;
                }

            }
        }

    }

    private boolean checkCounts() {
        boolean ok = true;
        // detect if any count is negative (that's a bug)
        for (final int count : mci.unmethylatedCCountPerGroup) {

            if (count < 0) {
                ok = false;
            }
        }
        for (final int count : mci.methylatedCCountPerGroup) {
            if (count < 0) {
                ok = false;
            }
        }
        return ok;
    }


}
