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

import edu.cornell.med.icb.goby.R.FisherExact;
import edu.cornell.med.icb.goby.R.GobyRengine;
import edu.cornell.med.icb.goby.algorithmic.data.GroupComparison;
import edu.cornell.med.icb.goby.alignments.*;
import edu.cornell.med.icb.goby.readers.vcf.ColumnType;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import edu.cornell.med.icb.goby.stats.DifferentialExpressionAnalysis;
import edu.cornell.med.icb.goby.stats.DifferentialExpressionCalculator;
import edu.cornell.med.icb.goby.stats.VCFWriter;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.rosuda.JRI.Rengine;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * A Variant Call Format output to compare genomic variations across groups. This format implements a fisher exact
 * test that compares the number of samples in each group that have each possible allele. This is the usual allelic
 * test of population genetics. Formally, we consider the matrix:
 * group-1   group-2
 * allele=A   C1A       C2A
 * allele=C
 * allele=G
 * allele=T
 * allele=N   C1N       C2N
 * <p/>
 * and estimate the chance that the number of samples in this N x 2 matrix are non-randomly distributed. The Cij in
 * the matrix are the number of samples in each group i that have at least one count for the j allele. We use N to
 * indicate any non base allele (such as deletion in the reference N='-').
 *
 * @author Fabien Campagne
 *         Date: April 4 2011
 *         Time: 2:38:13 AM
 */
public class CompareGroupsVCFOutputFormat implements SequenceVariationOutputFormat {

    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(CompareGroupsVCFOutputFormat.class);

    private VCFWriter statWriter;
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
    private int[] readerIndexToGroupIndex;
    private int[] refCountsPerGroup;
    private int[] variantsCountPerGroup;
    private int[] distinctReadIndexCountPerGroup;
    private IntSet[] distinctReadIndicesCountPerGroup;
    private float[] averageVariantQualityScorePerGroup;
    private int[] variantsCountPerSample;
    private int[] refCountsPerSample;
    private int numberOfSamples;
    private int biomartFieldIndex;
    private GenotypesOutputFormat genotypeFormatter;
    private int depthFieldIndex;
    private int varCountsIndex[];
    private int refCountsIndex[];
    /**
     * The array used for testing, pointing to pre-allocated or dynamically allocated version. Allele count is defined
     * as the number of times a specific allele is observed in a group of samples (each sample can contribute at most
     * one to the count).
     */
    private int[][] alleleCountsPerGroup;
    // pre-allocated array, when the site has no indels:
    private int[][] alleleCountsPerGroupBaseAllelesOnly;
    // pre-allocated vector, when the site has no indels:

    private int[] vectorBaseAllelesOnly;
    /**
     * The vector used for testing, pointing to pre-allocated or dynamically allocated version:
     */
    private int[] fisherVector;
    /**
     * The number of alleles, for sites that have no indels:
     */
    private final int numberOfAlleles = SampleCountInfo.BASE_MAX_INDEX;
    /**
     * The number of allele observed in group for the allele with most counts in the group.  For instance if counts are
     * group0 group1
     * allele=A    0        0
     * allele=C   10        0
     * allele=T    2        0
     * allele=G    4       15
     * allele=N    3        2
     * <p/>
     * alleleCountGroupMax will contain {10,15}
     */
    private int[] alleleCountGroupMax;

    /**
     * The vector that will be used for testing.
     */
    private long[] fisherVectorUnderTest;
    /**
     * Index of the VCF fields that hold allele counts.
     */
    private int[] alleleCountsGroupIndex;
    private MutableString[] alleleCountsMutableStrings;


    protected void setStatWriter(VCFWriter statWriter) {
        this.statWriter = statWriter;
    }

    ArrayList<GroupComparison> groupComparisons;

    public void defineColumns(final PrintWriter writer, final DiscoverSequenceVariantsMode mode) {
        deAnalyzer = mode.getDiffExpAnalyzer();
        deCalculator = mode.getDiffExpCalculator();
        groups = mode.getGroups();
        samples = mode.getSamples();
        readerIndexToGroupIndex = mode.getReaderIndexToGroupIndex();
        final ObjectArrayList<ReadIndexStats> readIndexStats = mode.getReadIndexStats();
        this.statWriter = new VCFWriter(writer);

        //activate R only if we need it:
        final Rengine rEngine = GobyRengine.getInstance().getRengine();
        fisherRInstalled = rEngine != null && rEngine.isAlive();

        if (groups.length < 1) {
            System.err.println("CompareGroupsVCFOutputFormat requires at least one group.");
            System.exit(1);
        }
        groupComparisons = mode.getGroupComparisons();
        this.readIndexStats = readIndexStats;

        log2OddsRatioColumnIndex = new int[groupComparisons.size()];
        log2OddsRatioStandardErrorColumnIndex = new int[groupComparisons.size()];
        log2OddsRatioZColumnIndex = new int[groupComparisons.size()];
        fisherExactPValueColumnIndex = new int[groupComparisons.size()];

        numberOfGroups = groups.length;
        biomartFieldIndex = statWriter.defineField("INFO", "BIOMART_COORDS", 1, ColumnType.String, "Coordinates for use with Biomart.");

        genotypeFormatter.defineInfoFields(statWriter);

        for (GroupComparison comparison : groupComparisons) {
            log2OddsRatioColumnIndex[comparison.index] = statWriter.defineField("INFO", String.format("LOD[%s/%s]",
                    comparison.nameGroup1, comparison.nameGroup2),
                    1, ColumnType.Float, String.format("Log2 of the odds-ratio of observing a variant in group %s versus group %s",
                    comparison.nameGroup1, comparison.nameGroup2));

            log2OddsRatioStandardErrorColumnIndex[comparison.index] = statWriter.defineField("INFO", String.format("LOD_SE[%s/%s]",
                    comparison.nameGroup1, comparison.nameGroup2),
                    1, ColumnType.Float, String.format("Standard Error of the log2 of the odds-ratio between group %s and group %s",
                    comparison.nameGroup1, comparison.nameGroup2));

            log2OddsRatioZColumnIndex[comparison.index] = statWriter.defineField("INFO", String.format("LOD_Z[%s/%s]",
                    comparison.nameGroup1, comparison.nameGroup2),
                    1, ColumnType.Float, String.format("Z value of the odds-ratio of observing a variant in group %s versus group %s",
                    comparison.nameGroup1, comparison.nameGroup2));
            fisherExactPValueColumnIndex[comparison.index] = statWriter.defineField("INFO", String.format("FisherP[%s/%s]",
                    comparison.nameGroup1, comparison.nameGroup2),
                    1, ColumnType.Float, String.format("Fisher exact P-value that the allelic frequencies (each sample contributes at least one toward the group count of each allele) differ between group %s and group %s.", comparison.nameGroup1, comparison.nameGroup2));
        }
        for (int i = 0; i < numberOfGroups; i++) {
            alleleCountsGroupIndex[i] = statWriter.defineField("INFO", String.format("AlleleCountsInGroup[%s]", groups[i]),
                    1, ColumnType.String, String.format("Count of specific alleles in group %s, in the format (allele=count[,])+. Allele count is defined " +
                    " as the number of times a specific allele is observed in a group of samples (each sample can contribute at most" +
                    " one to the count)", groups[i]));
        }
        varCountsIndex = new int[groups.length];
        refCountsIndex = new int[groups.length];

        for (int groupIndex = 0; groupIndex < groups.length; groupIndex++) {
            refCountsIndex[groupIndex] = statWriter.defineField("INFO", String.format("refCountGroup[%s]", groups[groupIndex]),
                    1, ColumnType.Integer, String.format("Number of reference allele called in group %s.", groups[groupIndex]));
            varCountsIndex[groupIndex] = statWriter.defineField("INFO", String.format("varCountGroup[%s]", groups[groupIndex]),
                    1, ColumnType.Integer, String.format("Number of variant allele(s) called in group %s.", groups[groupIndex]));
        }
        depthFieldIndex = statWriter.defineField("INFO", "DP",
                1, ColumnType.Integer, "Total depth of sequencing across groups at this site");

        statWriter.defineSamples(samples);
        genotypeFormatter.defineGenotypeField(statWriter);
        statWriter.writeHeader();
    }


    @Override
    public void allocateStorage
            (
                    final int numberOfSamples,
                    final int numberOfGroups) {
        this.numberOfGroups = numberOfGroups;
        this.numberOfSamples = numberOfSamples;
        refCountsPerGroup = new int[numberOfGroups];
        variantsCountPerGroup = new int[numberOfGroups];
        distinctReadIndexCountPerGroup = new int[numberOfGroups];
        averageVariantQualityScorePerGroup = new float[numberOfGroups];
        alleleCountsGroupIndex = new int[numberOfGroups];
        alleleCountsMutableStrings = new MutableString[numberOfGroups];
        for (int i = 0; i < numberOfGroups; i++) {
            alleleCountsMutableStrings[i] = new MutableString();
        }

        // used for allelic fisher exact test:
        alleleCountsPerGroupBaseAllelesOnly = new int[numberOfAlleles][numberOfGroups];
        vectorBaseAllelesOnly = new int[numberOfAlleles * numberOfGroups];

        refCountsPerSample = new int[numberOfSamples];
        variantsCountPerSample = new int[numberOfSamples];
        distinctReadIndicesCountPerGroup = new IntSet[numberOfGroups];
        for (int i = 0; i < numberOfGroups; i++) {
            distinctReadIndicesCountPerGroup[i] = new IntArraySet();
        }
        genotypeFormatter = new GenotypesOutputFormat();
        genotypeFormatter.allocateStorage(numberOfSamples, numberOfGroups);
    }

    public void writeRecord(final DiscoverVariantIterateSortedAlignments iterator,
                            final SampleCountInfo[] sampleCounts,
                            final int referenceIndex,
                            int position,
                            final ObjectArrayList<PositionBaseInfo> list,
                            final int groupIndexA,
                            final int groupIndexB) {

        // report 1-based positions
        position = position + 1;
        int totalCount = 0;
        int maxGenotypeIndexAcrossSamples = 0;

        for (int i = 0; i < numberOfGroups; i++) {
            alleleCountsMutableStrings[i].setLength(0);
        }
        for (int sampleIndex = 0; sampleIndex < numberOfSamples; sampleIndex++) {
            final SampleCountInfo sci = sampleCounts[sampleIndex];
            int sumInSample = 0;

            //for (final int baseCount : sci.counts) {
            final int maxGenotypeIndex = sci.getGenotypeMaxIndex();
            maxGenotypeIndexAcrossSamples = Math.max(maxGenotypeIndex, maxGenotypeIndexAcrossSamples);
            for (int genotypeIndex = 0; genotypeIndex < maxGenotypeIndex; genotypeIndex++) {
                final int genotypeCount = sci.getGenotypeCount(genotypeIndex);
                totalCount += genotypeCount;
                sumInSample += genotypeCount;
                assert genotypeCount >= 0 : "counts must not be negative.";


            }
            // must observe at least one base in each sample to write output for this position.
            if (sumInSample == 0) {
                return;
            }
        }
        if (totalCount == 0) {
            return;
        }

        statWriter.setInfo(depthFieldIndex, totalCount);

        fillVariantCountArrays(sampleCounts, maxGenotypeIndexAcrossSamples);
        // since genotypes have been aligned across all samples, we use the first one to find the genotype string
        // corresponding to each index:
        final int sampleIndex = 0;
        final SampleCountInfo sci = sampleCounts[sampleIndex];

        for (int groupIndex = 0; groupIndex < numberOfGroups; groupIndex++) {
            for (int genotypeIndex = 0; genotypeIndex < maxGenotypeIndexAcrossSamples; genotypeIndex++) {
                // build the alleleCount string for each group:
                alleleCountsMutableStrings[groupIndex].append(sci.getGenotypeString(genotypeIndex));
                alleleCountsMutableStrings[groupIndex].append('=');
                alleleCountsMutableStrings[groupIndex].append(alleleCountsPerGroup[genotypeIndex][groupIndex]);
                alleleCountsMutableStrings[groupIndex].append(',');
            }
        }
        final CharSequence currentReferenceId = iterator.getReferenceId(referenceIndex);

        statWriter.setChromosome(currentReferenceId);
        statWriter.setPosition(position);
        statWriter.setReferenceAllele(Character.toString(sampleCounts[0].referenceBase));

        // construct a biomart region span in the format chr:pos1:chr:pos
        final String biomartRegionSpan = String.format("%s:%s:%s", currentReferenceId, position,
                position);

        statWriter.setInfo(biomartFieldIndex, biomartRegionSpan);

        for (int groupIndex = 0; groupIndex < numberOfGroups; groupIndex++) {

            // clip the trailing ','
            alleleCountsMutableStrings[groupIndex].setLength(alleleCountsMutableStrings[groupIndex].length() - 1);
            // set the value to the VCFWriter:
            statWriter.setInfo(alleleCountsGroupIndex[groupIndex], alleleCountsMutableStrings[groupIndex]);

            statWriter.setInfo(refCountsIndex[groupIndex], refCountsPerGroup[groupIndex]);
            statWriter.setInfo(varCountsIndex[groupIndex], variantsCountPerGroup[groupIndex]);
        }
        for (final GroupComparison comparison : groupComparisons) {
            final int indexGroup1 = comparison.indexGroup1;
            final int indexGroup2 = comparison.indexGroup1;
            final double denominator = (double) (refCountsPerGroup[indexGroup1] + 1) * (double) (variantsCountPerGroup[indexGroup2] + 1);
            double oddsRatio = denominator == 0 ? Double.NaN :
                    ((double) (refCountsPerGroup[indexGroup2] + 1) * (double) (variantsCountPerGroup[indexGroup1] + 1)) /
                            denominator;
            final double logOddsRatioSE;

            if (variantsCountPerGroup[indexGroup1] < 10 ||
                    variantsCountPerGroup[indexGroup2] < 10 ||
                    refCountsPerGroup[indexGroup1] < 10 ||
                    refCountsPerGroup[indexGroup2] < 10)

            {
                // standard error estimation is unreliable when any of the counts are less than 10.
                logOddsRatioSE = Double.NaN;
            } else

            {
                logOddsRatioSE = Math.sqrt(1d / refCountsPerGroup[indexGroup2] +
                        1d / variantsCountPerGroup[indexGroup1] +
                        1d / variantsCountPerGroup[indexGroup2] +
                        1d / refCountsPerGroup[indexGroup1]);
            }

            double log2OddsRatio = Math.log(oddsRatio) / Math.log(2);
            final double log2OddsRatioZValue = log2OddsRatio / logOddsRatioSE;

            double fisherP = Double.NaN;


            if (fisherRInstalled) {

                if (checkCounts()) {
                    final FisherExact.Result result = FisherExact.fexact(fisherVector, maxGenotypeIndexAcrossSamples,
                            numberOfGroups,
                            FisherExact.AlternativeHypothesis.twosided, true);
                    fisherP = result.getPValue();

                }
            }

            statWriter.setInfo(log2OddsRatioColumnIndex[comparison.index], log2OddsRatio);
            statWriter.setInfo(log2OddsRatioStandardErrorColumnIndex[comparison.index], logOddsRatioSE);
            statWriter.setInfo(log2OddsRatioZColumnIndex[comparison.index], log2OddsRatioZValue);
            statWriter.setInfo(fisherExactPValueColumnIndex[comparison.index], fisherP);

        }
        genotypeFormatter.writeGenotypes(statWriter, sampleCounts, position);

        statWriter.writeRecord();
    }


    public void close() {
        statWriter.close();
    }

    @Override
    public void setGenome(RandomAccessSequenceInterface genome) {

    }


    private void fillVariantCountArrays(final SampleCountInfo[] sampleCounts, int maxGenotypeIndexAcrossSamples) {

        Arrays.fill(variantsCountPerGroup, 0);
        Arrays.fill(refCountsPerGroup, 0);
        Arrays.fill(distinctReadIndexCountPerGroup, 0);

        if (maxGenotypeIndexAcrossSamples == SampleCountInfo.BASE_MAX_INDEX) {
            // reuse the vector if we have only SNPs at that site:

            fisherVector = vectorBaseAllelesOnly;
            Arrays.fill(fisherVector, 0);
            alleleCountsPerGroup = alleleCountsPerGroupBaseAllelesOnly;
            for (int alleleIndex = 0; alleleIndex < numberOfAlleles; alleleIndex++) {
                Arrays.fill(alleleCountsPerGroup[alleleIndex], 0);
            }
        } else {
            // allocate a new custom-size array if the site has indels:

            fisherVector = new int[maxGenotypeIndexAcrossSamples * numberOfGroups];
            alleleCountsPerGroup = new int[maxGenotypeIndexAcrossSamples][numberOfGroups];

            for (int genotypeIndex = 0; genotypeIndex < maxGenotypeIndexAcrossSamples; genotypeIndex++) {
                alleleCountsPerGroup[genotypeIndex] = new int[numberOfGroups];
            }
        }


        for (int groupIndex = 0; groupIndex < numberOfGroups; groupIndex++) {
            distinctReadIndicesCountPerGroup[groupIndex].clear();

        }
        for (final SampleCountInfo csi : sampleCounts) {
            final int sampleIndex = csi.sampleIndex;

            variantsCountPerSample[sampleIndex] = csi.varCount;
            refCountsPerSample[sampleIndex] = csi.refCount;

            final int groupIndex = readerIndexToGroupIndex[sampleIndex];
            variantsCountPerGroup[groupIndex] += csi.varCount;
            refCountsPerGroup[groupIndex] += csi.refCount;
            distinctReadIndicesCountPerGroup[groupIndex].addAll(csi.distinctReadIndices);
            for (int genotypeIndex = 0; genotypeIndex < maxGenotypeIndexAcrossSamples; genotypeIndex++) {

                final int alleleCount = csi.getGenotypeCount(genotypeIndex) > 0 ? 1 : 0;
                alleleCountsPerGroup[genotypeIndex][groupIndex] += alleleCount;
            }
        }

        // reorder allelic counts for fisher exact test:
        int j = 0;
        for (int groupIndex = 0; groupIndex < numberOfGroups; ++groupIndex) {
            for (int genotypeIndex = 0; genotypeIndex < maxGenotypeIndexAcrossSamples; genotypeIndex++) {

                fisherVector[j] = alleleCountsPerGroup[genotypeIndex][groupIndex];
                ++j;
            }
        }

        //fisher.test(matrix(c(10,10,10,0,0,10,10,10,0,0),5,2), hybrid="TRUE")
        for (int groupIndex = 0; groupIndex < numberOfGroups; groupIndex++) {
            distinctReadIndexCountPerGroup[groupIndex] = distinctReadIndicesCountPerGroup[groupIndex].size();
        }
    }

    private boolean checkCounts() {
        boolean ok = true;
        int totalCount = 0;
        // detect if any count is negative
        for (final int count : refCountsPerGroup) {

            if (count < 0) ok = false;
            totalCount += count;
        }
        for (final int count : variantsCountPerGroup) {
            if (count < 0) ok = false;
            totalCount += count;
        }
        // skip fisher if all counts are zero:
        if (totalCount == 0) return false;
        return ok;
    }


}
