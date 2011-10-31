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
import edu.cornell.med.icb.goby.alignments.*;
import edu.cornell.med.icb.goby.readers.vcf.ColumnType;
import edu.cornell.med.icb.goby.stats.DifferentialExpressionAnalysis;
import edu.cornell.med.icb.goby.stats.DifferentialExpressionCalculator;
import edu.cornell.med.icb.goby.stats.FisherExactRCalculator;
import edu.cornell.med.icb.goby.stats.VCFWriter;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.rosuda.JRI.Rengine;

import java.io.PrintWriter;
import java.util.Arrays;

/**
 * A Variant Call Format output to estimate methylation rates for a set of samples and find methylation rate differences
 * between group.
 * TODO evaluate if supporting indel genotypes would make any difference.
 *
 * @author Fabien Campagne
 *         Date: April 4 2011
 *         Time: 2:38:13 AM
 */
public class MethylationRateVCFOutputFormat implements SequenceVariationOutputFormat {

    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(MethylationRateVCFOutputFormat.class);

    VCFWriter statWriter;
    private int refIdColumnIndex;
    private int positionColumnIndex;

    private boolean fisherRInstalled;
    private String[] groups;
    private String[] samples;
    private ObjectArrayList<ReadIndexStats> readIndexStats;
    private int log2OddsRatioColumnIndex;
    private int fisherExactPValueColumnIndex;
    private int numberOfGroups;

    private DifferentialExpressionAnalysis deAnalyzer;
    private DifferentialExpressionCalculator deCalculator;

    private int log2OddsRatioStandardErrorColumnIndex;
    private int log2OddsRatioZColumnIndex;
    int[] readerIndexToGroupIndex;
    private int[] unmethylatedCCountsPerGroup;
    private int[] methylatedCCountPerGroup;

    private IntSet[] distinctReadIndicesCountPerGroup;

    private int[] unmethylatedCCountPerSample;
    private int[] methylatedCCountsPerSample;

    private int numberOfSamples;
    private int biomartFieldIndex;
    private GenotypesOutputFormat genotypeFormatter;
    private int depthFieldIndex;
    private int[] methylatedCCountsIndex;
    private int[] notMethylatedCCountsIndex;
    private int eventCountAtSite;
    private int methylationRateFieldIndex;
    private char strandAtSite;
    private int strandFieldIndex;

    public void setMinimumEventThreshold(final int minimumEventThreshold) {
        this.minimumEventThreshold = minimumEventThreshold;
    }

    private int minimumEventThreshold = 1;


    public void defineColumns(final PrintWriter writer, final DiscoverSequenceVariantsMode mode) {
        deAnalyzer = mode.getDiffExpAnalyzer();
        deCalculator = mode.getDiffExpCalculator();
        groups = mode.getGroups();
        samples = mode.getSamples();
        readerIndexToGroupIndex = mode.getReaderIndexToGroupIndex();
        final ObjectArrayList<ReadIndexStats> readIndexStats = mode.getReadIndexStats();
        this.statWriter = new VCFWriter(writer);
        try {
            //activate R only if we need it:
            final Rengine rEngine = GobyRengine.getInstance().getRengine();
            fisherRInstalled = rEngine != null && rEngine.isAlive();
        } catch (java.lang.UnsatisfiedLinkError e) {
            System.out.println("Cannot initialize R");
            e.printStackTrace();
            throw e;
        }
        if (groups.length != 2) {
            System.err.println("CompareGroupsVCFOutputFormat requires exactly two groups.");
            System.exit(1);
        }

        this.readIndexStats = readIndexStats;


        log2OddsRatioColumnIndex = -1;
        fisherExactPValueColumnIndex = -1;
        numberOfGroups = groups.length;
        biomartFieldIndex = statWriter.defineField("INFO", "BIOMART_COORDS", 1, ColumnType.String, "Coordinates for use with Biomart.");
        strandFieldIndex = statWriter.defineField("INFO", "Strand", 1, ColumnType.String, "Strand of the cytosine site on the reference sequence.");


        log2OddsRatioColumnIndex = statWriter.defineField("INFO", String.format("LOD[%s/%s]", groups[0], groups[1]),
                1, ColumnType.Float, String.format("Log2 of the odds-ratio of observing methylation in  group %s versus group %s", groups[0], groups[1]));

        log2OddsRatioStandardErrorColumnIndex = statWriter.defineField("INFO", String.format("LOD_SE[%s/%s]", groups[0], groups[1]),
                1, ColumnType.Float, String.format("Standard Error of the log2 of the odds-ratio between group %s and group %s", groups[0], groups[1]));

        log2OddsRatioZColumnIndex = statWriter.defineField("INFO", String.format("LOD_Z[%s/%s]", groups[0], groups[1]),
                1, ColumnType.Float, String.format("Z value of the odds-ratio between group %s and group %s", groups[0], groups[1]));

        fisherExactPValueColumnIndex = statWriter.defineField("INFO", String.format("FisherP[%s/%s]", groups[0], groups[1]),
                1, ColumnType.Float, String.format("Fisher exact P-value of observing as large a difference by chance between group %s and group %s.", groups[0], groups[1]));

        methylatedCCountsIndex = new int[groups.length];
        notMethylatedCCountsIndex = new int[groups.length];

        for (int groupIndex = 0; groupIndex < groups.length; groupIndex++) {
            notMethylatedCCountsIndex[groupIndex] = statWriter.defineField("INFO", String.format("#C Group[%s]", groups[groupIndex]),
                    1, ColumnType.Integer, String.format("Number of unmethylated Cytosines at site in group %s.", groups[groupIndex]));
            methylatedCCountsIndex[groupIndex] = statWriter.defineField("INFO", String.format("#Cm Group[%s]", groups[groupIndex]),
                    1, ColumnType.Integer, String.format("Number of methylated Cytosines at site in group %s.", groups[groupIndex]));
        }
        depthFieldIndex = statWriter.defineField("INFO", "DP",
                1, ColumnType.Integer, "Total depth of sequencing across groups at this site");

        statWriter.defineSamples(samples);
        // define Genotype as first field (required by VCF specification):
        genotypeFormatter.defineGenotypeField(statWriter);
        methylationRateFieldIndex = statWriter.defineField("FORMAT", "MR", 1, ColumnType.Integer, "Methylation rate. 0-100%, 100% indicate fully methylated.");

        statWriter.writeHeader();
    }

    @Override
    public void allocateStorage(final int numberOfSamples, final int numberOfGroups) {
        this.numberOfGroups = numberOfGroups;
        this.numberOfSamples = numberOfSamples;
        unmethylatedCCountsPerGroup = new int[numberOfGroups];
        methylatedCCountPerGroup = new int[numberOfGroups];


        methylatedCCountsPerSample = new int[numberOfSamples];
        unmethylatedCCountPerSample = new int[numberOfSamples];
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
                            final ObjectArrayList<PositionBaseInfo> list,
                            final int groupIndexA,
                            final int groupIndexB) {

        position = position + 1;
        final char refBase = sampleCounts[0].referenceBase;
        if (refBase != 'C' && refBase != 'G') {
            return;
        }
        fillMethylationCountArrays(sampleCounts, list, position, refBase);
        if (eventCountAtSite < minimumEventThreshold) {
            return;
        }
        statWriter.setInfo(depthFieldIndex, list.size());
        final CharSequence currentReferenceId = iterator.getReferenceId(referenceIndex);

        statWriter.setChromosome(currentReferenceId);
        statWriter.setPosition(position);
        statWriter.setReferenceAllele(Character.toString(sampleCounts[0].referenceBase));

        // construct a biomart region span in the format chr:pos1:chr:pos
        final String biomartRegionSpan = String.format("%s:%s:%s", currentReferenceId, position,
                position);

        statWriter.setInfo(biomartFieldIndex, biomartRegionSpan);
        statWriter.setInfo(strandFieldIndex, Character.toString(strandAtSite));

        for (int groupIndex = 0; groupIndex < numberOfGroups; groupIndex++) {
            statWriter.setInfo(notMethylatedCCountsIndex[groupIndex], unmethylatedCCountsPerGroup[groupIndex]);
            statWriter.setInfo(methylatedCCountsIndex[groupIndex], methylatedCCountPerGroup[groupIndex]);
        }

        for (int sampleIndex = 0; sampleIndex < numberOfSamples; sampleIndex++) {
            final float numerator = methylatedCCountsPerSample[sampleIndex];
            final float denominator = numerator + unmethylatedCCountPerSample[sampleIndex];

            final float methylationRate = numerator * 100 / denominator;
            statWriter.setSampleValue(methylationRateFieldIndex, sampleIndex, Math.round(methylationRate));
        }

        final double denominator = (double) (unmethylatedCCountsPerGroup[groupIndexA]) * (double) (methylatedCCountPerGroup[groupIndexB]);
        final double oddsRatio = denominator == 0 ? Double.NaN :
                ((double) (unmethylatedCCountsPerGroup[groupIndexB]) * (double) (methylatedCCountPerGroup[groupIndexA])) /
                        denominator;
        final double logOddsRatioSE;

        if (methylatedCCountPerGroup[groupIndexA] < 10 ||
                methylatedCCountPerGroup[groupIndexB] < 10 ||
                unmethylatedCCountsPerGroup[groupIndexA] < 10 ||
                unmethylatedCCountsPerGroup[groupIndexB] < 10) {
            // standard error estimation is unreliable when any of the counts are less than 10.
            logOddsRatioSE = Double.NaN;
        } else {
            logOddsRatioSE = Math.sqrt(1d / unmethylatedCCountsPerGroup[groupIndexB] +
                    1d / methylatedCCountPerGroup[groupIndexA] +
                    1d / methylatedCCountPerGroup[groupIndexB] +
                    1d / unmethylatedCCountsPerGroup[groupIndexA]);
        }
        final double log2OddsRatio = Math.log(oddsRatio) / Math.log(2);
        final double log2OddsRatioZValue = log2OddsRatio / logOddsRatioSE;
        double fisherP = Double.NaN;
        if (eventCountAtSite >= 10) {
            // estimate Fisher only if we have seen at least 10 events.

            final boolean ok = checkCounts();
            if (ok) {
                fisherP = fisherRInstalled ? FisherExactRCalculator.getFisherPValue(
                        unmethylatedCCountsPerGroup[groupIndexB], methylatedCCountPerGroup[groupIndexB],
                        unmethylatedCCountsPerGroup[groupIndexA], methylatedCCountPerGroup[groupIndexA]) : Double.NaN;
            } else {
                System.err.printf("An exception was caught evaluating the Fisher Exact test P-value. Details are provided below%n" +
                        "referenceId=%s referenceIndex=%d position=%d %n" +
                        "unmethylatedCCountsPerGroup[1]=%d methylatedCCountPerGroup[1]=%d%n" +
                        "unmethylatedCCountsPerGroup[0]=%d, methylatedCCountPerGroup[0]=%d",
                        currentReferenceId, referenceIndex,
                        position,
                        unmethylatedCCountsPerGroup[groupIndexB], methylatedCCountPerGroup[groupIndexB],
                        unmethylatedCCountsPerGroup[groupIndexA], methylatedCCountPerGroup[groupIndexA]
                );
            }
        }

        statWriter.setInfo(log2OddsRatioColumnIndex, log2OddsRatio);
        statWriter.setInfo(log2OddsRatioStandardErrorColumnIndex, logOddsRatioSE);
        statWriter.setInfo(log2OddsRatioZColumnIndex, log2OddsRatioZValue);
        statWriter.setInfo(fisherExactPValueColumnIndex, fisherP);
        genotypeFormatter.writeGenotypes(statWriter, sampleCounts, position);

        statWriter.writeRecord();
    }


    public void close() {
        statWriter.close();
    }


    private void fillMethylationCountArrays(final SampleCountInfo[] sampleCounts, final ObjectArrayList<PositionBaseInfo> list,
                                            final int position, final char refBase) {

        Arrays.fill(methylatedCCountPerGroup, 0);
        Arrays.fill(unmethylatedCCountsPerGroup, 0);
        Arrays.fill(methylatedCCountsPerSample, 0);
        Arrays.fill(unmethylatedCCountPerSample, 0);
        eventCountAtSite = 0;
        strandAtSite = '?';

        if (refBase == 'C') {
            strandAtSite = '+';
        } else if (refBase == 'G') {
            strandAtSite = '-';
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
                    ++methylatedCCountsPerSample[sampleIndex];
                    ++methylatedCCountPerGroup[groupIndex];
                    ++eventCountAtSite;
                } else {
                    // C became T on forward strand (G->A on reverse) indicates that the Cytosine was not methylated and was converted.
                    ++unmethylatedCCountPerSample[sampleIndex];
                    ++unmethylatedCCountsPerGroup[groupIndex];
                    ++eventCountAtSite;
                }

            }
        }

    }

    private boolean checkCounts() {
        boolean ok = true;
        // detect if any count is negative (that's a bug)
        for (final int count : unmethylatedCCountsPerGroup) {

            if (count < 0) {
                ok = false;
            }
        }
        for (final int count : methylatedCCountPerGroup) {
            if (count < 0) {
                ok = false;
            }
        }
        return ok;
    }


}