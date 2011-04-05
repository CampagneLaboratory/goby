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

package edu.cornell.med.icb.goby.alignments;

import edu.cornell.med.icb.goby.modes.DiscoverSequenceVariantsMode;
import edu.cornell.med.icb.goby.modes.SequenceVariationOutputFormat;
import edu.cornell.med.icb.goby.readers.vcf.ColumnType;
import edu.cornell.med.icb.goby.stats.VCFWriter;
import edu.cornell.med.icb.goby.stats.TTestCalculator;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;

import java.io.PrintWriter;

import org.apache.commons.math.stat.inference.TTest;
import org.apache.commons.math.stat.inference.TTestImpl;
import org.apache.commons.math.MathException;

/**
 * @author Fabien Campagne
 *         Date: Mar 21, 2011
 *         Time: 2:37:43 PM
 */
public class AlleleFrequencyOutputFormat implements SequenceVariationOutputFormat {
    private int refIdColumnIndex;
    private int positionColumnIndex;
    private int numberOfGroups;
    private int numberOfSamples;
    private int[] refCountsPerSample;
    private int[] variantsCountPerSample;
    private VCFWriter statsWriter;
    String[] samples;
    private boolean outputVCF;
    private int refPropFieldIndex;
    private int biomartFieldIndex;
    GenotypesOutputFormat genotypeFormatter = new GenotypesOutputFormat();
    private int depthFieldIndex;

    private double[] valuesGroupsA;
    private double[] valuesGroupsB;
    private int pValueIndex;

    public void defineColumns(PrintWriter writer, DiscoverSequenceVariantsMode mode) {
        samples = mode.getSamples();
        statsWriter = new VCFWriter(writer);
        biomartFieldIndex = statsWriter.defineField("INFO", "BIOMART_COORDS", 1, ColumnType.String, "Coordinates for use with Biomart.");
        pValueIndex = statsWriter.defineField("INFO", "P", 1, ColumnType.Float, "P-values of a t-test comparing arcsin transformed values of the ratio between groups.");


        refPropFieldIndex = statsWriter.defineField("FORMAT", "RP", 1, ColumnType.Float,
                "Proportion of reference allele in the sample (count(ref)/(sum count other alleles)");
        depthFieldIndex = statsWriter.defineField("INFO", "DP",
                1, ColumnType.Integer, "Total depth of sequencing across groups at this site");
        genotypeFormatter.defineGenotypeField(statsWriter);

        readerIndexToGroupIndex = mode.getReaderIndexToGroupIndex();
        int countA = 0;
        int countB = 0;
        for (int groupIndex : readerIndexToGroupIndex) {
            countA += groupIndex == 1 ? 0 : 1;
            countB += groupIndex == 1 ? 1 : 0;
        }
        valuesGroupsA = new double[countA];
        valuesGroupsB = new double[countB];
        statsWriter.defineSamples(samples);
        statsWriter.writeHeader();

    }

    private final TTest mathCommonsTTest = new TTestImpl();

    public void allocateStorage(int numberOfSamples, int numberOfGroups) {
        this.numberOfGroups = numberOfGroups;
        this.numberOfSamples = numberOfSamples;

        refCountsPerSample = new int[numberOfSamples];
        variantsCountPerSample = new int[numberOfSamples];
        genotypeFormatter = new GenotypesOutputFormat();
        genotypeFormatter.allocateStorage(numberOfSamples, numberOfGroups);


    }

    public void writeRecord(DiscoverVariantIterateSortedAlignments iterator, SampleCountInfo[] sampleCounts,
                            int referenceIndex, int position, ObjectArrayList<PositionBaseInfo> list, int groupIndexA, int groupIndexB) {
        fillVariantCountArrays(sampleCounts);
        int totalCount = 0;
        for (int sampleIndex = 0; sampleIndex < numberOfSamples; sampleIndex++) {
            SampleCountInfo sci = sampleCounts[sampleIndex];
            int sumInSample = 0;
            for (int baseCount : sci.counts) {
                totalCount += baseCount;
                sumInSample += baseCount;
                assert baseCount >= 0 : "counts must not be negative.";
            }
            // must observe at least one base in each sample to write output for this position.
            if (sumInSample == 0) return;
        }
        if (totalCount == 0) return;
        statsWriter.setInfo(depthFieldIndex, totalCount);

        CharSequence currentReferenceId = iterator.getReferenceId(referenceIndex);

        statsWriter.setId(".");
        statsWriter.setInfo(biomartFieldIndex,
                String.format("%s:%d:%d", currentReferenceId, position,
                        position));
        statsWriter.setChromosome(currentReferenceId);
        statsWriter.setPosition(position);

        int valuesGroupAIndex = 0;
        int valuesGroupBIndex = 0;

        for (int sampleIndex = 0; sampleIndex < numberOfSamples; sampleIndex++) {
            int numAlleles = 0;
            totalCount = 0;
            for (int count : sampleCounts[sampleIndex].counts) {
                if (count > 0) numAlleles++;
                totalCount += count;
            }

            // estimate reference allele proportion:
            double refProportion = (double) refCountsPerSample[sampleIndex];
            refProportion /= (refCountsPerSample[sampleIndex] + variantsCountPerSample[sampleIndex]);
            statsWriter.setSampleValue(refPropFieldIndex, sampleIndex, refProportion);
            int groupIndex = readerIndexToGroupIndex[sampleIndex];
            final double transformedValue = StrictMath.asin(StrictMath.sqrt(refProportion));

            if (groupIndex == 0) {
                valuesGroupsA[valuesGroupAIndex++] = transformedValue;
            } else {
                valuesGroupsB[valuesGroupBIndex++] = transformedValue;
            }
        }
        assert valuesGroupAIndex + valuesGroupBIndex == numberOfSamples;
        genotypeFormatter.writeGenotypes(statsWriter, sampleCounts);
        double pValue = 1;
        try {
            pValue = mathCommonsTTest.homoscedasticTTest(valuesGroupsA, valuesGroupsB);

        } catch (MathException e) {
            pValue = 1;
        }
        statsWriter.setInfo(pValueIndex, pValue);
        statsWriter.writeRecord();
    }

    public void close() {
        statsWriter.close();
    }

    public void outputVCF(boolean state) {
        outputVCF = state;
    }

    int[] readerIndexToGroupIndex;

    private void fillVariantCountArrays(SampleCountInfo[] sampleCounts) {


        for (SampleCountInfo csi : sampleCounts) {
            final int sampleIndex = csi.sampleIndex;
            variantsCountPerSample[sampleIndex] = csi.varCount;
            refCountsPerSample[sampleIndex] = csi.refCount;
            int groupIndex = readerIndexToGroupIndex[sampleIndex];
        }

        for (SampleCountInfo csi : sampleCounts) {


        }
    }

}
