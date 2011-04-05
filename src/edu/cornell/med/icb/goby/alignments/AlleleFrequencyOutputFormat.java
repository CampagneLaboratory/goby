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
import it.unimi.dsi.fastutil.objects.ObjectArrayList;

import java.io.PrintWriter;

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

    public void defineColumns(PrintWriter writer, DiscoverSequenceVariantsMode mode) {
        samples = mode.getSamples();
        statsWriter = new VCFWriter(writer);
        biomartFieldIndex = statsWriter.defineField("INFO", "BIOMART_COORDS", 1, ColumnType.String, "Coordinates for use with Biomart.");


        refPropFieldIndex = statsWriter.defineField("FORMAT", "RP", 1, ColumnType.Float,
                "Proportion of reference allele in the sample (count(ref)/(sum count other alleles)");
        depthFieldIndex = statsWriter.defineField("INFO", "DP",
                1, ColumnType.Integer, "Total depth of sequencing across groups at this site");
        genotypeFormatter.defineGenotypeField(statsWriter);

        statsWriter.defineSamples(samples);
        statsWriter.writeHeader();

    }

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


        for (int sampleIndex = 0; sampleIndex < numberOfSamples; sampleIndex++) {
            int numAlleles = 0;
            totalCount = 0;
            for (int count : sampleCounts[sampleIndex].counts) {
                if (count > 0) numAlleles++;
                totalCount += count;
            }
            if (numAlleles == 2) {
                // need exactly two alleles to estimate reference allele imbalance:
                double refProportion = (double) refCountsPerSample[sampleIndex];
                refProportion /= refCountsPerSample[sampleIndex] + variantsCountPerSample[sampleIndex];
                statsWriter.setSampleValue(refPropFieldIndex, sampleIndex, refProportion);
            }
        }
        genotypeFormatter.writeGenotypes(statsWriter, sampleCounts);

        statsWriter.writeRecord();
    }

    public void close() {
        statsWriter.close();
    }

    public void outputVCF(boolean state) {
        outputVCF = state;
    }

    private void fillVariantCountArrays(SampleCountInfo[] sampleCounts) {


        for (SampleCountInfo csi : sampleCounts) {
            final int sampleIndex = csi.sampleIndex;
            variantsCountPerSample[sampleIndex] = csi.varCount;
            refCountsPerSample[sampleIndex] = csi.refCount;
        }

    }

}
