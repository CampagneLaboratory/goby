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

import edu.cornell.med.icb.goby.modes.SequenceVariationOutputFormat;
import edu.cornell.med.icb.goby.modes.DiscoverSequenceVariantsMode;
import edu.cornell.med.icb.goby.stats.StatisticsWriter;
import edu.cornell.med.icb.goby.readers.vcf.ColumnType;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.ints.IntArraySet;

import java.util.Arrays;

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
    private StatisticsWriter statWriter;
    String[] samples;
    private boolean outputVCF;

    public void defineColumns(StatisticsWriter statsWriter, DiscoverSequenceVariantsMode mode) {
        samples = mode.getSamples();
        refIdColumnIndex = statsWriter.defineColumn("chr:position:position");
        statsWriter.defineColumnAttributes("ID", 1, ColumnType.String, "Unique site id, in the format chr:start:end.", "chr:position:position");

        statsWriter.defineColumnSet(samples,
                "refProportion[%s]"
        );
        statsWriter.defineColumnSet(samples,
                "count[%s]"
        );
        statsWriter.defineColumnAttributes(1, ColumnType.Float, samples, "refProportion[%s]", "Zygosity[%s]");
        statsWriter.defineColumnAttributes(1, ColumnType.Integer, samples, "count[%s]");

        this.statWriter = statsWriter;
        statsWriter.setOutputVCF(outputVCF);
        statsWriter.writeHeader();
    }

    public void allocateStorage(int numberOfSamples, int numberOfGroups) {
        this.numberOfGroups = numberOfGroups;
        this.numberOfSamples = numberOfSamples;

        refCountsPerSample = new int[numberOfSamples];
        variantsCountPerSample = new int[numberOfSamples];
    }

    public void writeRecord(DiscoverVariantIterateSortedAlignments iterator, SampleCountInfo[] sampleCounts,
                            int referenceIndex, int position, ObjectArrayList<PositionBaseInfo> list, int groupIndexA, int groupIndexB) {
        fillVariantCountArrays(sampleCounts);


        CharSequence currentReferenceId = iterator.getReferenceId(referenceIndex);

        final int positionString = position + 1;
        statWriter.setValue(refIdColumnIndex, String.format("%s:%d:%d", currentReferenceId, positionString, positionString));

        for (int sampleIndex = 0; sampleIndex < numberOfSamples; sampleIndex++) {
            int numAlleles = 0;
            int totalCount = 0;
            for (int count : sampleCounts[sampleIndex].counts) {
                if (count > 0) numAlleles++;
                totalCount += count;
            }
            if (numAlleles == 2) {
                // need exactly two alleles to estimate reference allele imbalance:
                double refProportion = (double) refCountsPerSample[sampleIndex];
                refProportion /= refCountsPerSample[sampleIndex] + variantsCountPerSample[sampleIndex];
                statWriter.setValue(refProportion,
                        "refProportion[%s]", samples[sampleIndex]);
            } else {
                statWriter.setValue("",
                        "refProportion[%s]", samples[sampleIndex]);
            }
            statWriter.setValue(totalCount,
                    "count[%s]", samples[sampleIndex]);
        }
        statWriter.writeRecord();
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
