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
import it.unimi.dsi.fastutil.chars.CharArraySet;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.lang.MutableString;

import java.io.PrintWriter;

/**
 * @author Fabien Campagne
 *         Date: Mar 21, 2011
 *         Time: 2:37:43 PM
 */
public class GenotypesOutputFormat implements SequenceVariationOutputFormat {
    private int biomartCoordColumnIndex;
    private int positionColumnIndex;
    private int numberOfGroups;
    private int numberOfSamples;
    private int[] refCountsPerSample;
    private int[] variantsCountPerSample;
    private VCFWriter statsWriter;
    String[] samples;
    private int chromosomeColumnIndex;
    private int idColumnIndex;
    private int biomartFieldIndex;
    private int genotypeFieldIndex;
    private int zygFieldIndex;


    public void defineColumns(PrintWriter writer, DiscoverSequenceVariantsMode mode) {
        samples = mode.getSamples();
        this.statsWriter = new VCFWriter(writer);


        biomartFieldIndex = statsWriter.defineField("INFO", "BIOMART_COORDS", 1, ColumnType.String, "Coordinates for use with Biomart.");

        genotypeFieldIndex = statsWriter.defineField("FORMAT", "GT", 1, ColumnType.String, "Genotype");
        zygFieldIndex = statsWriter.defineField("FORMAT", "Zygosity", 1, ColumnType.String, "Zygosity");
        statsWriter.defineSamples(samples);
        statsWriter.writeHeader();
    }

    public void allocateStorage(int numberOfSamples, int numberOfGroups) {
        this.numberOfGroups = numberOfGroups;
        this.numberOfSamples = numberOfSamples;

        refCountsPerSample = new int[numberOfSamples];
        variantsCountPerSample = new int[numberOfSamples];
    }

    IntArrayList decreasingCounts = new IntArrayList();
    CharArraySet alleleSet = new CharArraySet();
    MutableString genotypeBuffer = new MutableString();

    public void writeRecord(DiscoverVariantIterateSortedAlignments iterator, SampleCountInfo[] sampleCounts,
                            int referenceIndex, int position, ObjectArrayList<PositionBaseInfo> list, int groupIndexA, int groupIndexB) {

        fillVariantCountArrays(sampleCounts);

        CharSequence currentReferenceId = iterator.getReferenceId(referenceIndex);

        statsWriter.setId(".");
        statsWriter.setInfo(biomartCoordColumnIndex,
                String.format("%s:%d:%d", currentReferenceId, position,
                        position));
        statsWriter.setChromosome(currentReferenceId);
        statsWriter.setPosition(position);
        for (int sampleIndex = 0; sampleIndex < numberOfSamples; sampleIndex++) {
            alleleSet.clear();
            SampleCountInfo sci = sampleCounts[sampleIndex];
            int totalCount = 0;
            for (int sampleCount : sci.counts) {
                totalCount += sampleCount;
            }

            int baseIndex = 0;
            int genotypeCount = 0;
            genotypeBuffer.setLength(0);
            char refBase = '.';
            for (int count : sci.counts) {
                final char base = sci.base(baseIndex);
                if (count > 0) {
                    alleleSet.add(base);
                    if (base != sci.referenceBase) {
                        genotypeBuffer.append(String.format("%c/%c,", sci.referenceBase, base));
                        genotypeCount++;
                    }
                }
                baseIndex++;
            }
            if (genotypeBuffer.length() > 1) {
                // trim the trailing coma:
                genotypeBuffer.setLength(genotypeBuffer.length() - 1);
            }


            String zygozity;
            switch (alleleSet.size()) {
                case 0:
                    zygozity = "not-typed";
                    break;
                case 1:
                    zygozity = "homozygous";
                    genotypeBuffer.setLength(0);
                    genotypeBuffer.append(String.format("%c/%c", sci.referenceBase, sci.referenceBase));
                    break;
                case 2:
                    zygozity = "heterozygous";
                    break;
                default:
                    /*     zygozity = String.format("mixture:%d:%d:%d:%d:%d",
                    sci.counts[SampleCountInfo.BASE_A_INDEX],
                    sci.counts[SampleCountInfo.BASE_T_INDEX],
                    sci.counts[SampleCountInfo.BASE_C_INDEX],
                    sci.counts[SampleCountInfo.BASE_G_INDEX],
                    sci.counts[SampleCountInfo.BASE_OTHER_INDEX])
                    ;*/
                    zygozity = "Mixture";
                    break;
            }

            statsWriter.setSampleValue(zygFieldIndex, sampleIndex, zygozity);
            
            statsWriter.setSampleValue(genotypeFieldIndex, sampleIndex, genotypeBuffer);

        }

        statsWriter.writeRecord();
    }

    public void close() {
        statsWriter.close();
    }


    private void fillVariantCountArrays(SampleCountInfo[] sampleCounts) {


        for (SampleCountInfo csi : sampleCounts) {
            final int sampleIndex = csi.sampleIndex;
            variantsCountPerSample[sampleIndex] = csi.varCount;
            refCountsPerSample[sampleIndex] = csi.refCount;
        }

    }

}