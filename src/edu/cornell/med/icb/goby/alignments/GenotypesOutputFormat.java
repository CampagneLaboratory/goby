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
    public int baseCountFieldIndex;
    private int zygFieldIndex;
    private int failBaseCountFieldIndex;
    private int goodBaseCountFieldIndex;

    public void defineColumns(PrintWriter writer, DiscoverSequenceVariantsMode mode) {
        samples = mode.getSamples();
        this.statsWriter = new VCFWriter(writer);


        biomartFieldIndex = statsWriter.defineField("INFO", "BIOMART_COORDS", 1, ColumnType.String, "Coordinates for use with Biomart.");

        defineGenotypeField(statsWriter);
        zygFieldIndex = statsWriter.defineField("FORMAT", "Zygosity", 1, ColumnType.String, "Zygosity");
        statsWriter.defineSamples(samples);
        statsWriter.writeHeader();
    }

    public void defineGenotypeField(VCFWriter statsWriter) {
        genotypeFieldIndex = statsWriter.defineField("FORMAT", "GT", 1, ColumnType.String, "Genotype");
        baseCountFieldIndex = statsWriter.defineField("FORMAT", "BC", 1, ColumnType.String, "Base counts in format A=?;T=?;C=?;G=?;N=?.");
        goodBaseCountFieldIndex = statsWriter.defineField("FORMAT", "GB", 1, ColumnType.String, "Number of bases that pass base filters in this sample.");
        failBaseCountFieldIndex = statsWriter.defineField("FORMAT", "FB", 1, ColumnType.String, "Number of bases that failed base filters in this sample.");
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

        position = position - 1;
        fillVariantCountArrays(sampleCounts);

        CharSequence currentReferenceId = iterator.getReferenceId(referenceIndex);

        statsWriter.setId(".");
        statsWriter.setInfo(biomartFieldIndex,
                String.format("%s:%d:%d", currentReferenceId, position,
                        position));
        statsWriter.setChromosome(currentReferenceId);
        statsWriter.setPosition(position);

        writeGenotypes(statsWriter, sampleCounts);

        writeZygozity(sampleCounts);

        statsWriter.writeRecord();
    }

    private void writeZygozity(SampleCountInfo[] sampleCounts) {
        for (int sampleIndex = 0; sampleIndex < numberOfSamples; sampleIndex++) {
            SampleCountInfo sci = sampleCounts[sampleIndex];

            String zygozity;
            switch (alleleSet.size()) {
                case 0:
                    zygozity = "not-typed";
                    break;
                case 1:
                    zygozity = "homozygous";
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

        }
    }

    public void writeGenotypes(VCFWriter statsWriter, SampleCountInfo[] sampleCounts) {

        for (int sampleIndex = 0; sampleIndex < numberOfSamples; sampleIndex++) {

            alleleSet.clear();
            SampleCountInfo sci = sampleCounts[sampleIndex];
            int totalCount = 0;
            for (int sampleCount : sci.counts) {
                totalCount += sampleCount;
            }
            //  System.out.printf("totalCount %d failedCount %d%n",totalCount,sci.failedCount);
            statsWriter.setSampleValue(goodBaseCountFieldIndex, sampleIndex, totalCount);
            statsWriter.setSampleValue(failBaseCountFieldIndex, sampleIndex, sci.failedCount);

            int baseIndex = 0;
            int genotypeCount = 0;
            genotypeBuffer.setLength(0);
            char refBase = '.';
            String referenceAllele = ".";
            boolean siteObserved = false;

            for (int count : sci.counts) {
                final char base = sci.base(baseIndex);
                if (count > 0) {
                    siteObserved = true;
                    alleleSet.add(base);
                    if (base != sci.referenceBase) {
                        statsWriter.addAlternateAllele(Character.toString(base));

                        genotypeBuffer.append(String.format("%c/", base));
                        genotypeCount++;
                    } else {
                        referenceAllele = Character.toString(sci.referenceBase);
                        genotypeBuffer.append(String.format("%c/", base));
                    }
                }
                baseIndex++;
            }

            if (alleleSet.size() == 1) {
                // when homozygous genotype 0/ write 0/0/ (last / will be removed when length is adjusted)
                genotypeBuffer.append(genotypeBuffer.copy());
            }
            MutableString baseCountString = new MutableString();
            baseIndex = 0;
            for (int count : sci.counts) {
                final char base = sci.base(baseIndex);
                baseCountString.append(base);
                baseCountString.append('=');
                baseCountString.append(Integer.toString(count));
                baseIndex++;
                baseCountString.append(',');
            }
            baseCountString.setLength(baseCountString.length() - 1);
            statsWriter.setSampleValue(baseCountFieldIndex, sampleIndex, baseCountString);

            if (siteObserved) {

                if (genotypeBuffer.length() > 1) {
                    // trim the trailing coma:
                    genotypeBuffer.setLength(genotypeBuffer.length() - 1);
                }
                statsWriter.setReferenceAllele(referenceAllele);
                statsWriter.setSampleValue(genotypeFieldIndex, sampleIndex, statsWriter.codeGenotype(genotypeBuffer.toString()));

            } else {
                statsWriter.setSampleValue(genotypeFieldIndex, sampleIndex, "./.");

            }
        }
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