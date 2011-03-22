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
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.chars.CharArrayList;
import it.unimi.dsi.fastutil.chars.CharArraySet;
import it.unimi.dsi.lang.MutableString;

import java.util.Arrays;
import java.util.Collections;

/**
 * @author Fabien Campagne
 *         Date: Mar 21, 2011
 *         Time: 2:37:43 PM
 */
public class GenotypesOutputFormat implements SequenceVariationOutputFormat {
    private int refIdColumnIndex;
    private int positionColumnIndex;
    private int numberOfGroups;
    private int numberOfSamples;
    private int[] refCountsPerSample;
    private int[] variantsCountPerSample;
    private StatisticsWriter statWriter;
    String[] samples;

    public void defineColumns(StatisticsWriter statsWriter, DiscoverSequenceVariantsMode mode) {
        samples = mode.getSamples();
        refIdColumnIndex = statsWriter.defineColumn("chr:position:chr:position");

        statsWriter.defineColumnSet(samples,
                "Genotype[%s]"
        );
        statsWriter.defineColumnSet(samples,
                "Zygosity[%s]"
        );
        this.statWriter = statsWriter;
        statsWriter.writeHeader();
    }

    public void allocateStorage(int numberOfSamples, int numberOfGroups) {
        this.numberOfGroups = numberOfGroups;
        this.numberOfSamples = numberOfSamples;

        refCountsPerSample = new int[numberOfSamples];
        variantsCountPerSample = new int[numberOfSamples];
    }

    IntArrayList decreasingCounts = new IntArrayList();
    CharArraySet baseIndicesByDecreasingCounts = new CharArraySet();
    MutableString genotypeBuffer = new MutableString();

    public void writeRecord(DiscoverVariantIterateSortedAlignments iterator, SampleCountInfo[] sampleCounts,
                            int referenceIndex, int position, ObjectArrayList<PositionBaseInfo> list, int groupIndexA, int groupIndexB) {
     
        fillVariantCountArrays(sampleCounts);
        CharSequence currentReferenceId = iterator.getReferenceId(referenceIndex);

        final int positionString = position + 1;
        statWriter.setValue(refIdColumnIndex, String.format("%s:%d:%s:%d", currentReferenceId, positionString, currentReferenceId, positionString));

        for (int sampleIndex = 0; sampleIndex < numberOfSamples; sampleIndex++) {
            decreasingCounts.size(0);
            baseIndicesByDecreasingCounts.clear();
            
            SampleCountInfo sci = sampleCounts[sampleIndex];
            boolean hasMoreCounts = true;
            while (hasMoreCounts) {
                int maxCount = 0;
                int maxCountIndex = -1;
                for (int baseIndex = 0; baseIndex < SampleCountInfo.BASE_MAX_INDEX; baseIndex++) {
                    final int count = sci.counts[baseIndex];
                    if (count > maxCount) {
                        maxCount = count;
                        maxCountIndex = baseIndex;
                    }
                }
                if (maxCount != 0) {
                    decreasingCounts.add(maxCount);
                    baseIndicesByDecreasingCounts.add(sci.base(maxCountIndex));
                    sci.counts[maxCountIndex] = 0;
                } else {
                    hasMoreCounts = false;
                }
            }
            String zygozity;
            if (decreasingCounts.size() == 1) {
                zygozity = "Homozygous";
            } else if (decreasingCounts.size() == 2) {
                zygozity = "Heterozygous";
            } else {
                zygozity = "Mixture";
            }
            statWriter.setValue(zygozity,
                    "Zygosity[%s]", samples[sampleIndex]);

            final CharArrayList charArrayList = CharArrayList.wrap(baseIndicesByDecreasingCounts.toCharArray());
            Collections.sort(charArrayList);

            genotypeBuffer.setLength(0);
            int i = 0;
            int max = baseIndicesByDecreasingCounts.size()-1;
            for (char base : charArrayList) {
                genotypeBuffer.append(base);
                if (i++ < max) {
                    genotypeBuffer.append('/');
                }
            }

            statWriter.setValue(genotypeBuffer.toString(),
                    "Genotype[%s]", samples[sampleIndex]);
        }

        statWriter.writeRecord();
    }


    private void fillVariantCountArrays(SampleCountInfo[] sampleCounts) {


        for (SampleCountInfo csi : sampleCounts) {
            final int sampleIndex = csi.sampleIndex;
            variantsCountPerSample[sampleIndex] = csi.varCount;
            refCountsPerSample[sampleIndex] = csi.refCount;
        }

    }

}