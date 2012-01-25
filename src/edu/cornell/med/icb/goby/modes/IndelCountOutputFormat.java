/*
 * Copyright (C) 2009-2012 Institute for Computational Biomedicine,
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

import edu.cornell.med.icb.goby.Release1_9_7_2;
import edu.cornell.med.icb.goby.algorithmic.data.EquivalentIndelRegion;
import edu.cornell.med.icb.goby.alignments.DiscoverVariantIterateSortedAlignments;
import edu.cornell.med.icb.goby.alignments.DiscoverVariantPositionData;
import edu.cornell.med.icb.goby.alignments.SampleCountInfo;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * @author Fabien Campagne
 *         Date: 1/25/12
 *         Time: 12:39 PM
 */
public class IndelCountOutputFormat implements SequenceVariationOutputFormat {
    private int[] sampleIndelCounts;
    private int[] groupIndelCounts;
    int[] readerIndexToGroupIndex;
    private String[] sampleIds;
    private String[] groupIds;
    private boolean allocated;

    @Override
    public void defineColumns(PrintWriter statsWriter, DiscoverSequenceVariantsMode mode) {
        readerIndexToGroupIndex = mode.getReaderIndexToGroupIndex();
        sampleIds = mode.getSamples();
        groupIds = mode.getGroups();
        // force Goby to call indels:
        Release1_9_7_2.callIndels = true;
    }

    @Override
    public void allocateStorage(int numberOfSamples, int numberOfGroups) {
        if (!allocated) {
            sampleIndelCounts = new int[numberOfSamples];
            groupIndelCounts = new int[numberOfGroups];
            allocated = true;
        }
    }

    @Override
    public void writeRecord(final DiscoverVariantIterateSortedAlignments iterator, final SampleCountInfo[] sampleCounts,
                            final int referenceIndex, final int position,
                            final DiscoverVariantPositionData list, final int groupIndexA,
                            final int groupIndexB) {

        if (list.hasCandidateIndels()) {
            for (final EquivalentIndelRegion indel : list.getIndels()) {

                if (indel.getFrequency() >= Math.max(5,list.size()/3)) {
                   // System.out.printf("sample %d referenceIndex %d position: %d %s %n", indel.sampleIndex, referenceIndex, position, indel);
                    sampleIndelCounts[indel.sampleIndex]++;
                    final int groupIndex = readerIndexToGroupIndex[indel.sampleIndex];
                    groupIndelCounts[groupIndex]++;
                }
            }
        }
    }

    @Override
    public void close() {
        try {
            PrintWriter output = new PrintWriter(new FileWriter("indel-counts.tsv"));
            int sampleIndex = 0;
            for (String sample : sampleIds) {
                output.write(String.format("SAMPLE\t%s\t%d%n", sample, sampleIndelCounts[sampleIndex]));
                sampleIndex++;
            }
            int groupIndex = 0;
            for (String group : groupIds) {
                output.write(String.format("GROUP\t%s\t%d%n", group, groupIndelCounts[readerIndexToGroupIndex[groupIndex]]));
                groupIndex++;
            }
            output.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Override
    public void setGenome(RandomAccessSequenceInterface genome) {

    }
}
