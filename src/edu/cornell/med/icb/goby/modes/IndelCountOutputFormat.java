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
import edu.cornell.med.icb.goby.algorithmic.data.MethylCountInfo;
import edu.cornell.med.icb.goby.alignments.DiscoverVariantIterateSortedAlignments;
import edu.cornell.med.icb.goby.alignments.DiscoverVariantPositionData;
import edu.cornell.med.icb.goby.alignments.SampleCountInfo;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

/**
 * Format that estimates the number of called indels per hundred thousands bases observed.
 *
 * @author Fabien Campagne
 *         Date: 1/25/12
 *         Time: 12:39 PM
 */
public class IndelCountOutputFormat implements SequenceVariationOutputFormat {
    private int[] sampleIndelCounts;
    private int[] groupIndelCounts;
    private long[] sampleSitesObserved;
    private long[] groupSitesObserved;
    int[] readerIndexToGroupIndex;
    private String[] sampleIds;
    private String[] groupIds;
    private boolean allocated;
    private int minPosition = Integer.MAX_VALUE;
    private int minRefIndex = Integer.MAX_VALUE;
    private int maxPosition;
    private int maxRefIndex;
    private CharSequence minRefId;
    private int previousMaxRefIndex = -1;
    private CharSequence maxRefId;
    private static final int MIN_COVERAGE_THRESHOLD = 5;
    private boolean headerWritten;
    private PrintWriter output;
    private int numSites;
    private int maxSitesPerAccumulation = 10000000;

    @Override
    public void defineColumns(PrintWriter statsWriter, DiscoverSequenceVariantsMode mode) {
        readerIndexToGroupIndex = mode.getReaderIndexToGroupIndex();
        sampleIds = mode.getSamples();
        groupIds = mode.getGroups();
        output = statsWriter;
        assert mode.getCallIndels() : "indel calling must be active.";

    }

    @Override
    public void allocateStorage(int numberOfSamples, int numberOfGroups) {
        if (!allocated) {
            sampleIndelCounts = new int[numberOfSamples];
            groupIndelCounts = new int[numberOfGroups];
            sampleSitesObserved = new long[numberOfSamples];
            groupSitesObserved = new long[numberOfGroups];
            allocated = true;
        }
    }

    @Override
    public void writeRecord(final DiscoverVariantIterateSortedAlignments iterator, final SampleCountInfo[] sampleCounts,
                            final int referenceIndex, final int position,
                            final DiscoverVariantPositionData list, final int groupIndexA,
                            final int groupIndexB) {
        writeHeader();
        minPosition = Math.min(position, minPosition);
        minRefIndex = Math.min(referenceIndex, minRefIndex);
        maxPosition = Math.max(position, maxPosition);
        maxRefIndex = Math.max(referenceIndex, maxRefIndex);
        if (maxRefIndex != previousMaxRefIndex || maxRefId == null) {
            maxRefId = iterator.getReferenceId(maxRefIndex);
        }
        previousMaxRefIndex = maxRefIndex;
        if (minRefId == null) {
            minRefId = iterator.getReferenceId(minRefIndex);
        }
        for (SampleCountInfo sci : sampleCounts) {
            int totalCount = 0;
            for (int count : sci.counts) {
                totalCount += count;
            }

            if (totalCount >= MIN_COVERAGE_THRESHOLD) {
                sampleSitesObserved[sci.sampleIndex]++;
                groupSitesObserved[readerIndexToGroupIndex[sci.sampleIndex]]++;
            }
        }
        if (list.hasCandidateIndels()) {
            for (final EquivalentIndelRegion indel : list.getIndels()) {

                if (indel.getFrequency() >= Math.max(MIN_COVERAGE_THRESHOLD, list.size() / 3)) {
                    // frequency must be at least 5 or a third of the number of bases at position, whichever is smaller.
                    // System.out.printf("sample %d referenceIndex %d position: %d %s %n", indel.sampleIndex, referenceIndex, position, indel);
                    sampleIndelCounts[indel.sampleIndex]++;
                    final int groupIndex = readerIndexToGroupIndex[indel.sampleIndex];
                    groupIndelCounts[groupIndex]++;
                }
            }
        }
        if (numSites++ > maxSitesPerAccumulation) {
            flushToDisk();
            numSites = 0;
            minPosition = -1;
            maxPosition = -1;
            minRefIndex = -1;
            maxRefIndex = -1;
            minRefId = null;
            maxRefId = null;
        }
    }

    @Override
    public void close() {


        writeHeader();
        flushToDisk();
        output.close();

    }

    private void flushToDisk() {

        int sampleIndex = 0;
        for (String sample : sampleIds) {
            output.write(String.format("SAMPLE\t%s\t%s:%d-%s:%d\t%d\t%d\t%g%n", sample,
                    minRefId, minPosition, maxRefId, maxPosition,
                    sampleIndelCounts[sampleIndex],
                    sampleSitesObserved[sampleIndex],
                    100000d * fraction(sampleIndelCounts[sampleIndex], sampleSitesObserved[sampleIndex])));
            sampleIndex++;
        }
        int groupIndex = 0;
        for (String group : groupIds) {

            output.write(String.format("GROUP\t%s\t%s:%d-%s:%d\t%d\t%d\t%g%n", group,
                    minRefId, minPosition, maxRefId, maxPosition,
                    groupIndelCounts[groupIndex],
                    groupSitesObserved[groupIndex],
                    100000d * fraction(groupIndelCounts[groupIndex], groupSitesObserved[groupIndex])));
            groupIndex++;
        }
        output.flush();
        Arrays.fill(sampleIndelCounts, 0);
        Arrays.fill(groupIndelCounts, 0);
        Arrays.fill(sampleSitesObserved, 0);
        Arrays.fill(groupSitesObserved, 0);
    }

    private void writeHeader() {
        if (!headerWritten) {

            output.write("STAT-TYPE\tlabel\tslice-id\tindel-count\t#sites-observed\tindels/100k-bases\n");
            headerWritten = true;

        }
    }


    private double fraction(int a, long b) {
        return ((double) a) / ((double) b);
    }

    @Override
    public void setGenome(RandomAccessSequenceInterface genome) {

    }

    public static void fillMethylationCountArrays(SampleCountInfo[] sampleCounts, DiscoverVariantPositionData list, int position, char refBase, MethylCountInfo mci, int[] readerIndexToGroupIndex) {
        // don't use threshold on events at site for indel rates:
        mci.eventCountAtSite=0;
        for (SampleCountInfo sci : sampleCounts) {
            int totalCount = 0;
            for (int count : sci.counts) {
                totalCount += count;
            }

            if (totalCount >= MIN_COVERAGE_THRESHOLD) {

                mci.unmethylatedCCountPerSample[sci.sampleIndex]++;
                mci.unmethylatedCCountsPerGroup[readerIndexToGroupIndex[sci.sampleIndex]]++;
                mci.eventCountAtSite+=totalCount;
                if (list.hasCandidateIndels()) {
                    for (final EquivalentIndelRegion indel : list.getIndels()) {

                        if (indel.getFrequency() >= Math.max(MIN_COVERAGE_THRESHOLD, list.size() / 3)) {
                            // frequency must be at least 5 or a third of the number of bases at position, whichever is smaller.
                            // System.out.printf("sample %d referenceIndex %d position: %d %s %n", indel.sampleIndex, referenceIndex, position, indel);

                            final int groupIndex = readerIndexToGroupIndex[indel.sampleIndex];
                            mci.methylatedCCountsPerSample[sci.sampleIndex]++;
                            mci.methylatedCCountPerGroup[groupIndex]++;
                        }
                    }
                }
            }
        }
    }
}
