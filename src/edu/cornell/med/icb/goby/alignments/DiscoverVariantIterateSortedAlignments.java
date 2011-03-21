/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 *  This file is part of the Goby IO API.
 *
 *     The Goby IO API is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     The Goby IO API is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with the Goby IO API.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.cornell.med.icb.goby.alignments;

import edu.cornell.med.icb.goby.modes.DiscoverSequenceVariantsMode;
import edu.cornell.med.icb.goby.modes.SequenceVariationOutputFormat;
import edu.cornell.med.icb.goby.stats.StatisticsWriter;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.PrintWriter;

/**
 * Helper class to implement the logic of discovering sequence variations in and across groups of samples.
 * Implements most of the work done by DiscoverSequenceVariantsMode
 *
 * @author Fabien Campagne
 *         Date: Sep 7, 2010
 *         Time: 2:14:38 PM
 * @see edu.cornell.med.icb.goby.modes.DiscoverSequenceVariantsMode
 */
public class DiscoverVariantIterateSortedAlignments
        extends IterateSortedAlignmentsListImpl {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(DiscoverVariantIterateSortedAlignments.class);

    private int thresholdDistinctReadIndices = 10;
    private int minimumVariationSupport = 3;

    private SequenceVariationOutputFormat format;
    private int numberOfGroups;
    private int[] readerIndexToGroupIndex;

    public void setMinimumVariationSupport(int minimumVariationSupport) {
        this.minimumVariationSupport = minimumVariationSupport;
    }

    public void setThresholdDistinctReadIndices(int thresholdDistinctReadIndices) {
        this.thresholdDistinctReadIndices = thresholdDistinctReadIndices;
    }

    public DiscoverVariantIterateSortedAlignments(SequenceVariationOutputFormat format) {
        this.format = format;
    }


    StatisticsWriter statWriter;
    String[] samples;

    public void initialize(DiscoverSequenceVariantsMode mode,
                           PrintWriter outWriter) {
        readerIndexToGroupIndex = mode.getReaderIndexToGroupIndex();
        statWriter = new StatisticsWriter(outWriter);
        format.defineColumns(statWriter, mode);

    }

    public void finish() {
        statWriter.close();

    }

    public class PositionBaseInfo {
        public int readIndex;
        public int readerIndex;
        public byte qualityScore;
        public boolean matchesReference;
        public char from;
        public char to;
        public int position;
    }

    public void allocateStorage(int numberOfSamples, int numberOfGroups) {
        this.numberOfSamples = numberOfSamples;
        this.numberOfGroups = numberOfGroups;
        sampleCounts = new SampleCountInfo[numberOfSamples];
        for (int i = 0; i < numberOfSamples; i++) {
            sampleCounts[i] = new SampleCountInfo();
        }
        format.allocateStorage(numberOfSamples, numberOfGroups);
    }

    public class SampleCountInfo {
        public int countBaseA;
        public int countBaseT;
        public int countBaseC;
        public int countBaseG;
        public int countBaseOther;
        public char referenceBase;
        public IntSet distinctReadIndices = new IntArraySet();
        public int sampleIndex;
        public int varCount;
        public int refCount;
    }

    private SampleCountInfo[] sampleCounts;

    private int numberOfSamples;


    public void processPositions(int referenceIndex, int position,
                                 ObjectArrayList<IterateSortedAlignmentsListImpl.PositionBaseInfo> list) {
        int sumVariantCounts = 0;


        for (int sampleIndex = 0; sampleIndex < numberOfSamples; sampleIndex++) {
            sampleCounts[sampleIndex].countBaseA = 0;
            sampleCounts[sampleIndex].countBaseT = 0;
            sampleCounts[sampleIndex].countBaseC = 0;
            sampleCounts[sampleIndex].countBaseG = 0;
            sampleCounts[sampleIndex].countBaseOther = 0;
            sampleCounts[sampleIndex].referenceBase = '?';
            sampleCounts[sampleIndex].distinctReadIndices.clear();
            sampleCounts[sampleIndex].sampleIndex = sampleIndex;
            sampleCounts[sampleIndex].varCount = 0;
            sampleCounts[sampleIndex].refCount = 0;
        }

        if (list != null) {
            IntSet distinctReadIndices = new IntArraySet();
            for (IterateSortedAlignmentsListImpl.PositionBaseInfo info : list) {
                final int sampleIndex = info.readerIndex;
                distinctReadIndices.add(info.readIndex);
                if (info.matchesReference) {

                    sampleCounts[sampleIndex].referenceBase = info.from;
                    sampleCounts[sampleIndex].refCount++;

                } else {
                    sampleCounts[sampleIndex].varCount++;
                    sumVariantCounts++;

                    sampleCounts[sampleIndex].distinctReadIndices.add(info.readIndex);
                    switch (info.to) {
                        case 'A':
                            sampleCounts[sampleIndex].countBaseA += 1;
                            break;
                        case 'T':
                            sampleCounts[sampleIndex].countBaseT += 1;
                            break;
                        case 'C':
                            sampleCounts[sampleIndex].countBaseC += 1;
                            break;
                        case 'G':
                            sampleCounts[sampleIndex].countBaseG += 1;
                            break;
                        default:
                            sampleCounts[sampleIndex].countBaseOther += 1;
                            break;
                    }
                }
            }

            if (distinctReadIndices.size() >= thresholdDistinctReadIndices && sumVariantCounts > minimumVariationSupport) {
                int groupIndexA = 0;
                int groupIndexB = 1;
                // Do not write statistics for positions in the start flap. The flap start is used to accumulate
                // base counts for reads that can overlap with the window under consideration.

                if (!isWithinStartFlap(referenceIndex, position)) {
                    format.writeRecord(this, sampleCounts, referenceIndex, position, list, groupIndexA, groupIndexB);
                }
            }
        }
    }
}