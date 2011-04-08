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

import edu.cornell.med.icb.goby.stats.TSVWriter;
import edu.cornell.med.icb.goby.stats.VCFWriter;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectIterator;
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

    private BaseFilter[] baseFilters;

    public void setMinimumVariationSupport(int minimumVariationSupport) {
        this.minimumVariationSupport = minimumVariationSupport;
    }

    public void setThresholdDistinctReadIndices(int thresholdDistinctReadIndices) {
        this.thresholdDistinctReadIndices = thresholdDistinctReadIndices;
    }

    public DiscoverVariantIterateSortedAlignments(SequenceVariationOutputFormat format) {
        this.format = format;
    }


    Object statWriter;
    String[] samples;

    public void initialize(DiscoverSequenceVariantsMode mode,
                           PrintWriter outWriter) {
        readerIndexToGroupIndex = mode.getReaderIndexToGroupIndex();

        if (mode.outputVCF()) {
            statWriter = new VCFWriter(outWriter);
        } else {
            statWriter = new TSVWriter(outWriter);
        }

        format.defineColumns(outWriter, mode);

        if (mode.getDiffExpAnalyzer().eval("filter")) {
            baseFilters = new BaseFilter[]{

                    new QualityScoreFilter(),
                    //     new FisherBaseFilter(mode.getReadIndexStats()),
                    new LeftOverFilter()
            };
            System.out.println("Filtering reads that have these criteria:");
            for (BaseFilter filter : baseFilters) {

                System.out.println(filter.describe());
            }

        } else {
            baseFilters = new BaseFilter[0];
        }
    }

    public void finish() {
        format.close();

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

    private SampleCountInfo[] sampleCounts;

    private int numberOfSamples;


    public void processPositions(int referenceIndex, int position,
                                 ObjectArrayList<edu.cornell.med.icb.goby.alignments.PositionBaseInfo> list) {
        int sumVariantCounts = 0;


        for (int sampleIndex = 0; sampleIndex < numberOfSamples; sampleIndex++) {
            sampleCounts[sampleIndex].counts[SampleCountInfo.BASE_A_INDEX] = 0;
            sampleCounts[sampleIndex].counts[SampleCountInfo.BASE_T_INDEX] = 0;
            sampleCounts[sampleIndex].counts[SampleCountInfo.BASE_C_INDEX] = 0;
            sampleCounts[sampleIndex].counts[SampleCountInfo.BASE_G_INDEX] = 0;
            sampleCounts[sampleIndex].counts[SampleCountInfo.BASE_OTHER_INDEX] = 0;
            sampleCounts[sampleIndex].referenceBase = '?';
            sampleCounts[sampleIndex].distinctReadIndices.clear();
            sampleCounts[sampleIndex].sampleIndex = sampleIndex;
            sampleCounts[sampleIndex].varCount = 0;
            sampleCounts[sampleIndex].refCount = 0;
            sampleCounts[sampleIndex].failedCount = 0;
        }

        if (list != null) {
            IntSet distinctReadIndices = new IntArraySet();
            char refBase = setReferenceAllele(list);
            for (edu.cornell.med.icb.goby.alignments.PositionBaseInfo info : list) {
                final int sampleIndex = info.readerIndex;
                distinctReadIndices.add(info.readIndex);
                if (info.matchesReference) {

                    sampleCounts[sampleIndex].referenceBase = info.from;
                    sampleCounts[sampleIndex].refCount++;
                    incrementBaseCounter(info.from, sampleIndex);

                } else {
                    sampleCounts[sampleIndex].varCount++;
                    sumVariantCounts++;
                    sampleCounts[sampleIndex].referenceBase = refBase;
                    sampleCounts[sampleIndex].distinctReadIndices.add(info.readIndex);
                    incrementBaseCounter(info.to, sampleIndex);
                }
            }

            if (distinctReadIndices.size() >= thresholdDistinctReadIndices && sumVariantCounts > minimumVariationSupport) {
                int groupIndexA = 0;
                int groupIndexB = 1;
                // Do not write statistics for positions in the start flap. The flap start is used to accumulate
                // base counts for reads that can overlap with the window under consideration.

                if (!isWithinStartFlap(referenceIndex, position)) {

                    if (baseFilters.length != 0) {
                        ObjectArrayList<edu.cornell.med.icb.goby.alignments.PositionBaseInfo> filteredList =
                                new ObjectArrayList<edu.cornell.med.icb.goby.alignments.PositionBaseInfo>();
                        for (BaseFilter filter : baseFilters) {
                            filter.filterBases(list, sampleCounts, filteredList);
                            //       System.out.printf("filter %s removed %3g %% %n", filter.getName(), filter.getPercentFilteredOut());
                        }
                        CountFixer fixer = new CountFixer();
                        fixer.fix(list, sampleCounts, filteredList);
                    }
                }
                format.writeRecord(this, sampleCounts, referenceIndex, position, list, groupIndexA, groupIndexB);
            }
        }
    }

    private char setReferenceAllele(ObjectArrayList<edu.cornell.med.icb.goby.alignments.PositionBaseInfo> list) {
        final ObjectIterator<edu.cornell.med.icb.goby.alignments.PositionBaseInfo> iterator = list.iterator();
        char refBase = '\0';
        // find the reference base from any variant:
        while (iterator.hasNext()) {
            edu.cornell.med.icb.goby.alignments.PositionBaseInfo positionBaseInfo = iterator.next();
            if (!positionBaseInfo.matchesReference) {
                refBase = positionBaseInfo.from;
                break;
            }
        }
        // set on elements that match the reference:
        for (edu.cornell.med.icb.goby.alignments.PositionBaseInfo elem : list) {
            if (elem.matchesReference) {
                elem.from = refBase;
            }
        }
        return refBase;
    }

    private void incrementBaseCounter(char base, int sampleIndex) {
        switch (base) {
            case 'A':
                sampleCounts[sampleIndex].counts[SampleCountInfo.BASE_A_INDEX] += 1;
                break;
            case 'T':
                sampleCounts[sampleIndex].counts[SampleCountInfo.BASE_T_INDEX] += 1;
                break;
            case 'C':
                sampleCounts[sampleIndex].counts[SampleCountInfo.BASE_C_INDEX] += 1;
                break;
            case 'G':
                sampleCounts[sampleIndex].counts[SampleCountInfo.BASE_G_INDEX] += 1;
                break;
            default:
                sampleCounts[sampleIndex].counts[SampleCountInfo.BASE_OTHER_INDEX] += 1;
                break;
        }
    }


}