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

import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.ints.IntArraySet;

import java.util.Map;

import org.apache.log4j.Logger;

/**
 * Created by IntelliJ IDEA.
 * User: kdorff
 * Date: Apr 13, 2011
 * Time: 12:56:01 PM
 * To change this template use File | Settings | File Templates.
 */
public class IterateSortedAlignmentsTester extends IterateSortedAlignments<Object> {

    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(IterateSortedAlignmentsTester.class);

    public PositionToBasesMap<PerQueryAlignmentData> queryIndexToAlignmentDataMap;

    public IterateSortedAlignmentsTester() {
        queryIndexToAlignmentDataMap = new PositionToBasesMap<PerQueryAlignmentData>();
    }

    @Override
    public void observeReferenceBase(ConcatSortedAlignmentReader sortedReaders, Alignments.AlignmentEntry alignmentEntry,
                                     PositionToBasesMap<Object> positionToBases,
                                     int currentReferenceIndex, int currentRefPosition, int currentReadIndex) {
        if (LOG.isDebugEnabled()) {
            LOG.debug(String.format("RB: queryIndex=%d\tref_position=%d\tread_index=%d",
                alignmentEntry.getQueryIndex(), currentRefPosition, currentReadIndex));
        }
        int queryIndex = alignmentEntry.getQueryIndex();
        if (currentReadIndex >= 1) {
            PerQueryAlignmentData alignmentData = queryIndexToAlignmentDataMap.get(queryIndex);
            if (alignmentData == null) {
                alignmentData = new PerQueryAlignmentData();
                queryIndexToAlignmentDataMap.put(queryIndex, alignmentData);
            }
            if (alignmentData.firstReadIndex == -1) {
                alignmentData.firstReadIndex = currentReadIndex;
                alignmentData.queryPosition = alignmentEntry.getQueryPosition();
                alignmentData.targetPosition = alignmentEntry.getPosition();
                alignmentData.queryLength = alignmentEntry.getQueryLength();
                alignmentData.queryAlignedLength = alignmentEntry.getQueryAlignedLength();
                alignmentData.targetAlignedLength = alignmentEntry.getTargetAlignedLength();
                alignmentData.reverseStrand = alignmentEntry.getMatchingReverseStrand();
            }
            alignmentData.observe(currentRefPosition, currentReadIndex);
        } else {
            throw new RuntimeException(String.format("queryIndex=%d readIndex=%d should be >=1",
                    queryIndex, currentReadIndex));
        }
    }

    @Override
    public void observeVariantBase(
            ConcatSortedAlignmentReader sortedReaders, Alignments.AlignmentEntry alignmentEntry,
            PositionToBasesMap<Object> positionToBases, Alignments.SequenceVariation var,
            char toChar, char fromChar, byte toQual, int currentReferenceIndex,
            int currentRefPosition, int currentReadIndex) {
        if (LOG.isDebugEnabled()) {
            LOG.debug(String.format("VB: queryIndex=%d\tref_position=%d\tread_index=%d\tfromChar=%c\ttoChar=%c",
                    alignmentEntry.getQueryIndex(), currentRefPosition, currentReadIndex, fromChar, toChar));
        }

        int queryIndex = alignmentEntry.getQueryIndex();
        if (currentReadIndex >= 1) {
            PerQueryAlignmentData alignmentData = queryIndexToAlignmentDataMap.get(queryIndex);
            if (alignmentData == null) {
                alignmentData = new PerQueryAlignmentData();
                queryIndexToAlignmentDataMap.put(queryIndex, alignmentData);
            }
            if (alignmentData.firstReadIndex == -1) {
                alignmentData.firstReadIndex = currentReadIndex;
                alignmentData.queryPosition = alignmentEntry.getQueryPosition();
                alignmentData.targetPosition = alignmentEntry.getPosition();
                alignmentData.queryLength = alignmentEntry.getQueryLength();
                alignmentData.queryAlignedLength = alignmentEntry.getQueryAlignedLength();
                alignmentData.targetAlignedLength = alignmentEntry.getTargetAlignedLength();
                alignmentData.reverseStrand = alignmentEntry.getMatchingReverseStrand();
            }
            alignmentData.observe(currentRefPosition, currentReadIndex, fromChar, toChar);
        } else {
            throw new RuntimeException(String.format("queryIndex=%d readIndex=%d should be >=1",
                    queryIndex, currentReadIndex));
        }
    }

    @Override
    public void processPositions(int referenceIndex, int intermediatePosition, Object positionBaseInfos) {
    }

    /**
     * Remove any items from the map that don't have sequence variations.
     */
    public void removeWithoutSeqvars() {
        IntSet toRemoveQueryIndexes = new IntArraySet(queryIndexToAlignmentDataMap.size());
        for (Map.Entry<Integer, PerQueryAlignmentData> entry : queryIndexToAlignmentDataMap.entrySet()) {
            if (entry.getValue().refPositionReadIndexToBaseMap.size() == 0) {
                toRemoveQueryIndexes.add(entry.getKey());
            }
        }
        for (int queryIndex : toRemoveQueryIndexes) {
            queryIndexToAlignmentDataMap.remove(queryIndex);
        }
    }
}
