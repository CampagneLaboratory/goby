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

import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.apache.log4j.Logger;

/**
 * @author Fabien Campagne
 *         Date: Sep 7, 2010
 *         Time: 2:14:38 PM
 */
public abstract class IterateSortedAlignmentsListImpl
        extends IterateSortedAlignments<ObjectArrayList<PositionBaseInfo>> {
    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(IterateSortedAlignmentsListImpl.class);

    public abstract void processPositions(int referenceIndex, int intermediatePosition, ObjectArrayList<PositionBaseInfo> positionBaseInfos);

    public void observeReferenceBase(ConcatSortedAlignmentReader sortedReaders,
                                     Alignments.AlignmentEntry alignmentEntry,
                                     Int2ObjectMap<ObjectArrayList<PositionBaseInfo>> positionToBases,
                                     int currentReferenceIndex, int currentRefPosition, int currentReadIndex) {
        if (LOG.isDebugEnabled()) {
            LOG.debug(String.format("RB: queryIndex=%d\tref_position=%d\tread_index=%d",
                alignmentEntry.getQueryIndex(), currentRefPosition, currentReadIndex));
        }

        PositionBaseInfo info = new PositionBaseInfo();

        info.readerIndex = sortedReaders.activeIndex;
   //     System.out.printf("observing ref readerIndex=%d%n",info.readerIndex);
        info.readIndex = currentReadIndex;
        info.from = '\0';
        info.to = '\0';
        info.matchesReference = true;
        info.position = currentRefPosition;
        info.qualityScore = 40;

        addToFuture(positionToBases, info);
    }


    public void observeVariantBase(ConcatSortedAlignmentReader sortedReaders,
                                   Alignments.AlignmentEntry alignmentEntry,
                                   Int2ObjectMap<ObjectArrayList<PositionBaseInfo>> positionToBases,
                                   Alignments.SequenceVariation var,
                                   char toChar, char fromChar,
                                   byte toQual, int currentReferenceIndex,
                                   int currentRefPosition,
                                   int currentReadIndex) {

        if (LOG.isDebugEnabled()) {
            LOG.debug(String.format("VB: queryIndex=%d\tref_position=%d\tread_index=%d\tfromChar=%c\ttoChar=%c",
                    alignmentEntry.getQueryIndex(), currentRefPosition, currentReadIndex, fromChar, toChar));
        }

        PositionBaseInfo info = new PositionBaseInfo();
        info.readerIndex = sortedReaders.activeIndex;
    //    System.out.printf("observing var readerIndex=%d%n",info.readerIndex);

        info.readIndex = currentReadIndex;
        info.from = fromChar;
        info.to = toChar;
        info.matchesReference = false;
        info.position = currentRefPosition;
        info.qualityScore=toQual;

        addToFuture(positionToBases, info);
    }


    private void addToFuture(Int2ObjectMap<ObjectArrayList<PositionBaseInfo>> positionToBases,
                             PositionBaseInfo info) {

        ObjectArrayList<PositionBaseInfo> list = positionToBases.get(info.position);
        if (list == null) {
            list = new ObjectArrayList<PositionBaseInfo>();
            positionToBases.put(info.position, list);
        }
        list.add(info);
    }


}
