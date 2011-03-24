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
import it.unimi.dsi.fastutil.objects.ObjectArrayList;

/**
 * @author Fabien Campagne
 *         Date: Jan 28 2011
 *         Time: 11:53 AM
 */
public abstract class IterateSortedAlignmentsListImpl2
        extends IterateSortedAlignments<ObjectArrayList<IterateSortedAlignmentsListImpl2.PositionBaseInfo2>> {

    public abstract void processPositions(int referenceIndex, int intermediatePosition, ObjectArrayList<PositionBaseInfo2> positionBaseInfos);

    public void observeReferenceBase(ConcatSortedAlignmentReader sortedReaders,
                                     Alignments.AlignmentEntry alignmentEntry,
                                     Int2ObjectMap<ObjectArrayList<PositionBaseInfo2>> positionToBases,
                                     int currentReferenceIndex, int currentRefPosition, int currentReadIndex) {
        PositionBaseInfo2 info = new PositionBaseInfo2();
        info.readerIndex = sortedReaders.activeIndex;
        info.readIndex = currentReadIndex;
        info.from = '\0';
        info.to = '\0';
        info.matchesReference = true;
        info.position = currentRefPosition;
        info.qualityScore = 40;
        info.alignmentEntryQueryIndex = alignmentEntry.getQueryIndex();

        addToFuture(positionToBases, info);
    }


    public void observeVariantBase(ConcatSortedAlignmentReader sortedReaders,
                                   Alignments.AlignmentEntry alignmentEntry, Int2ObjectMap<ObjectArrayList<PositionBaseInfo2>> positionToBases,
                                   Alignments.SequenceVariation var,
                                   char toChar, char fromChar, byte toQual, int currentReferenceIndex, int currentRefPosition, int currentReadIndex) {

        PositionBaseInfo2 info = new PositionBaseInfo2();
        info.readerIndex = sortedReaders.activeIndex;
        info.readIndex = currentReadIndex;
        info.from = fromChar;
        info.to = toChar;
        info.matchesReference = false;
        info.position = currentRefPosition;
        info.alignmentEntryQueryIndex = alignmentEntry.getQueryIndex();
        info.variationLength=Math.max(var.getFrom().length(), var.getTo().length());
        addToFuture(positionToBases, info);
    }


    private void addToFuture(Int2ObjectMap<ObjectArrayList<IterateSortedAlignmentsListImpl2.PositionBaseInfo2>> positionToBases,
                             PositionBaseInfo2 info) {

        ObjectArrayList<PositionBaseInfo2> list = positionToBases.get(info.position);
        if (list == null) {
            list = new ObjectArrayList<PositionBaseInfo2>();
            positionToBases.put(info.position, list);
        }
        list.add(info);
    }

    public class PositionBaseInfo2 {
        public int readIndex;
        public int readerIndex;
        public byte qualityScore;
        public boolean matchesReference;
        public char from;
        public char to;
        public int position;
        public int alignmentEntryQueryIndex;
        public int variationLength;
    }

}