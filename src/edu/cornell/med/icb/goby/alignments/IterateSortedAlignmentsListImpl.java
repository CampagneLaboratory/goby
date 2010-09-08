/*
 * Copyright (C) 2010 Institute for Computational Biomedicine,
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

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;

/**
 * @author Fabien Campagne
 *         Date: Sep 7, 2010
 *         Time: 2:14:38 PM
 */
public abstract class IterateSortedAlignmentsListImpl extends IterateSortedAlignments<ObjectArrayList<IterateSortedAlignmentsListImpl.PositionBaseInfo>> {

    public abstract void processPositions(int position, ObjectArrayList<PositionBaseInfo> positionBaseInfos);

    public void observeReferenceBase(ConcatSortedAlignmentReader sortedReaders, Alignments.AlignmentEntry alignmentEntry,
                                     Int2ObjectMap<ObjectArrayList<PositionBaseInfo>> positionToBases,
                                     int currentRefPosition, int currentReadIndex) {
        PositionBaseInfo info = new PositionBaseInfo();
        info.readerIndex = sortedReaders.activeIndex;
        info.readIndex = currentReadIndex;
        info.from = '\0';
        info.to = '\0';
        info.matchesReference = true;
        info.position = currentRefPosition;
        info.qualityScore = 40;

        addToFuture(positionToBases, info);
    }


    public void observeVariantBase(ConcatSortedAlignmentReader sortedReaders,
                                   Int2ObjectMap<ObjectArrayList<PositionBaseInfo>> positionToBases,
                                   Alignments.SequenceVariation var,
                                   char toChar, char fromChar, int currentRefPosition, int currentReadIndex) {

        PositionBaseInfo info = new PositionBaseInfo();
        info.readerIndex = sortedReaders.activeIndex;
        info.readIndex = currentReadIndex;
        info.from = fromChar;
        info.to = toChar;
        info.matchesReference = false;
        info.position =  currentRefPosition;


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

    public class PositionBaseInfo {
        public int readIndex;
        public int readerIndex;
        public byte qualityScore;
        public boolean matchesReference;
        public char from;
        public char to;
        public int position;
    }

}
