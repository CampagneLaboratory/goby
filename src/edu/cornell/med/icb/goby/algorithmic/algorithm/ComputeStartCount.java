/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                         Weill Medical College of Cornell University
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

package edu.cornell.med.icb.goby.algorithmic.algorithm;

import edu.cornell.med.icb.goby.counts.CountsWriter;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntSet;

import java.io.IOException;
import java.util.Collections;

/**
 * Support to write read start positions in a count writer. Strands of the reference sequence can be inspected together
 * or separately.
 *
 * @author Fabien Campagne
 *         Date: Jun 4, 2009
 *         Time: 6:37:41 PM
 */
public class ComputeStartCount extends ComputeCount {
    /**
     * Indicates which strand(s) to collect.
     */
    int focusOnStrand;

    // TODO: Change to enum
    public static int POSITIVE_STRAND_ONLY = 0;
    public static int REVERSE_STRAND_ONLY = 1;
    public static int BOTH_STRAND = 2;

    /**
     * @param focusOnStrand Indicates which strand(s) should be inspected.
     */
    public ComputeStartCount(final int focusOnStrand) {
        super();
        this.focusOnStrand = focusOnStrand;
        starts = new Int2IntOpenHashMap();
    }

    /**
     * Number of reads that start on this position on the forward, reverse, or on both strands.
     */
    protected Int2IntMap starts;

    @Override
    public void accumulate() {
        // do nothing in this implementation.
    }

    @Override
    public void populate(final int startPosition, final int endPosition, final boolean forwardStrand) {

        //    assert forwardStrand && startPosition <= endPosition ||
        //           !forwardStrand && endPosition <= startPosition : "start (resp. end) must occur before end (start) on forward (reverse) strand.";
        final int beginPosition;

        if (!forwardStrand) {
            beginPosition = endPosition;

        } else {
            beginPosition = startPosition;

        }
        if (forwardStrand && focusOnStrand == BOTH_STRAND || focusOnStrand == POSITIVE_STRAND_ONLY) {
            final int count = starts.get(beginPosition) + 1;
            starts.put(beginPosition, count);

        } else if (!forwardStrand && (focusOnStrand == BOTH_STRAND || focusOnStrand == REVERSE_STRAND_ONLY)) {

            final int count = starts.get(beginPosition) + 1;
            starts.put(beginPosition, count);

        }
    }

    /**
     * Accumulate the number of reads that start at a given position. Write the resulting histogram to the counts writer.
     *
     * @param writer
     * @throws java.io.IOException
     */
    @Override
    public void baseCount(final CountsWriter writer) throws IOException {


        final IntSet startPositions = starts.keySet();
        if (starts.size() == 0) {
            writer.close();
            return;
        }
        final IntList sortedStartPositions = new IntArrayList();
        sortedStartPositions.addAll(startPositions);
        Collections.sort(sortedStartPositions);

        final int firstPosition = sortedStartPositions.get(0);
        writer.appendCount(0, firstPosition);
        int lengthConstant = 1;
        int count;
        int prevCount = 0;

        int currentIndex = 0;
        final int maxPosition = Integer.MAX_VALUE;
        int startPosition;
        for (currentIndex = 0; currentIndex < sortedStartPositions.size(); ++currentIndex) {

            startPosition = sortedStartPositions.get(currentIndex);
            count = starts.get(startPosition);
            //determine the stretch of constant values after this start position:
            lengthConstant = 1;
            for (int i = startPosition + 1; i < maxPosition; ++i) {
                if (starts.get(i) == count) {
                    lengthConstant++;
                }
                else {
                    currentIndex += lengthConstant - 1;
                    break;
                }
            }

            writer.appendCount(count, lengthConstant);
            prevCount = count;
            //
            if (starts.get(startPosition + lengthConstant) == 0) {
                // going back to zero
                // for how many positions?
                final int nextStart;
                final int length;
                if (currentIndex + lengthConstant >= sortedStartPositions.size()) {
                    // past the end, just append one more zero.
                    length = 1;
                } else {
                    // end not reached, find out how many positions are filled with zero between
                    // here and the next count!=0 position.
                    nextStart = sortedStartPositions.get(currentIndex + lengthConstant);
                    length = nextStart - startPosition - 1;
                }
                writer.appendCount(0, length);
                prevCount = 0;
            }

        }

        writer.close();
    }
}
