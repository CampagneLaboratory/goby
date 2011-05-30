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

package edu.cornell.med.icb.goby.counts;


import it.unimi.dsi.fastutil.ints.IntArraySet;

import java.io.IOException;
import java.util.Arrays;
import java.util.NoSuchElementException;

/**
 * Skeleton for algorithm development discussion.
 *
 * @author Fabien Campagne
 *         Date: 5/26/11
 *         Time: 10:14 PM
 */
public class UnionAlgorithmSkeleton implements CountsAggregatorI {
    private int numReaders;
    private CountsReaderI[] readers;
    private boolean hasNextTransition;
    private int length;
    private int position = 0;
    private int first;
    private int second;
    private int[] positions;
    private int[] startPositions;
    private int[] endPositions;
    private boolean[] finished;
    private int previousPosition = -1;

    public UnionAlgorithmSkeleton(CountsReaderI... readers) {
        numReaders = readers.length;
        this.readers = readers;
        counts = new int[this.numReaders];
        positions = new int[this.numReaders];
        startPositions = new int[this.numReaders];
        endPositions = new int[this.numReaders];
        Arrays.fill(startPositions, Integer.MAX_VALUE);
        Arrays.fill(endPositions, Integer.MAX_VALUE);
        finished = new boolean[numReaders];
    }

    public int getPosition() {
        return position;
    }

    int counts[];
    IntArraySet startAndEndPositions = new IntArraySet();

    public boolean hasNextTransition() throws IOException {
        if (hasNextTransition) {
            return true;
        }
        hasNextTransition = false;
        position = first;
        first = first(startAndEndPositions);

        for (int readerIndex = 0; readerIndex < numReaders; ++readerIndex) {
            final CountsReaderI reader = readers[readerIndex];
            if (needsLoading(readerIndex)) {
                if (reader.hasNextTransition()) {

                    reader.nextTransition();
                    System.out.printf("loading transition for reader[%d] position=%d length=%d count=%d %n", readerIndex, reader.getPosition(), reader.getLength(), reader.getCount());
                    final int startPosition = reader.getPosition();
                    final int endPosition = startPosition + reader.getLength();
                    startPositions[readerIndex] = startPosition;
                    endPositions[readerIndex] = endPosition;
                    startAndEndPositions.add(startPosition);
                    startAndEndPositions.add(endPosition);
                    counts[readerIndex] = reader.getCount();
                    positions[readerIndex] = startPosition;
                } else {
                    finished[readerIndex] = true;
                }
            }
        }
        first = first(startAndEndPositions);
        second = second(startAndEndPositions, first);
        length = second - first;
        if (second == Integer.MAX_VALUE) {
            length = 0;
        }
        startAndEndPositions.rem(first);
        position = first;
        hasNextTransition = length > 0;
        previousPosition = position;
        return hasNextTransition;
    }

    /**
     * Determine if we should load the next transition for a specific reader.
     *
     * @param readerIndex
     * @return
     */
    private boolean needsLoading(int readerIndex) {
        if (!finished[readerIndex]) {
            return endPositions[readerIndex] == Integer.MAX_VALUE ||
                    first + 1 >endPositions[readerIndex];
        } else {
            return false;
        }
    }

    /**
     * Returns the  minimum value in array excluding first. Returns Integer.MAX_VALUE if the array
     * is empty or contains first.
     *
     * @param array
     * @param first
     * @return
     */
    private int second(final IntArraySet array, final int first) {
        int min = Integer.MAX_VALUE;
        for (final int value : array) {
            if (value != first) min = Math.min(value, min);
        }
        return min;
    }

    private int first(final IntArraySet array) {
        int min = Integer.MAX_VALUE;
        for (final int value : array) {
            min = Math.min(value, min);
        }
        return min;
    }


    public void nextTransition() throws IOException {
        if (!hasNextTransition()) {
            throw new NoSuchElementException("no elements left in reader.");
        }
        hasNextTransition = false;
    }


    public void skipTo(int position) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public int getLength() {
        return length;
    }

    public void close() {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * Return the sum of counts over the readers that have non zero counts at the current position.
     */
    public int getCount() {
        int count = 0;
        for (int i = 0; i < numReaders; i++) {
            count += getCount(i);
        }
        return count;
    }

    public final CountsReaderI[] getReaders() {
        return readers;
    }

    /**
     * Return the count for a specific reader.
     *
     * @param readerIndex Index ((zero-based) of the reader when provided as parameter to the constructor
     * @return count for the reader identified by readerIndex.
     */
    public final int getCount(final int readerIndex) {
        return isReaderInRange(readerIndex) ? counts[readerIndex] : 0;

    }

    /**
     * Determine if the position of the reader partially overlaps with the range [first-second[
     *
     * @param readerIndex Index of the reader
     * @return True if the position of the reader overlaps
     */
    private boolean isReaderInRange(final int readerIndex) {
        final int readerStart = startPositions[readerIndex];
        final int readerEnd = endPositions[readerIndex];
        /*
      xxxx
         xxxx
       xxxxxx reader range
        <    first
          >  second
        */
        if (readerStart == Integer.MAX_VALUE) {
            return false;
        }
        boolean result = readerStart == first ||
                readerStart <= position && readerEnd > position;

        //  System.out.printf("Position=%d reader[%d]: [%d-%d[ in range=%b count=%d%n", position, readerIndex,
        //         readerStart, readerEnd, result, counts[readerIndex]);
        return result;
    }
}
