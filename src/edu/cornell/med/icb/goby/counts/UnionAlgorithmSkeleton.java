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


import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntArraySet;

import java.io.IOException;
import java.util.NoSuchElementException;

/**
 * Skeleton for algorithm development discussion.
 *
 * @author Fabien Campagne
 *         Date: 5/26/11
 *         Time: 10:14 PM
 */
public class UnionAlgorithmSkeleton implements CountsReaderI {
    private int numReaders;
    private CountsReaderI[] readers;
    private boolean hasNextTransition;
    private int length;
    private int position = 0;
    private int first;
    private int second;
    private int[] positions;

    public UnionAlgorithmSkeleton(CountsReaderI... readers) {
        this.numReaders = readers.length;
        this.readers = readers;
        this.counts = new int[this.numReaders];
        this.positions = new int[this.numReaders];
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

        for (int readerIndex = 0; readerIndex < numReaders; ++readerIndex) {
            final CountsReaderI reader = readers[readerIndex];
            if (reader.hasNextTransition()) {
                reader.nextTransition();
                startAndEndPositions.add(reader.getPosition());
                startAndEndPositions.add(reader.getPosition() + reader.getLength());
                counts[readerIndex] = reader.getCount();
                positions[readerIndex] = reader.getPosition();
            }
        }

        first = first(startAndEndPositions);
        second = second(startAndEndPositions, first);
        //TODO handle if Integer.MAX_VALUE is back.
        length = second - first;
        startAndEndPositions.rem(first);
        position = first;
        hasNextTransition=startAndEndPositions.size() >= 2;
        return hasNextTransition;
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
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
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
        return isReaderInRange(readerIndex) ? counts[readerIndex]:0;

    }

    private boolean isReaderInRange(final int readerIndex) {
        final int readerPosition = positions[readerIndex];
        if (readerPosition >=first  && readerPosition<second) return true;
        else return false;
    }
}
