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

package edu.cornell.med.icb.goby.counts;

import it.unimi.dsi.fastutil.ints.IntAVLTreeSet;
import it.unimi.dsi.fastutil.ints.IntSortedSet;

import java.io.IOException;
import java.util.Arrays;

/**
 * Iterates through a set of count readers, returning a transition whenever one of the
 * underlying readers would report a transition at the position. A transition is never
 * triggered at positions where the underlying readers are all constant, making it
 * efficient to compare counts across readers.
 *
 * @author Fabien Campagne
 *         Date: Jun 13, 2009
 *         Time: 2:00:21 PM
 */
public class AnyTransitionCountsIterator implements CountsReaderI {

    private final CountsReaderI[] readers;
    /**
     * The position each reader is at:
     */
    private final int[] position;
    /**
     * The count the reader sees at the current position.
     */
    private final int[] count;
    /**
     * The length each reader has, at its own position.
     */
    private int[] length;
    /**
     * Contains the start and end position of the last entry obtained from each active reader. An extremity is
     * removed when a reader is no longer in scope.
     */
    private IntSortedSet extremities = new IntAVLTreeSet();
    private boolean nextTransitionLoaded;
    private final int numReaders;
    private int currentPosition;
    private int currentLength;
    /**
     * An array whose element is true when the reader is completely exhausted.
     */
    private boolean[] finished;
    /**
     * For readers that are exhausted, the following array records the last position the reader returned.
     */
    private int[] lastPosition;


    public AnyTransitionCountsIterator(final CountsReaderI... countReader) {
        super();
        readers = countReader;
        numReaders = countReader.length;
        position = new int[numReaders];
        Arrays.fill(position, -1);
        count = new int[numReaders];
        currentPosition = 0;
        finished = new boolean[numReaders];
        lastPosition = new int[numReaders];
        length = new int[numReaders];
    }

    public boolean hasNextTransition() throws IOException {
        if (nextTransitionLoaded) {
            return true;
        }
       
//        extremities.clear();
        // advance by at least one position:
        //   currentPosition++;
        // load count for each reader that transitions at this position:
        int countReadersFinished = 0;
        for (int i = 0; i < numReaders; i++) {


            // passed the previous transition, must advance this reader.
            if (!readers[i].hasNextTransition()) {
                if (!finished[i]) {
                    lastPosition[i] = position[i];
                    //       System.out.printf("reader %d finished at position %d %n", i, position[i]);
                }
                // reader has no more transitions.

                finished[i] = true;
                countReadersFinished++;
                position[i]++;
                if (currentPosition < lastPosition[i]) {
                    // last position is still in scope:
                    extremities.add(lastPosition[i]);
                }
            } else {
                // load the next transition:
                if (currentPosition >= position[i]) {
                    readers[i].nextTransition();
                    position[i] = readers[i].getPosition();
                    length[i] = readers[i].getLength();
                    extremities.add(position[i]);
                    extremities.add(position[i] + length[i] - 1);
                }
            }
        }

        // determine the new current position, since we may have advanced all the readers in a big jump:

        // length is current position minus previous position.
        if (!extremities.isEmpty()) {
            currentLength = extremities.firstInt() - currentPosition; 
            currentPosition = extremities.firstInt();
        }   else {
            currentPosition++;
            currentLength=0;
        }
        for (int i = 0; i < numReaders; i++) {
            // set the reader count if current position is within the window of the reader's transition:

            final int readerPosition = position[i];
            if (currentPosition >= (readerPosition - length[i]) && currentPosition <= readerPosition) {
                count[i] = readers[i].getCount();
            } else {
                count[i] = 0;
            }
     //       System.out.printf("reader[%d] start=%d end=%d position=%d count=%d%n", i, readerPosition - length[i], readerPosition, currentPosition, count[i]);

            if (currentPosition > readerPosition) {
                extremities.rem(readerPosition);
                extremities.rem(readerPosition + length[i] - 1);
            }
        }
        if (countReadersFinished == numReaders) {
            // no more transitions in any reader.
            return false;
        }
        //    System.out.println("currentPosition=" + currentPosition);
        extremities.rem(currentPosition);
        nextTransitionLoaded = true;
        return true;
    }



    /**
     * Advance to the next transition. After this method has been called successfully, position,
     * length, deltaCount and currentCount are available through getters of this reader.
     *
     * @throws IOException
     */
    public void nextTransition() throws IOException {
        if (hasNextTransition()) {
            nextTransitionLoaded = false;

        } else {
            throw new IllegalStateException("next cannot be called when hasNext would return false.");
        }
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

    public final int getCount(final int readerIndex) {
        return count[readerIndex];
    }

    public void close() throws IOException {
        for (final CountsReaderI reader : readers) {
            reader.close();
        }
    }

    public void skipTo(final int position) throws IOException {
        throw new UnsupportedOperationException("skipTo is not currently supported by this implementation.");
    }

    public int getLength() {
        // System.out.println("currentLength=" + currentLength);
        return currentLength;
    }

    /**
     * Position on the sequence before the count transition is observed.
     *
     * @return
     */
    public int getPosition() {

        final int value = currentPosition;
        //      System.out.println("currentPosition=" + value);
        return value;
    }
}
