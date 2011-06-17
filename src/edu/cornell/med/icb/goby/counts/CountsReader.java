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

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.io.BinIO;
import it.unimi.dsi.fastutil.io.RepositionableStream;
import it.unimi.dsi.io.InputBitStream;
import org.bdval.io.compound.CompoundDataInput;

import java.io.DataInput;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;

/**
 * Reads counts encoded written by {@link CountsWriter}.
 *
 * @author Fabien Campagne
 *         Date: May 6, 2009
 *         Time: 3:25:30 PM
 */
public class CountsReader implements CountsReaderI {
    private final InputBitStream input;
    protected static final int END_OF_DATA_MARKER = 277492431;
    private boolean endOfStream;
    private int deltaCount;
    private int currentCount;
    private int length = -1;
    final int[] positions;
    final int[] offsets;
    final int[] counts;
    /**
     * True iff {@link #hasNextTransition()} will fetch the next transition or will find
     * the end of file.
     */
    private boolean nextTransitionLoaded;
    private int count;

    /**
     * Position on the sequence before the count transition is observed.
     */
    private int position = -1;
    /**
     * This boolean is true if the reader has an index.
     */
    private boolean hasIndex;


    /**
     * Return the position along the sequence where the count is observed.
     *
     * @return
     */
    public int getPosition() {
        return position;
    }

    public CountsReader(final InputStream inputStream) throws IOException {
        input = new InputBitStream(inputStream);
        currentCount = input.readDelta() - 1;
        count = currentCount;
        counts = offsets = positions = null;

    }

    public CountsReader(final InputBitStream inputBitStream) throws IOException {
        input = inputBitStream;
        currentCount = input.readDelta() - 1;
        count = currentCount;
        counts = offsets = positions = null;
    }


    public CountsReader(InputStream inputStream, DataInput indexInputStream) throws IOException {
        assert inputStream instanceof RepositionableStream : "inputStream must be repositionable.";
        input = new InputBitStream(inputStream);
        currentCount = input.readDelta() - 1;
        count = currentCount;
        if (indexInputStream != null) {
            final int length = indexInputStream.readInt();

            positions = new int[length];
            BinIO.loadInts(indexInputStream, positions);
            offsets = new int[length];
            BinIO.loadInts(indexInputStream, offsets);
            counts = new int[length];
            BinIO.loadInts(indexInputStream, counts);
            hasIndex = true;
        } else {
            counts = offsets = positions = null;
        }
    }

    /**
     * Determines if the reader has data about another transition.
     *
     * @return True when a call to nextTransition() will succeed, False otherwise.
     * @throws IOException
     */
    public boolean hasNextTransition() throws IOException {
        if (nextTransitionLoaded) {
            return true;
        }
        if (endOfStream) {
            return false;
        }
        final int deltaCount = input.readGamma();


        if (deltaCount == END_OF_DATA_MARKER) {
            endOfStream = true;
            return false;
        }
        position += Math.max(1, length);
        length = input.readGamma();

        final int decodedDeltaCount = decodeDeltaCount(deltaCount);
        this.deltaCount = decodedDeltaCount;
        count += decodedDeltaCount;
        assert count>=0:"Count must never be negative! now at position ="+position;

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

    public int getDeltaCount() {
        return deltaCount;
    }

    /**
     * Check if count information is available for another position.
     *
     * @return
     * @throws IOException
     */
    public boolean hasNextPosition() throws IOException {
        if (length > 0) {
            return true;
        }
        if (endOfStream) {
            return false;
        }
        final int deltaCount = input.readGamma();

        if (deltaCount == END_OF_DATA_MARKER) {
            endOfStream = true;
            return false;
        }

        length = input.readGamma();

        final int decodedDeltaCount = decodeDeltaCount(deltaCount);
        currentCount += decodedDeltaCount;
        return true;
    }

    /**
     * Returns count for the next position.
     *
     * @return
     * @throws IOException
     */
    public int nextCountAtPosition() throws IOException {
        if (hasNextPosition()) {
            --length;
            ++position;
            return currentCount;
        } else {
            throw new IllegalStateException("next cannot be called when hasNext would return false.");
        }
    }

    protected static int decodeDeltaCount(final int deltaCount) {
        if (deltaCount % 2 == 1) {
            // odd encoded value, recover the negative value:
            return -((deltaCount - 1) / 2);
        } else {
            return deltaCount / 2;
        }
    }

    public final int getLength() {
        return length;
    }

    public final int getCount() {
        return count;
    }

    /**
     * {@inheritDoc}
     */
    public void close() throws IOException {
        input.close();
    }

    /**
     * Advance up to or past the specified position. The reader is advanced until the position returned by getPosition()
     * is at least equal, or greater to the specified position.
     *
     * @param position
     * @throws IOException
     */
    public void skipTo(final int position) throws IOException {
        if (position < this.position) {
            if (hasNextTransition()) {
                nextTransition();
            }
            return;
        }
        if (hasIndex) {
            reposition(position);
        } else {
            // skip to the specified position

            while (hasNextTransition()) {
                nextTransition();
                if (getPosition() >= position) {
                    break;
                }
            }
        }

    }
    /**
     * Reposition the reader on a genomic position. In contrast to skipTo, this method reposition to any position,
     * including position that occured before the position returned last by getPosition().
     *
     * @param position Position to seek to.
     * @throws IOException If an error occurs accessing the index or counts data.
     */
    public void reposition(final int position) throws IOException {
        if (!hasIndex) {
            throw new IllegalStateException("The Counts must have an index to use the reposition method.");
        }
        final int r = Arrays.binarySearch(positions, position);
        // position at the position if found, or immediately before if the position was not found in the index
        //   System.out.println("CountsReader.reposition " + position);
        final int ip = r >= 0 ? r : -(r + 1);
        final int index = r >= 0 ? r : Math.max(0, ip);
        int priorIndex = index - 1;


        if (priorIndex < 0 || index >= positions.length) {
            // the index does not contain the position, go back to the beginning of the file.
            this.position = position;
            this.count = 0;
            deltaCount = 0;
            nextTransitionLoaded = false;
            input.position(0);
            endOfStream = false;
            return;
            //    throw new IllegalStateException(String.format("Positions array (length=%d) is too small for index %d.",
            //          positions.length, index));
        }
        if (positions[index] == position) {
            this.position = position;
            this.count = counts[index];
            deltaCount = 0;
            nextTransitionLoaded = false;

            input.position(offsets[index]);
            endOfStream = false;
            input.readGamma(); // advance delta count, ignore the value, we know the count already from the index.
            length = input.readGamma();
            return;

        } else {
            // initialize read data structure at the position of the index
            this.position = positions[priorIndex];
            count = priorIndex < 0 ? 0 : counts[priorIndex];
            input.position(offsets[priorIndex]);
            endOfStream = false;
            nextTransitionLoaded = false;
            currentCount = count;
            deltaCount = 0;
            length = 0;
            input.readGamma(); // advance delta count, ignore the value, we know the count already from the index.
            length = input.readGamma();


            // now iterate until we meet the skipTo position condition:
            while (hasNextTransition()) {
                nextTransition();
                //   System.out.printf("CountsReader.reposition seeking %d currentPosition=%d %n", position,  getPosition());
                if (getPosition() >= position) {
                    break;
                }
            }
        }
    }

    public boolean isPositionInIndex(int i) {
        return Arrays.binarySearch(positions, i)>=0;
    }
}
