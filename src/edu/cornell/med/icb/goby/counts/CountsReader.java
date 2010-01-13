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

package edu.cornell.med.icb.goby.counts;

import it.unimi.dsi.io.InputBitStream;

import java.io.IOException;
import java.io.InputStream;

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
     * Return the position along the sequence where the count is observed.
     * @return
     */
    public int getPosition() {
        return position;
    }

    public CountsReader(final InputStream inputStream) throws IOException {
        input = new InputBitStream(inputStream);
        currentCount = input.readDelta() - 1;
        count = currentCount;
    }

    /**
     * Determines if the reader has data about another transition.
     * @return True when a call to nextTransition() will succeed, False otherwise.
     *
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
     * @param position
     * @throws IOException
     */
    public void skipTo(final int position) throws IOException {
        // skip to the specified position
        while (hasNextTransition()) {
            nextTransition();
            if (getPosition() >= position) {
                break;
            }
        }
    }
}
