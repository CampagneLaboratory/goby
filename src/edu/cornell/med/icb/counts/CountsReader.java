/*
 * Copyright (C) 2009 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 * WEILL MEDICAL COLLEGE OF CORNELL UNIVERSITY MAKES NO REPRESENTATIONS
 * ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY PURPOSE. IT IS PROVIDED
 * "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY. THE WEILL MEDICAL COLLEGE
 * OF CORNELL UNIVERSITY SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY
 * THE USERS OF THIS SOFTWARE.
 */

package edu.cornell.med.icb.counts;

import it.unimi.dsi.io.InputBitStream;

import java.io.Closeable;
import java.io.IOException;
import java.io.InputStream;

/**
 * Reads counts encoded written by CountsWrite.
 *
 * @author Fabien Campagne
 *         Date: May 6, 2009
 *         Time: 3:25:30 PM
 */
public class CountsReader implements Closeable, CountsReaderI {
    private InputBitStream input;
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
    public void skipTo(int position) throws IOException {
        // skip to the specified position
        while (hasNextTransition()) {
            nextTransition();
            if (getPosition() >= position) break;
        }
    }
}
