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

package edu.cornell.med.icb.reads;

import com.google.protobuf.GeneratedMessage;
import it.unimi.dsi.fastutil.io.FastBufferedInputStream;
import org.apache.commons.io.IOUtils;

import java.io.DataInputStream;
import java.io.IOException;

/**
 * Reads from a stream produced with {@link edu.cornell.med.icb.reads.MessageChunksWriter}.
 *
 * @author Fabien Campagne
 *         Date: Apr 30, 2009
 *         Time: 5:06:55 PM
 */
public class FastBufferedMessageChunksReader extends MessageChunksReader {
    private final long end;
    private FastBufferedInputStream input;

    /**
     * Support for splitting the entries on the file system.
     * Seek the input to start and start looking for the beginning of a new collection.
     * When found, return all entries in the collection through hasNext, next().
     * Will possibly return additional collection of entries, but will stop
     * returning new entries if the position in the input stream is past end.
     *
     * @param start The start index for the split
     * @param end The end index for the split
     * @param input The input stream containing the data
     * @throws IOException if there is a problem reading from the stream
     */
    public FastBufferedMessageChunksReader(final long start, long end,
                                           final FastBufferedInputStream input) throws IOException {
        super();

        if (start < 0L) {
            throw new IllegalArgumentException("Start position ("
                    + start + ") must not be less than zero");
        }
        if (end < 0L) {
            throw new IllegalArgumentException("End position ("
                    + end + ") must not be less than zero");
        }
        if (start > end) {
            throw new IllegalArgumentException("Start position ("
                    + start + ") must not be greater than the end position (" + end + ")");
        }

        input.position(start);
        this.end = end;
        if (end != Long.MAX_VALUE) {
            end += MessageChunksWriter.DELIMITER_LENGTH + 4;
        }

        int b;
        int contiguousZeroBytes = 0;
        int skipped = 0;
        long position = 0;

        // search though the input stream until a delimiter chunk or end of stream is reached
        while (position < end && (b = input.read()) != -1) {
            if (b == MessageChunksWriter.DELIMITER_CONTENT) {
                contiguousZeroBytes++;
            } else {
                contiguousZeroBytes = 0;
            }
            ++skipped;
            if (contiguousZeroBytes == MessageChunksWriter.DELIMITER_LENGTH) {
                // a delimiter was found, start reading data from here
                this.input = input;
                in = new DataInputStream(input);
                final long seekPosition = start + skipped - contiguousZeroBytes;
                input.position(seekPosition);
                break;
            }
            position = start + skipped;
        }
    }

    /**
     * Returns true if the input has more entries.
     *
     * @param collection     The current collection, or null if no collection has been read yet.
     * @param collectionSize The size of the current collection (can be zero).
     * @return True if the input has more entries, False otherwise.
     */
    @Override
    public boolean hasNext(final GeneratedMessage collection, final int collectionSize) {
        // do not read a new collection if we are past the end of the file split allocated to us
        if (collection == null || entryIndex >= collectionSize) {
            if (input != null) {
                try {
                    if (input.position() >= end) {
                        return false;
                    }
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }
            return in != null && super.hasNext(collection, collectionSize);
        } else {
            uncompressStream = null;
        }

        return entryIndex < collectionSize;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() {
        super.close();
        IOUtils.closeQuietly(input);
    }
}
