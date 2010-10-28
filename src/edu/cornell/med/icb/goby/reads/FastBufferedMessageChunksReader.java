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

package edu.cornell.med.icb.goby.reads;

import com.google.protobuf.GeneratedMessage;
import edu.cornell.med.icb.goby.exception.GobyRuntimeException;
import it.unimi.dsi.fastutil.io.FastBufferedInputStream;
import org.apache.commons.io.IOUtils;

import java.io.DataInputStream;
import java.io.IOException;

/**
 * Reads from a stream produced with {@link edu.cornell.med.icb.goby.reads.MessageChunksWriter}.
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
     * @param end   The end index for the split
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
        if (end != Long.MAX_VALUE && start != end) {
            end += MessageChunksWriter.DELIMITER_LENGTH + 4;
        }
        this.end = end;
        this.input = input;
        reposition(start, end);
    }

    private void reposition(final long start, final long end) throws IOException {
        input.position(start);


        int b;
        int contiguousZeroBytes = 0;
        long skipped = 0;
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
                in = new DataInputStream(input);
                final long seekPosition = start + skipped - contiguousZeroBytes;
                input.position(seekPosition);
                break;
            }
            position = start + skipped;
        }
    }

    /**
     * Seek to the given position in the compact file.
     *
     * @param position Position where to seek to.
     * @throws IOException If an error occurs reading this file.
     */
    public void seek(final long position) throws IOException {
        reposition(position, Long.MAX_VALUE);
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
                    throw new GobyRuntimeException(e);
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

    public long position() throws IOException {
        return input.position();
    }
}
