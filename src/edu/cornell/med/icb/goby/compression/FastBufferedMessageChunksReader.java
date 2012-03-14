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

package edu.cornell.med.icb.goby.compression;

import com.google.protobuf.GeneratedMessage;
import edu.cornell.med.icb.goby.compression.MessageChunksReader;
import edu.cornell.med.icb.goby.compression.MessageChunksWriter;
import edu.cornell.med.icb.goby.exception.GobyRuntimeException;
import it.unimi.dsi.fastutil.io.FastBufferedInputStream;
import org.apache.commons.io.IOUtils;

import java.io.DataInputStream;
import java.io.IOException;

/**
 * Reads from a stream produced with {@link edu.cornell.med.icb.goby.compression.MessageChunksWriter}.
 *
 * @author Fabien Campagne
 *         Date: Apr 30, 2009
 *         Time: 5:06:55 PM
 */
public class FastBufferedMessageChunksReader extends MessageChunksReader {
    private final long end;
    private final FastBufferedInputStream input;
    /**
     * Start offset of the slice in the file, in bytes.
     */
    private long startOffset;
    /**
     * End offset of the slice in the file, in bytes.
     */
    private long endOffset;

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
        this.startOffset = start;
        this.endOffset = end;
        if (start < 0L) {
            throw new IllegalArgumentException("Start position ("
                    + start + ") must not be less than zero");
        }
        if (end != Long.MAX_VALUE && end < 0L) {
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
        assert end >= start : "end must be larger than start ";
        input.position(start);

        int b;
        int contiguousDelimiterBytes = 0;
        long skipped = 0;
        long position = 0;

        // search though the input stream until a delimiter chunk or end of stream is reached
        while ((b = input.read()) != -1) {
            final byte c = (byte) b;

            if (c == MessageChunksWriter.DELIMITER_CONTENT) {
                contiguousDelimiterBytes++;
            } else {
                contiguousDelimiterBytes = 0;
            }
            ++skipped;
            if (contiguousDelimiterBytes == MessageChunksWriter.DELIMITER_LENGTH) {
                if (hasFF(input)) {
                    skipped++;
                }
                if (skipped >= MessageChunksWriter.DELIMITER_LENGTH + 1) {
                    // make sure we have seen the delimited AND the codec registration code since start, otherwise continue looking
                    // a delimiter was found, start reading data from here
                    in = new DataInputStream(input);
                    final long seekPosition = start + skipped - MessageChunksWriter.DELIMITER_LENGTH - 1; // positions  before the codec registration code.
                    input.position(seekPosition);
                    break;
                }
            }
            position = start + skipped;
        }
    }

    private boolean hasFF(FastBufferedInputStream input) throws IOException {
        if (input.available() >= 1) {
            int b = input.read();
            byte code = (byte) b;
            if (code == -1) return true;
        }
        return false;
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
            compressedBytes = null;
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
