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
import edu.cornell.med.icb.goby.exception.GobyRuntimeException;
import it.unimi.dsi.fastutil.bytes.ByteSet;
import it.unimi.dsi.fastutil.io.FastBufferedInputStream;
import org.apache.commons.io.IOUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

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
    private static final Log LOG = LogFactory.getLog(FastBufferedMessageChunksReader.class);
    private final long end;
    private final FastBufferedInputStream input;

    /**
     * Start offset of the slice in the file, in bytes.
     */
    private final long startOffset;
    /**
     * End offset of the slice in the file, in bytes.
     */
    private final long endOffset;


    private ByteSet supportedCodecRegistrationCodes;
    private boolean withinSlice = true;

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
        startOffset = start;
        endOffset = end;
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
        this.in = new DataInputStream(input);
        supportedCodecRegistrationCodes = ChunkCodecHelper.registrationCodes();
        reposition(start, end);
    }

    byte lastCodecCodeSeen = 0;

    private void reposition(final long start, final long end) throws IOException {
        assert end >= start : "end must be larger than start ";
        if (start == 0 && input.position() == 0) {

            return;
        }

      /*  if (start > input.length()) {
            withinSlice = false;
            return;
        }*/
        input.position(start);
        if (input.position()!=start) {
          // must have happened because we are past the end of stream.
            withinSlice = false;
            return;
        }
        in = new DataInputStream(input);
        int contiguousDelimiterBytes = 0;
        long skipped = 0;
        long position = 0;

        withinSlice = true;
        boolean codecSeen = false;
        // search though the input stream until a delimiter chunk, end of stream, or end of slice is reached
        int b;
        while ((b=in.read())!=-1 && input.position() < end) {

            final byte c = (byte)b;
            //     System.out.printf("%2X(%d) ", c, contiguousDelimiterBytes);
            //    System.out.flush();
            if (!codecSeen && hasValidCodecCode(c)) {
                skipped++;
                contiguousDelimiterBytes++;
                codecSeen = true;
                continue;
            }
            if (codecSeen && c == MessageChunksWriter.DELIMITER_CONTENT) {
                contiguousDelimiterBytes++;

            } else {
                if (contiguousDelimiterBytes == MessageChunksWriter.DELIMITER_LENGTH + 1) {
                    chunkCodec = ChunkCodecHelper.withRegistrationCode(lastCodecCodeSeen);
                    // position exactly after the 7th 0xFF byte, past the first byte of the size:
                    // the first byte of size was already read and is provided in c.
                    long positionBeforeValidation = input.position();
                    if (!chunkCodec.validate(c, in)) {
                        LOG.warn(String.format("Found spurious boundary around position %d ", input.position()));
                        contiguousDelimiterBytes = 0;
                        chunkCodec = null;
                        lastCodecCodeSeen = 0;
                        codecSeen = true;
                        continue;
                    }

                    final long newPosition = positionBeforeValidation - (MessageChunksWriter.DELIMITER_LENGTH + 2);
                    input.position(newPosition);
                    return;

                }
                if (skipped > MessageChunksWriter.DELIMITER_LENGTH + 1) {
                    contiguousDelimiterBytes = 0;
                    codecSeen = false;
                }


            }
            ++skipped;
        }
        if (b==-1||input.position() >= end) {
            withinSlice = false;
        }
        position = start + skipped;
        streamPositionAtStart = input.position();


    }


    private boolean hasValidCodecCode(byte registrationCode) throws IOException {

        if (supportedCodecRegistrationCodes.contains(registrationCode)) {
            lastCodecCodeSeen = registrationCode;
            return true;
        } else return false;


    }


    /**
     * Seek to the given position in the compact file.
     *
     * @param position Position where to seek to.
     * @throws IOException If an error occurs reading this file.
     */
    public void seek(final long position) throws IOException {
        input.flush();
        reposition(position, Long.MAX_VALUE);
        // invalidate any bytes already read:
        compressedBytes = null;

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
            return withinSlice && in != null && super.hasNext(collection, collectionSize);
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


    /**
     * Flush buffer so that content will be read from the input again.
     */
    public void flush() {

        input.flush();
    }
}
