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
import org.apache.commons.io.IOUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.ByteArrayInputStream;
import java.io.Closeable;
import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.zip.GZIPInputStream;

/**
 * Reads from a stream produced with {@link edu.cornell.med.icb.goby.reads.MessageChunksWriter}.
 *
 * @author Fabien Campagne
 *         Date: Apr 30, 2009
 *         Time: 5:06:55 PM
 */
public class MessageChunksReader implements Closeable {
    private static final Log LOG = LogFactory.getLog(MessageChunksReader.class);
    protected DataInputStream in;
    protected int entryIndex;
    protected GZIPInputStream uncompressStream;

    protected MessageChunksReader() {
        super();
    }

    public MessageChunksReader(final InputStream input) {
        super();
        assert input != null : "The input stream must not be null";
        in = new DataInputStream(input);
    }

    /**
     * Returns true if the input has more entries.
     *
     * @param collection     The current collection, or null if no collection has been read yet.
     * @param collectionSize The size of the current collection (can be zero).
     * @return True if the input has more entries, False otherwise.
     */
    public boolean hasNext(final GeneratedMessage collection,
                           final int collectionSize) {
        if (collection == null || entryIndex >= collectionSize) {
            if (in == null) {
                return false;
            }

            try {
                final long numSkipped = in.skip(MessageChunksWriter.DELIMITER_LENGTH);
                if (numSkipped != MessageChunksWriter.DELIMITER_LENGTH) {
                    LOG.warn("Skip returned " + numSkipped);
                    uncompressStream = null;
                    return false;
                }

                // read the number of compressed bytes to follow:
                final int numBytes = in.readInt();
                if (numBytes == 0) {
                    uncompressStream = null;
                    return false;
                }

                // read the compressed stream:
                final byte[] bytes = new byte[numBytes];
                in.read(bytes);
                // uncompress the bytes on the fly and read into collection:
                uncompressStream = new GZIPInputStream(new ByteArrayInputStream(bytes));
                // read the 8 bytes of separation:

                entryIndex = 0;
                return true;
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        } else {
            uncompressStream = null;
        }

        return entryIndex < collectionSize;
    }

    /**
     * Get the compressed representation of the next available collection. Is null when the
     * current collection is not expired.
     *
     * @return a compressed stream from where to deserialize the next collection.
     */
    public GZIPInputStream getUncompressStream() {
        return uncompressStream;
    }

    /**
     * Returns the current entry index and increment.
     *
     * @return The current entry index
     */
    public int incrementEntryIndex() {
        return entryIndex++;
    }

    /**
     * {@inheritDoc}
     */
    public void close() {
        IOUtils.closeQuietly(uncompressStream);
        IOUtils.closeQuietly(in);
    }
}
