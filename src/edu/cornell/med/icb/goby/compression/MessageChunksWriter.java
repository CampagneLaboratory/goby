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

import edu.cornell.med.icb.goby.alignments.AlignmentCollectionHandler;
import edu.cornell.med.icb.goby.util.dynoptions.DynamicOptionClient;
import edu.cornell.med.icb.goby.util.dynoptions.DynamicOptionRegistry;
import edu.cornell.med.icb.goby.util.dynoptions.RegisterThis;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.*;

/**
 * Helper class to write many messages concatenated to a large output file. This helper
 * compresses each message before it is written to the output stream, and interleaves
 * messages with boundaries and size information. Boundaries make
 * it possible to split the file efficiently (e.g., see Hadoop FileSplit mechanism).
 *
 * @author Fabien Campagne
 *         Date: Apr 24, 2009
 *         Time: 4:32:35 PM
 */
public class MessageChunksWriter {
    private static final Log LOG = LogFactory.getLog(MessageChunksWriter.class);

    public static final byte DELIMITER_CONTENT = (byte) 0xFF;
    public static final int DELIMITER_LENGTH = 7;
    public static final int SIZE_OF_MESSAGE_LENGTH = 4;
    private ChunkCodec chunkCodec = null;
    /**
     * Default number of entries per chunk.
     */
    private int numEntriesPerChunk = 10000;
    private final DataOutputStream out;

    /**
     * The number of messages appended in a chunk.
     */
    private int numAppended;

    /**
     * The total number of logical entries written to the output. Multiplicity governs how
     * many logical entries are written per message.
     */
    private long totalEntriesWritten;
    private long totalBytesWritten;
    private long currentChunkStartOffset;
    private long writtenBytes = 0;
    private final boolean compressingCodec;
    @RegisterThis
    public static final DynamicOptionClient doc = new DynamicOptionClient(MessageChunksWriter.class,
            "compressing-codec:boolean, when true compress protocol buffers with new chunk codec.:false",
            "template-compression:boolean, when true use template compression.:true",
            "codec:string, name of the chunk codec to use.:gzip",
            "chunk-size:integer, the number of entries per chunk. :10000");

    public static DynamicOptionClient doc() {
        DynamicOptionRegistry.register(AlignmentCollectionHandler.doc());

        return doc;
    }

    private boolean useTemplateCompression;


    /**
     * Specify the maximum number of entries to store in any given chunk.
     *
     * @param numEntriesPerChunk maximum number of entries per chunk.
     */
    public void setNumEntriesPerChunk(final int numEntriesPerChunk) {
        this.numEntriesPerChunk = numEntriesPerChunk;
    }

    public MessageChunksWriter(final OutputStream output) {
        this.out = new DataOutputStream(output);
        compressingCodec = doc.getBoolean("compressing-codec");
        numEntriesPerChunk = doc.getInteger("chunk-size");
        final String codecName = doc.getString("codec");
        chunkCodec = ChunkCodecHelper.load(codecName);
        useTemplateCompression = doc.getBoolean("template-compression");
    }

    /**
     * Write the entry collection as needed to the output stream. When the number of entries
     * per chunk is reached, the chunk is written to disk and the collection cleared. Clients
     * can just keep adding to the collection and call writeAsNeeded for every entry.
     *
     * @param collectionBuilder The builder prepared with the growing collection of entries.
     * @throws IOException if there was an error writing the entries
     */
    public void writeAsNeeded(final com.google.protobuf.GeneratedMessage.Builder collectionBuilder)
            throws IOException {

        writeAsNeeded(collectionBuilder, 1);
    }

    /**
     * Write the entry collection as needed to the output stream. When the number of entries
     * per chunk is reached, the chunk is written to disk and the collection cleared. Clients
     * can just keep adding to the collection and call writeAsNeeded for every entry.
     *
     * @param collectionBuilder The builder prepared with the growing collection of entries.
     * @param multiplicity      Indicates how many logical entries are included in the message that
     *                          was just appended.
     * @throws IOException if there was an error writing the entries
     */
    public long writeAsNeeded(final com.google.protobuf.GeneratedMessage.Builder collectionBuilder,
                              final int multiplicity) throws IOException {
        totalEntriesWritten += Math.max(1, multiplicity);
        if (++numAppended >= numEntriesPerChunk) {
            flush(collectionBuilder);
        }
        return currentChunkStartOffset;
    }

    /**
     * Return the offset of the beginning of the current chunk (in byte, from position zero in the file).
     *
     * @return offset of the beginning of the current chunk
     */
    public long getCurrentChunkStartOffset() {
        return currentChunkStartOffset;
    }

    /**
     * Force the writing of the collection to the output stream.
     *
     * @param collectionBuilder The builder prepared with the growing collection of entries.
     * @throws IOException if there was an error writing the entries
     */
    public void flush(final com.google.protobuf.GeneratedMessage.Builder collectionBuilder)
            throws IOException {
        // Write the separation between two chunks: eight bytes with value 0xFF.

        // If we are flushing a completely empty file, that's OK, the flush() should occur.
        // Otherwise, only flush if we've appended entries.
        if (totalEntriesWritten == 0 || numAppended > 0) {
            // the position just before this chunk is written is recorded:
            currentChunkStartOffset = writtenBytes;

            assert out.size() == Integer.MAX_VALUE || out.size() == writtenBytes;

            //     System.out.println("Writting new chunk at position "+currentChunkStartOffset);
            if (LOG.isTraceEnabled()) {
                LOG.trace("writing zero bytes length=" + DELIMITER_LENGTH);
            }

            out.writeByte(chunkCodec.registrationCode());
            writtenBytes += 1;
            for (int i = 0; i < DELIMITER_LENGTH; i++) {
                out.writeByte(DELIMITER_CONTENT);
                writtenBytes += 1;
            }
            final com.google.protobuf.Message protobuffCollection = collectionBuilder.clone().build();
            // compress the read collection:

            final ByteArrayOutputStream compressedBytes = chunkCodec.encode(protobuffCollection);
            final int serializedSize = compressedBytes.size();

            if (LOG.isTraceEnabled()) {
                LOG.trace("serialized compressed size: " + serializedSize);
            }

            // write the compressed size followed by the compressed stream:
            out.writeInt(serializedSize);
            writtenBytes += 4;
            final byte[] bytes = compressedBytes.toByteArray();
            out.write(bytes);
            writtenBytes += bytes.length;
            compressedBytes.close();
            totalBytesWritten += serializedSize + 4 + DELIMITER_LENGTH;
            if (LOG.isTraceEnabled()) {
                LOG.trace("current offset: " + totalBytesWritten);

            }
            out.flush();
            numAppended = 0;
            collectionBuilder.clear();
        }
    }

    /**
     * Flush and release resources.
     *
     * @param collectionBuilder The builder prepared with the growing collection of entries.
     * @throws IOException if there is a problem closing the stream unerlying stream
     */
    public void close(final com.google.protobuf.GeneratedMessage.Builder collectionBuilder)
            throws IOException {
        flush(collectionBuilder);
        out.writeByte(chunkCodec.registrationCode());
        writtenBytes += 1;
        for (int i = 0; i < DELIMITER_LENGTH; i++) {
            out.writeByte(DELIMITER_CONTENT);
            writtenBytes += 1;
        }
        out.writeInt(0); // last collection is empty
        writtenBytes += 4;
        out.flush();
        // we do not own the output stream, so we do not close it.
    }

    /**
     * Returns the number of entries written to output.
     *
     * @return The total number of entries were written
     */
    public long getTotalEntriesWritten() {
        return totalEntriesWritten;
    }

    /**
     * Returns the number of bytes written to output.
     *
     * @return The total number of bytes that have been written
     */
    public long getTotalBytesWritten() {
        return totalBytesWritten;
    }

    /**
     * Print statistics.
     *
     * @param writer the writer used to print the statistics
     */
    public void printStats(final PrintWriter writer) {
        writer.println("Total logical entries written: " + totalEntriesWritten);
        writer.println("Total bytes written: " + totalBytesWritten);
        writer.println("Average bytes/logical entry: "
                + (float) totalBytesWritten / (float) totalEntriesWritten);
        writer.flush();
    }

    /**
     * Print statistics.
     *
     * @param out Where to print.
     */
    public void printStats(final PrintStream out) {
        printStats(new PrintWriter(out));
    }

    /**
     * The number of entries appended in the current chunk. Zero indicates the start of a new chunk.
     *
     * @return
     */
    public int getAppendedInChunk() {
        return numAppended;
    }

    public void setParser(final ProtobuffCollectionHandler protobuffCollectionHandler) {
        if (chunkCodec == null) {
            if (protobuffCollectionHandler instanceof AlignmentCollectionHandler) {

                chunkCodec = compressingCodec ? new HybridChunkCodec1() : new GZipChunkCodec();
            } else {
                chunkCodec = new GZipChunkCodec();
            }
        }
        protobuffCollectionHandler.setUseTemplateCompression(useTemplateCompression);
        chunkCodec.setHandler(protobuffCollectionHandler);

        //     chunkCodec = new GZipChunkCodec();

    }


}
