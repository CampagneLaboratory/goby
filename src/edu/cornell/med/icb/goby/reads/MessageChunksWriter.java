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

package edu.cornell.med.icb.goby.reads;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.ByteArrayOutputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.zip.GZIPOutputStream;

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

    protected static final int DELIMITER_CONTENT = 0xFF;
    protected static final int DELIMITER_LENGTH = 8;

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
    }

    /**
     * Write the entry collection as needed to the output stream. When the number of entries
     * per chunk is reached, the chunk is written to disk and the collection cleared. Clients
     * can just keep adding to the collection and call writeAsNeeded for every entry.
     *
     * @param collectionBuilder The builder prepared with the growing collection of entries.
     * @throws IOException
     */
    public void writeAsNeeded(final com.google.protobuf.GeneratedMessage.Builder collectionBuilder) throws IOException {
        ++totalEntriesWritten;
        if (++numAppended >= numEntriesPerChunk) {
            flush(collectionBuilder);
        }
    }

    /**
     * Write the entry collection as needed to the output stream. When the number of entries
     * per chunk is reached, the chunk is written to disk and the collection cleared. Clients
     * can just keep adding to the collection and call writeAsNeeded for every entry.
     *
     * @param collectionBuilder The builder prepared with the growing collection of entries.
     * @param multiplicity Indicates how many logical entries are included in the message that
     * was just appended.
     * @throws IOException
     */
    public void writeAsNeeded(final com.google.protobuf.GeneratedMessage.Builder collectionBuilder,
                              final int multiplicity) throws IOException {
            totalEntriesWritten += multiplicity;
            if (++numAppended >= numEntriesPerChunk) {
                flush(collectionBuilder);
            }

        }

    /**
     * Force the writting of the collection to the output stream.
     *
     * @param collectionBuilder The builder prepared with the growing collection of entries.
     * @throws IOException If an error occurs writting to the output stream.
     */
    public void flush(final com.google.protobuf.GeneratedMessage.Builder collectionBuilder) throws IOException {
        // Write the separation between two chunks: eight bytes with value zero.
        if (LOG.isTraceEnabled()) {
            LOG.trace("writting zero bytes {" + DELIMITER_LENGTH);
        }
        for (int i = 0; i < DELIMITER_LENGTH; i++) {
            out.writeByte(DELIMITER_CONTENT);
        }
        final com.google.protobuf.Message readCollection = collectionBuilder.build();

        // compress the read collection:
        final ByteArrayOutputStream compressedStream = new ByteArrayOutputStream();
        final OutputStream byteArrayOutputStream = new GZIPOutputStream(compressedStream);
        readCollection.writeTo(byteArrayOutputStream);
        byteArrayOutputStream.close();

        final int serializedSize = compressedStream.size();
        if (LOG.isTraceEnabled()) {
            LOG.trace("serializedSize: " + serializedSize);
        }

        // write the compressed size followed by the compressed stream:
        out.writeInt(serializedSize);
        out.write(compressedStream.toByteArray());

        totalBytesWritten += serializedSize + 4 + DELIMITER_LENGTH;

        out.flush();
        numAppended = 0;
        collectionBuilder.clear();
    }

    /**
     * Flush and release resources.
     *
     * @param collectionBuilder The builder prepared with the growing collection of entries.
     * @throws IOException
     */
    public void close(final com.google.protobuf.GeneratedMessage.Builder collectionBuilder) throws IOException {
        flush(collectionBuilder);
        for (int i = 0; i < DELIMITER_LENGTH; i++) {
            out.writeByte(0); //last collection is empty
        }
        out.writeInt(0); //last collection is empty
        out.flush();
        // we do not own the output stream, so we do not close it.
    }

    /**
     * Returns the number of entries written to output.
     *
     * @return
     */
    public long getTotalEntriesWritten() {
        return totalEntriesWritten;
    }

    /**
     * Returns the number of bytes written to output.
     *
     * @return
     */
    public long getTotalBytesWritten() {
        return totalBytesWritten;
    }

    /**
     * Print statistics.
     *
     * @param out Where to print.
     */
    public void printStats(final PrintWriter out) {
        out.println("Total logical entries written: " + totalEntriesWritten);
        out.println("Total bytes written: " + totalBytesWritten);
        out.println("Average bytes/logical entry: " + (float) totalBytesWritten / (float) totalEntriesWritten);
        out.flush();
    }

    /**
     * Print statistics.
     *
     * @param out Where to print.
     */
    public void printStats(final PrintStream out) {
        printStats(new PrintWriter(out));
    }
}
