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

import com.google.protobuf.ByteString;

import java.io.Closeable;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

/**
 * Write reads to the compact format.
 *
 * @author Fabien Campagne
 *         Date: Apr 24, 2009
 *         Time: 4:32:35 PM
 */
public class ReadsWriter implements Closeable {
    private final Reads.ReadCollection.Builder collectionBuilder;

    private CharSequence sequence;
    private CharSequence description;
    private CharSequence identifier;
    private byte[] qualityScores;

    private int readIndex;
    private int previousReadLength;
    private long sequenceBasesWritten;
    private final MessageChunksWriter messageChunkWriter;

    private byte[] byteBuffer = new byte[100];

    public ReadsWriter(final OutputStream output) {
        collectionBuilder = Reads.ReadCollection.newBuilder();
        messageChunkWriter = new MessageChunksWriter(output);
        readIndex = 0;
    }

    public synchronized void setQualityScores(final byte[] qualityScores) {
        this.qualityScores = qualityScores;
    }

    public synchronized void setDescription(final CharSequence description) {
        this.description = description;
    }

    public synchronized void setSequence(final CharSequence sequence) {
        this.sequence = sequence;
    }

    public synchronized void appendEntry(final CharSequence description,
                                         final CharSequence sequence,
                                         final byte[] qualityScores) throws IOException {
        this.description = description;
        this.sequence = sequence;
        this.qualityScores = qualityScores;
        appendEntry();
    }

    public synchronized void appendEntry(final CharSequence description,
                                         final CharSequence sequence) throws IOException {
        this.description = description;
        this.sequence = sequence;
        appendEntry();
    }

    public synchronized void appendEntry(final CharSequence sequence) throws IOException {
        this.sequence = sequence;
        appendEntry();
    }

    /**
     * {@inheritDoc}
     */
    public void close() throws IOException {
        messageChunkWriter.close(collectionBuilder);
    }

    public synchronized void appendEntry() throws IOException {
        final Reads.ReadEntry.Builder entryBuilder = Reads.ReadEntry.newBuilder();

        entryBuilder.setReadIndex(readIndex++);
        if (description != null) {
            entryBuilder.setDescription(description.toString());
            description = null;
        }
        if (identifier != null) {
            entryBuilder.setReadIdentifier(identifier.toString());
            identifier = null;
        }
        if (sequence != null) {
            entryBuilder.setSequence(encodeSequence(sequence));
            sequence = null;
        }
        entryBuilder.setReadLength(previousReadLength);
        if (qualityScores != null) {
            entryBuilder.setQualityScores(ByteString.copyFrom(qualityScores));
            qualityScores = null;
        }

        collectionBuilder.addReads(entryBuilder.build());
        messageChunkWriter.writeAsNeeded(collectionBuilder);
    }

    private synchronized ByteString encodeSequence(final CharSequence sequence) {
        final int length = sequence.length();
        if (length > byteBuffer.length) {
            byteBuffer = new byte[length];
        }
        final byte[] bytes = byteBuffer;
        for (int i = 0; i < length; i++) {
            bytes[i] = (byte) sequence.charAt(i);
            ++sequenceBasesWritten;
        }
        previousReadLength = length;
        return ByteString.copyFrom(bytes, 0, length);
    }

    public void setNumEntriesPerChunk(final int numEntriesPerChunk) {
        messageChunkWriter.setNumEntriesPerChunk(numEntriesPerChunk);
    }

    public synchronized void setIdentifier(final CharSequence identifier) {
        this.identifier = identifier;
    }

    public synchronized long getSequenceBasesWritten() {
        return sequenceBasesWritten;
    }

    public synchronized void printStats(final PrintStream out) {
        messageChunkWriter.printStats(out);
        out.println("Number of bits/base " +
                (messageChunkWriter.getTotalBytesWritten() * 8.0f) / (float) sequenceBasesWritten);
    }
}
