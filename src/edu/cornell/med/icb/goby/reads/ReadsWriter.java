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
    private int barcodeIndex = -1;
    private CharSequence pairSequence;
    private byte[] qualityScoresPair;


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

    public synchronized void setPairSequence(final CharSequence sequence) {
        this.pairSequence = sequence;
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

    /**
     * Append an entry with the next available readindex.
     *
     * @throws IOException If an error occurs while writing the file.
     */
    public synchronized void appendEntry() throws IOException {
        appendEntry(readIndex);
        readIndex++;
    }

    /**
     * Append an entry with a specific read index.
     *
     * @param readIndex Index of the read that will be written
     * @throws IOException If an error occurs while writing the file.
     */
    public synchronized void appendEntry(final int readIndex) throws IOException {

        final Reads.ReadEntry.Builder entryBuilder = Reads.ReadEntry.newBuilder();

        entryBuilder.setReadIndex(readIndex);

        // set current read index to enable interleaving calls to appendEntry(readIndex)/appendEntry().
        this.readIndex = readIndex;
        if (barcodeIndex != -1) {
            entryBuilder.setBarcodeIndex(barcodeIndex);
        }
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
        if (pairSequence != null) {
            entryBuilder.setSequencePair(encodeSequence(pairSequence));
            pairSequence = null;
            entryBuilder.setReadLengthPair(previousReadLength);
        }

        if (qualityScores != null) {
            entryBuilder.setQualityScores(ByteString.copyFrom(qualityScores));
            qualityScores = null;
        }

        if (qualityScoresPair != null) {
            entryBuilder.setQualityScoresPair(ByteString.copyFrom(qualityScoresPair));
            qualityScoresPair = null;
        }

        collectionBuilder.addReads(entryBuilder.build());
        messageChunkWriter.writeAsNeeded(collectionBuilder);
        barcodeIndex = -1;
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

    public void setBarcodeIndex(final int barcodeIndex) {
        this.barcodeIndex = barcodeIndex;
    }

    /**
     * Set quality scores for the second sequence in a pair.
     * @param qualityScores quality Scores in Phred scale.
     */
    public void setQualityScoresPair(final byte[] qualityScores) {
        this.qualityScoresPair=qualityScores;
    }
}
