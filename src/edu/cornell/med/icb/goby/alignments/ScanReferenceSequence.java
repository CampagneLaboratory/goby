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

package edu.cornell.med.icb.goby.alignments;

import edu.cornell.med.icb.goby.reads.ReadSet;
import edu.cornell.med.icb.goby.reads.ReadsLoader;
import edu.cornell.med.icb.goby.reads.SequenceDigests;
import edu.cornell.med.icb.goby.reads.SequenceEncoder;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.logging.ProgressLogger;

import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: Jun 20, 2009
 *         Time: 12:11:09 PM
 */
public class ScanReferenceSequence {
    private byte[] sequence;
    private ProgressLogger progress;
    private int readLength;
    private long trueMatches;
    private byte[] byteBuffer;
    private SequenceDigests[] readDigests;
    private ObjectList<byte[]> compressedReads;
    private long potentialMatches;
    private int[] hitsPerRead;
    private final SequenceEncoder encoder = new SequenceEncoder();
    private int referenceSequenceIndex;
    private AlignmentWriter writer;
    private ReadSet readIndexFilter;

    public void setWriter(final AlignmentWriter writer) {
        this.writer = writer;
    }

    public ScanReferenceSequence() {
        matchingReferenceIndex = new Int2IntOpenHashMap();
    }

    public Int2IntMap matchingReferenceIndex;

    public void setAlphabet(final String alphabet) {
        this.alphabet = alphabet;
    }

    private String alphabet;

    public void setReadOccurenceThreshold(final int readOccurenceThreshold) {
        this.readOccurenceThreshold = readOccurenceThreshold;
    }

    private int readOccurenceThreshold;

    public void set(final int referenceSequenceIndex,
                    final byte[] sequence,
                    final ReadsLoader loader,
                    final int[] hitsPerRead) {

        this.sequence = sequence;
        this.byteBuffer = loader.getByteBuffer();
        this.readDigests = loader.getDigests();
        this.compressedReads = loader.getCompressedReads();
        this.hitsPerRead = hitsPerRead;
        this.referenceSequenceIndex = referenceSequenceIndex;
    }

    public long getTrueMatches() {
        return trueMatches;
    }

    public long getPotentialMatches() {
        return potentialMatches;
    }

    public ScanReferenceSequence scan() throws IOException {

        progress.itemsName = "x 1000 positions";

        final int length = sequence.length;
        final int end = length - readLength;
        for (int referencePosition = 0; referencePosition < end; ++referencePosition) {

            for (final SequenceDigests sd : readDigests) {
                // the reads have been reverse-complemented so we don't need to reverse complement the reference.
                final int digest = sd.digestDirectStrandOnly(sequence, referencePosition, readLength);
                final int[] readIndices = sd.getReadIndices(digest);
                for (final int readIndex : readIndices) {
                    if (readIndex != -1) {
                        if (hitsPerRead[readIndex] <= readOccurenceThreshold) {

                            final byte[] readCompressedSeq = compressedReads.get(readIndex);
                            if (sd.confirmMatch(readCompressedSeq, sequence, referencePosition, readLength)) {

                                ++trueMatches;
                                hitsPerRead[readIndex] += 1;

                                if (hitsPerRead[readIndex] > readOccurenceThreshold) {
                                    // too many hits already for this read, remove from consideration for further potential hits:
                                    // the digest may also be produced by a sequence that does not have too many hits. We cannot
                                    // remove without checking that the digest is not produced by another read.
                                    sd.remove(digest, readIndex);
                                } else {
                                    final Alignments.AlignmentEntry.Builder builder =
                                            Alignments.AlignmentEntry.newBuilder();
                                    builder.setQueryIndex(readIndex);
                                    builder.setTargetIndex(referenceSequenceIndex);
                                    builder.setPosition(referencePosition);
                                    builder.setMatchingReverseStrand(!sd.isMatchingPositiveStrand());

                                    builder.setScore(readLength);
                                    builder.setNumberOfIndels(0);

                                    builder.setQueryAlignedLength(readLength);
                                    int multiplicity = 1;
                                    // multiplicity of a read is the number of times the sequence
                                    // of the read is identically repeated across a sample file.
                                    // When the sequence is exactly repeated, the alignment would
                                    // yield exactly the same result. In such cases, we do not do
                                    // the alignment, but just keep repeating the alignment
                                    // multiplicity times.

                                    if (readIndexFilter != null) {
                                        // we have a multiplicity filter. Use it to determine multiplicity.
                                        multiplicity = readIndexFilter.getMultiplicity(readIndex);
                                    }
                                    builder.setMultiplicity(multiplicity);
                                    final Alignments.AlignmentEntry entry = builder.build();
                                    writer.appendEntry(entry);
                                }
                            }
                        }
                        ++potentialMatches;
                    }
                }
                //     System.out.printf("A read may perfectly match reference %d at position %d %n", referenceEntry.getReadIndex(), referencePosition);
            }

            if ((referencePosition % 1000) == 1) {
                progress.lightUpdate();
            }
            //progress.lightUpdate();

        }
        return this;
    }


    public void setProgress(final ProgressLogger progress) {
        this.progress = progress;
    }

    public void setReadLength(final int readLength) {
        this.readLength = readLength;
    }

    public void setReadIndexFilter(final ReadSet readIndexFilter) {
        this.readIndexFilter = readIndexFilter;
    }
}
