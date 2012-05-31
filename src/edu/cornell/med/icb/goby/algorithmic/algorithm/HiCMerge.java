/*
 * Copyright (C) 2009-2012 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
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

package edu.cornell.med.icb.goby.algorithmic.algorithm;

import edu.cornell.med.icb.goby.alignments.*;
import edu.cornell.med.icb.goby.stats.EmpiricalPValueEstimator;
import edu.cornell.med.icb.goby.util.HeaderUtil;
import edu.cornell.med.icb.goby.util.WarningCounter;
import it.unimi.dsi.fastutil.ints.*;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Merge two alignments that are part of a HiC experiment. The query index of the reads must match between
 * the two alignments, which are typically created
 *
 * @author Fabien Campagne
 *         Date: 5/8/12
 *         Time: 2:49 PM
 */
public class HiCMerge {
    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(HiCMerge.class);

    public void merge(List<File> inputFiles, String outputFile) throws IOException {
        if (inputFiles.size() != 2) {
            throw new UnsupportedOperationException("exactly two input files were expected.");
        }
        File fileA = inputFiles.get(0);
        File fileB = inputFiles.get(1);
        AlignmentReader readerA = null;
        AlignmentReader readerB = null;
        AlignmentTooManyHitsReader tmhA = null;
        AlignmentTooManyHitsReader tmhB = null;
        final String basenameA = fileA.getAbsolutePath();
        final String basenameB = fileB.getAbsolutePath();

        try {

            readerA = new AlignmentReaderImpl(basenameA);
            readerB = new AlignmentReaderImpl(basenameB);
        } catch (IOException e) {
            LOG.error("Unable to access an input alignment file.", e);
            throw e;
        }
        AlignmentWriter writer = new AlignmentWriterImpl(outputFile);
        merge(readerA, readerB, writer);
        transferHeader(basenameA, basenameB, writer);
        writer.close();
    }

    protected void transferHeader(String basenameA, String basenameB, AlignmentWriter writer) throws IOException {
        final AlignmentReaderFactory alignmentReaderFactory = new DefaultAlignmentReaderFactory();


        final ConcatAlignmentReader alignmentReader =
                new ConcatAlignmentReader(alignmentReaderFactory, false, new String[]{basenameA, basenameB});
        HeaderUtil.copyHeader(alignmentReader, writer);


    }

    WarningCounter atMost10 = new WarningCounter();

    protected void merge(AlignmentReader readerA, AlignmentReader readerB, AlignmentWriter writer) throws IOException {


        boolean done = false;
        Int2ObjectMap<Alignments.AlignmentEntry> mapA = new Int2ObjectAVLTreeMap<Alignments.AlignmentEntry>();
        Int2ObjectMap<Alignments.AlignmentEntry> mapB = new Int2ObjectAVLTreeMap<Alignments.AlignmentEntry>();
        int N = 1000;
        int k = 0;
        while (!done) {

            done |= readN(readerA, mapA, N);
            done |= readN(readerB, mapB, N);
            writeFlush(mapA, mapB, writer);
            k++;
            if (k % 100000 == 0) {
                System.out.printf("mapA.size=%d mapB.size=%d %n", mapA.size(), mapB.size());
            }
        }
        // write the last ones:
        writeFlush(mapA, mapB, writer);

    }

    /**
     * Find entries that exist in both maps A and B, and writes them to the output.
     *
     * @param mapA   map of query index to entry, for the primary read.
     * @param mapB   map of query index to entry, for the mate.
     * @param writer where output is written.
     * @throws IOException If an error occurs reading the input or writing the output.
     */
    private void writeFlush(final Int2ObjectMap<Alignments.AlignmentEntry> mapA,
                            final Int2ObjectMap<Alignments.AlignmentEntry> mapB,
                            final AlignmentWriter writer) throws IOException {
        final IntSortedSet commonQueryIndices = new IntAVLTreeSet();
        commonQueryIndices.addAll(mapA.keySet());
        commonQueryIndices.retainAll(mapB.keySet());
        for (final int queryIndex : commonQueryIndices) {
            final Alignments.AlignmentEntry.Builder entryA = Alignments.AlignmentEntry.newBuilder(mapA.get(queryIndex));
            final Alignments.AlignmentEntry.Builder entryB = Alignments.AlignmentEntry.newBuilder(mapB.get(queryIndex));
            if (entryA.getFragmentIndex() != entryB.getFragmentIndex()) {
                // OK
            } else {
                atMost10.warn(LOG, "fragment indices must differ for alignment entries in the input files. Forcing different indices");
                entryA.setFragmentIndex(0);
                entryB.setFragmentIndex(1);
            }
            final Alignments.RelatedAlignmentEntry.Builder linkA2B = Alignments.RelatedAlignmentEntry.newBuilder();
            final Alignments.RelatedAlignmentEntry.Builder linkB2A = Alignments.RelatedAlignmentEntry.newBuilder();
            // link A to B:
            linkA2B.setTargetIndex(entryB.getTargetIndex());
            linkA2B.setPosition(entryB.getPosition());
            linkA2B.setFragmentIndex(entryB.getFragmentIndex());
            entryA.setPairAlignmentLink(linkA2B);
            entryA.setPairFlags(EntryFlagHelper.firstInPair() | EntryFlagHelper.paired() | (entryA.getMatchingReverseStrand() ? EntryFlagHelper.readReverseStrand() : 0));
            // link B to A:
            linkB2A.setTargetIndex(entryA.getTargetIndex());
            linkB2A.setPosition(entryA.getPosition());
            linkB2A.setFragmentIndex(entryA.getFragmentIndex());
            entryA.setPairFlags(EntryFlagHelper.secondInPair() | EntryFlagHelper.paired() | (entryB.getMatchingReverseStrand() ? EntryFlagHelper.mateReverseStrand() : 0));
            entryB.setPairAlignmentLink(linkB2A);

            writer.appendEntry(entryA.build());
            writer.appendEntry(entryB.build());
            mapA.remove(queryIndex);
            mapB.remove(queryIndex);
        }

    }

    private boolean readN(final AlignmentReader readerA, final Int2ObjectMap<Alignments.AlignmentEntry> mapA, int N) {
        int count = 0;

        while (count < N && readerA.hasNext()) {

            final Alignments.AlignmentEntry entry = readerA.next();

            if (entry.hasAmbiguity() && entry.getAmbiguity() > 1) {
                // ignore ambiguous entries.
            } else {
                mapA.put(entry.getQueryIndex(), entry);
                count++;
            }

        }
        if (count == N) {
            // not done
            return false;
        } else {
            return true;
        }
    }
}
