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

package edu.cornell.med.icb.goby.alignments.perms;

import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.io.FastBufferedOutputStream;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.*;
import java.util.Arrays;
import java.util.Collections;

/**
 * Write query index to small index permutations. Goby replaces large query indices by small indices
 * to improve compression of sorted alignment entries files. The permutation writer class takes care
 * of writing the reverse permutation smallIndex->queryIndex, to make it possible to retrieve the read
 * index for these applications that need it (i.e., to retrieve the reads that did not align against a
 * genome).
 * Since many uses of an alignment only need entries, index and header files, this strategy makes it
 * possible to transfer only the relevant data files for the task at hand.
 * @author Fabien Campagne
 *         Date: 3/8/12
 *         Time: 12:41 PM
 */
public class PermutationWriter implements Closeable {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(PermutationWriter.class);

    private final String basename;

    private int lastSmallIndexWritten = -1;
    private final DataOutputStream output;
    private int numWritten;
    private int chunkSize = 1000;

    public PermutationWriter(final String basename) {
        this.basename = AlignmentReaderImpl.getBasename(basename);
        DataOutputStream o = null;
        try {
            o = new DataOutputStream(new FastBufferedOutputStream(new FileOutputStream(basename + ".perm")));

        } catch (FileNotFoundException e) {
            LOG.error(e);
            o = null;
        }
        output = o;
    }

    public void close() {

        edu.cornell.med.icb.goby.util.IOUtil.closeQuietly(output);

    }

    /**
     * Append a piece of the permutation.
     *
     * @param permPiece map from query indices to small indices.
     * @throws IOException
     */
    public void append(final Int2IntMap permPiece) throws IOException {

        if (permPiece.isEmpty()) {
            return;
        }
        //final Int2IntSortedMap reverse=new Int2IntAVLTreeMap();

        final IntArrayList smallIndices = new IntArrayList();
        smallIndices.addAll(permPiece.values());
        Collections.sort(smallIndices);
        final int size = smallIndices.size();
        final int minSmallIndexInChunk = smallIndices.get(0);
        final int maxSmallIndexInChunk = smallIndices.get(size - 1);
        IntArrayList largeIndices = new IntArrayList();
        largeIndices.size(permPiece.size());
        for (int queryIndex : permPiece.keySet()) {
            final int smallIndex = permPiece.get(queryIndex);
            final int indexOf = Collections.binarySearch(smallIndices,smallIndex);
            largeIndices.set(indexOf, queryIndex);
        }
        chunk(smallIndices, minSmallIndexInChunk, maxSmallIndexInChunk, largeIndices);
        output.flush();
    }

    private int previousWritten = 0;

    private void chunk(final IntArrayList smallIndices, final int minSmallIndexInChunk,
                       final int maxSmallIndexInChunk, final IntArrayList largeIndices) throws IOException {
        final int max = smallIndices.size();
        int breakpointIndex = -1;
        while (previousWritten < max) {
            breakpointIndex = getBreakPoint(Math.max(0, previousWritten), smallIndices);

            write(previousWritten, smallIndices, breakpointIndex, largeIndices);
            //  lastSmallIndexWritten = smallIndices.get(breakpointIndex - 1);
            previousWritten = breakpointIndex;
        }
    }

    private void write(int previousWritten, IntArrayList smallIndices, int breakpointIndex, IntArrayList largeIndices) throws IOException {
        final int firstSmallIndex = smallIndices.get(previousWritten);
        writeOnePiece(previousWritten, firstSmallIndex, largeIndices, breakpointIndex - previousWritten);

        output.flush();
    }

    public void writeOnePiece(int firstSmallIndex, IntArrayList largeIndices) throws IOException {

        writeOnePiece(0, firstSmallIndex, largeIndices, largeIndices.size());
        output.flush();
    }

    public void writeOnePiece(int previousWritten, int firstSmallIndex, IntArrayList largeIndices, int n) throws IOException {

        output.writeInt(n);
        output.writeInt(firstSmallIndex);

        final int pastLast = n + previousWritten;
        for (int i = previousWritten; i < pastLast; i++) {

            output.writeInt(largeIndices.get(i));
        }
        numWritten += n;

    }

    protected int getBreakPoint(int previousWritten, final IntArrayList smallIndices) throws IOException {
        return getBreakPoint(previousWritten, smallIndices, chunkSize);
    }

    /**
     * Return the index in  smallIndices where there is a gap in sequence, chunkSize,
     * whichever is smaller.
     *
     * @param previousWritten index of the position where to start looking for a breakpoint.
     * @param smallIndices    list of indices.
     * @param chunkSize
     * @return
     * @throws IOException
     */
    protected int getBreakPoint(final int previousWritten, final IntArrayList smallIndices, final int chunkSize) throws IOException {
        if (previousWritten >= smallIndices.size()) {
            return smallIndices.size();
        }
        int index = 0;
        int smallIndex;
        int previousSmallIndex = smallIndices.get(previousWritten);
        int numInChunk = 0;
        int n = Math.min(chunkSize, smallIndices.size() - previousWritten);
        while (numInChunk < n) {
            final int offset = index + previousWritten;
            smallIndex = smallIndices.get(offset);
            if (numInChunk > 0 && smallIndex != previousSmallIndex + 1) {
                // the index of the small index element immediately following the break:
                return offset;
            }

            previousSmallIndex = smallIndex;

            index++;
            numInChunk++;
        }
        return index + previousWritten;
    }
}
