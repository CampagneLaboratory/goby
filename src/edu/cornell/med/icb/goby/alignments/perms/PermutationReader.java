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
import it.unimi.dsi.fastutil.io.FastBufferedInputStream;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.*;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;

/**
 * Reads permutations associated with Goby alignments. Permutations store the mapping between small query indices
 * and large query indices (rank of the read in the original read file).
 *
 * @author Fabien Campagne
 *         Date: 3/9/12
 *         Time: 2:52 PM
 * @see PermutationWriter
 */
public class PermutationReader implements PermutationReaderInterface {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(PermutationWriter.class);

    private FastBufferedInputStream input;
    private String basename;
    private static final Comparator<? super Block> SMALL_INDEX_COMPARATOR = new CompareSmallIndex();
    private DataInputStream dataInput;

    public PermutationReader(String basename) throws IOException {
        this.basename = AlignmentReaderImpl.getBasename(basename);
        FastBufferedInputStream inputStream = null;
        final String filename = basename + ".perm";

        try {

            inputStream = new FastBufferedInputStream(new FileInputStream(filename));
            dataInput = new DataInputStream(inputStream);
            input = inputStream;
            makeIndex(inputStream);
        } catch (FileNotFoundException e) {
            final String error = String.format("A permutation file called %s could not be found, but is required to reconstruct original query indices to complete this task.",
                    filename);
            LOG.error(error, e);
            inputStream = null;
            assert inputStream!=null:error;
        }

    }

    private void makeIndex(FastBufferedInputStream inputStream) throws IOException {
        input.position(0);
        final ObjectArrayList<Block> blocks = new ObjectArrayList<Block>();

        final DataInputStream dataInput = new DataInputStream(new FastBufferedInputStream(new FileInputStream(basename + ".perm")));
        try {
            long offset = 0;

            while (dataInput.available() > 0) {

                final Block block = new Block();
                block.offset = offset;
                block.n = dataInput.readInt();
                block.firstSmallIndex = dataInput.readInt();
                dataInput.skip(block.n * 4L);
                blocks.add(block);
                offset += block.n * 4L + 8L;
            }
            Collections.sort(blocks, SMALL_INDEX_COMPARATOR);
            indexBlocks = blocks.toArray(new Block[blocks.size()]);
        } finally {
            dataInput.close();
        }
    }

    // these fields implement an index into the disk permutation data structure

    private Block[] indexBlocks;
    // used to search indexedBlocks for a block that contains a small index:
    private Block singletonQuery = new Block();

    static class Block {
        // the first small index in the block
        int firstSmallIndex;
        int n;
        // offset into the perm file.
        long offset;
        // available on disk at offset:
        // queryIndices for the n small indices    [firstSmallIndex firstSmallIndex+1 ... firstSmallIndex+n]
        // int []queryIndices;
    }

    /**
     * Return the query index associated with a small index, or -1 if the association was not defined.
     *
     * @param smallIndex
     * @return
     * @throws IOException
     */
    @Override
    public int getQueryIndex(int smallIndex) throws IOException {
        //   Util.invertPermutation()
        singletonQuery.firstSmallIndex = smallIndex;
        final int ip = Arrays.binarySearch(indexBlocks, singletonQuery, SMALL_INDEX_COMPARATOR);
        final Block block;
        if (ip >= 0) {
            block = indexBlocks[ip];
            seek(block.offset + 4);
            final int first = dataInput.readInt();
            assert first == smallIndex : "assertion failed at smallIndex=" + smallIndex;
            final int queryIndex = dataInput.readInt();
            return queryIndex;
        } else {
            final int index = -(ip + 1);
            final int off = index - 1;
            if (off < 0) {
                return -1;
            }
            block = indexBlocks[off];
            final int soffset = smallIndex - block.firstSmallIndex;
            if (soffset < 0 || block.firstSmallIndex + block.n < smallIndex) {
                // not in the block
                return -1;
            }
            seek(block.offset + 8 + soffset * 4L);
            final int queryIndex = dataInput.readInt();
            return queryIndex;
        }

    }

    private void seek(long offset) throws IOException {
        input.position(offset);
    }

    private static class CompareSmallIndex implements Comparator<Block> {
        @Override
        public int compare(Block a, Block b) {
            return a.firstSmallIndex - b.firstSmallIndex;
        }
    }
}
