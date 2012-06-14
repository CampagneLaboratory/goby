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

import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
import org.apache.commons.io.FileUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.File;
import java.io.IOException;

/**
 * Constructs a global query index permutation for the concatenation of source alignments.  Assume that some alignments
 * to be concatenated have query index permutations. It is likely that the small indices in each input alignment will
 * conflict with each other. When concatenating alignments, we need to resolve small indices to the original query
 * index of each alignment, then build a new global permutation from the original query index of each entry to
 * a new set of small indices (arranged in sort order for the output alignment). This class takes care of permuting
 * to a new set of indices when it is necessary.
 *
 * @author Fabien Campagne
 *         Date: 6/13/12
 *         Time: 5:48 PM
 */
public class ConcatenatePermutations extends QueryIndexPermutation {
    /**
     * Used to log informational and debug messages.
     */
    private static final Log LOG = LogFactory.getLog(ConcatenatePermutations.class);

    private PermutationReaderInterface[] permReaders;
    private boolean outputNeedsPermutation;
    private String outputFilename;


    public void concatenate(String destinationBasename) {
        if (!outputNeedsPermutation) {
            new File(outputFilename).delete();
        } else {
            close();
            final String destinationFilename = destinationBasename + ".perm";
            try {

                FileUtils.moveFile(new File(outputFilename), new File(destinationFilename));
            } catch (IOException e) {
                LOG.error(String.format("Unable to move temporary permutation file %s to destination %s ", outputFilename, destinationFilename));
            }
        }

    }

    public ConcatenatePermutations(String[] basenames) throws IOException {
        this(basenames, File.createTempFile("permutation", ".perm").getAbsolutePath());
    }

    public ConcatenatePermutations(String[] basenames, String outputFilename) throws IOException {
        super(outputFilename);
        this.outputFilename = outputFilename;
        int i = 0;
        alignmentHasPermutation = new boolean[basenames.length];
        permReaders = new PermutationReaderInterface[basenames.length];
        for (final String basename : basenames) {
            final AlignmentReader reader = new AlignmentReaderImpl(basename);

            try {
                reader.readHeader();
                alignmentHasPermutation[i] = reader.getQueryIndicesWerePermuted();
            } finally {
                reader.close();
            }

            permReaders[i] = alignmentHasPermutation[i] ? new PermutationReader(basenames[i]) : new NoOpPermutationReader();
            outputNeedsPermutation |= alignmentHasPermutation[i];
            i++;
        }


    }

    /**
     * Boolean at index is true when the reader at the given index has query index permutation data.
     */
    private final boolean[] alignmentHasPermutation;

    /**
     * Get the query index corresponding to smallIndex in a source reader, and generate a new global permutation.
     *
     * @param readerIndex Index of the source reader the small index was retried from.
     * @param smallIndex  small index value (or query index value when the source reader has no permutation)
     * @return The value of the new small index
     * @throws IOException
     */
    public final int combine(final int readerIndex, final int smallIndex) throws IOException {

        // Get the original query index from the source reader:
        final int queryIndex = alignmentHasPermutation[readerIndex] ? permReaders[readerIndex].getQueryIndex(smallIndex) : smallIndex;
        // permutate in the global order of the output alignment:
        return permutate(queryIndex);
    }

    /**
     * @return True when at least one input alignment had a query index permutation file (.perm extension).
     */
    public boolean needsPermutation() {
        return outputNeedsPermutation;
    }


}
