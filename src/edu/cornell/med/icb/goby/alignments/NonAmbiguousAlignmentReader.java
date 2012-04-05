/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.alignments;

import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import java.io.IOException;
import java.util.BitSet;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Properties;

/**
 * An alignment reader that returns only non-ambiguous entries. This class reads the TMH information associated with the
 * alignment basename, generates a bitset to keep track of ambigious query indices, then returns only alignments for
 * entries that the aligner did not consider ambiguous.
 *
 * @author Fabien Campagne
 *         Date: Apr 9, 2011
 *         Time: 4:09:55 PM
 */
public class NonAmbiguousAlignmentReader implements AlignmentReader {
    private static final Logger LOG = Logger.getLogger(NonAmbiguousAlignmentReader.class);

    private Alignments.AlignmentEntry possibleEntry;
    private boolean allQueriesHaveQueryIndexOccurences;
    private int maxLocations = 1;

    public void setMaxLocations(int maxLocations) {
        this.maxLocations = maxLocations;
    }

    public boolean isSorted() {
        return delegate.isSorted();
    }

    public boolean isIndexed() {
        return delegate.isIndexed();
    }

    public String basename() {
        return delegate.basename();
    }

    private final AlignmentReader delegate;
    int numQueries;
    private BitSet ambiguousQueryIndices;

    public NonAmbiguousAlignmentReader(final String basename,
                                       final int startReferenceIndex, final int startPosition,
                                       final int endReferenceIndex, final int endPosition) throws IOException {
        delegate = new AlignmentReaderImpl(basename, startReferenceIndex, startPosition, endReferenceIndex, endPosition);
        readTmh(basename);
    }

    private void readTmh(String basename) throws IOException {
        LOG.debug("start reading TMH for " + basename);
        final ProgressLogger pg = new ProgressLogger(LOG);

        {
            final AlignmentTooManyHitsReader tmh = new AlignmentTooManyHitsReader(basename);
            try {
                LOG.debug("finished loading TMH for " + basename);
                numQueries = tmh.getQueryIndices().size();
                ambiguousQueryIndices = new BitSet(numQueries);
                pg.expectedUpdates = tmh.getQueryIndices().size();
                pg.priority = Level.DEBUG;
                pg.info = "processed %d items";
                pg.start("start reading TMH for " + basename);


                for (int i = 0; i < numQueries; i++) {
                    if (tmh.isQueryAmbiguous(i)) {
                        ambiguousQueryIndices.set(i);
                    }
                    pg.lightUpdate();
                }
                pg.stop("done reading TMH for " + basename);
            } finally {
                tmh.close();
            }
        }
    }

    public NonAmbiguousAlignmentReader(final long startOffset,
                                       final long endOffset,
                                       final String basename) throws IOException {
        delegate = new AlignmentReaderImpl(startOffset, endOffset, basename);
        delegate.readHeader();
        if (!delegate.getHasQueryIndexOccurrences()) {
            LOG.debug("Alignment does not have query index occurrences, using TMH file.");
            readTmh(basename);
        } else {
            LOG.debug("Alignment entries have query index occurrences.");
            allQueriesHaveQueryIndexOccurences = true;
        }
    }

    public NonAmbiguousAlignmentReader(final String basename) throws IOException {
        delegate = new AlignmentReaderImpl(basename);
        delegate.readHeader();
        if (!delegate.getHasQueryIndexOccurrences()) {
            // We have to rely on the TMH (pre-goby 2.0):
            LOG.debug("Alignment does not have query index occurrences, using TMH file.");
            readTmh(basename);
            allQueriesHaveQueryIndexOccurences = false;
        } else {
            LOG.debug("Alignment entries have query index occurrences.");
            allQueriesHaveQueryIndexOccurences = true;
        }
    }

    /**
     * Determine if the reader has another entry that is not ambiguous.
     *
     * @return True if such an entry exists, false otherwise.
     */
    @Override
    public final boolean hasNext() {
        if (possibleEntry != null) {
            return true;
        }

        while (delegate.hasNext()) {
            possibleEntry = delegate.next();
            if (possibleEntry.hasQueryIndexOccurrences()) {
                if (possibleEntry.getQueryIndexOccurrences() <= maxLocations) {
                    return true;
                }

            } else if (!ambiguousQueryIndices.get(possibleEntry.getQueryIndex())) {
                // the query index is not ambiguous.
                return true;
            }
        }
        return false;

    }

    /**
     * Returns the next non-ambiguous entry.
     *
     * @return the next non-ambiguous entry.
     */
    public final Alignments.AlignmentEntry next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }

        final Alignments.AlignmentEntry returnValue = possibleEntry;
        possibleEntry = null;
        return returnValue;

    }

    public Alignments.AlignmentEntry skipTo(int targetIndex, int position) throws IOException {
        if (possibleEntry != null) {
            if (possibleEntry.getTargetIndex() >= targetIndex && possibleEntry.getPosition() >= position) {
                return possibleEntry;
            } else {
                possibleEntry = null;
            }
        }

        while (true) {
            possibleEntry = delegate.skipTo(targetIndex, position);
            if (possibleEntry == null) {
                return null;
            }
            if (allQueriesHaveQueryIndexOccurences) {
                if (possibleEntry.getQueryIndexOccurrences() <= maxLocations) {
                    final Alignments.AlignmentEntry tmp = possibleEntry;
                    possibleEntry = null;
                    // the query index is not ambiguous.
                    return tmp;
                } else {
                    // the last entry was ambiguous, look for a new one that is not.
                    possibleEntry = null;
                }
            } else {
                if (!ambiguousQueryIndices.get(possibleEntry.getQueryIndex())) {
                    final Alignments.AlignmentEntry tmp = possibleEntry;
                    possibleEntry = null;
                    // the query index is not ambiguous.
                    return tmp;
                } else {
                    // the last entry was ambiguous, look for a new one that is not.
                    possibleEntry = null;
                }
            }
        }
    }

    public void reposition(int targetIndex, int position) throws IOException {
        delegate.reposition(targetIndex, position);
    }

    public void remove() {
        delegate.remove();
    }

    public void readHeader() throws IOException {
        delegate.readHeader();
    }

    public void readIndex() throws IOException {
        delegate.readIndex();
    }

    public void close() {
        delegate.close();
        if (!allQueriesHaveQueryIndexOccurences) {
            this.ambiguousQueryIndices.clear();
        }
    }

    public Iterator<Alignments.AlignmentEntry> iterator() {
        return delegate.iterator();
    }

    public Properties getStatistics() {
        return delegate.getStatistics();
    }

    public int getNumberOfAlignedReads() {
        return delegate.getNumberOfAlignedReads();
    }


    public ObjectList<ReferenceLocation> getLocations(int modulo) throws IOException {
        return delegate.getLocations(modulo);
    }

    public boolean isQueryLengthStoredInEntries() {
        return delegate.isQueryLengthStoredInEntries();
    }

    public String getAlignerName() {
        return delegate.getAlignerName();
    }

    public String getAlignerVersion() {
        return delegate.getAlignerVersion();
    }

    public int getSmallestSplitQueryIndex() {
        return delegate.getSmallestSplitQueryIndex();
    }

    public int getLargestSplitQueryIndex() {
        return delegate.getLargestSplitQueryIndex();
    }

    public IndexedIdentifier getTargetIdentifiers() {
        return delegate.getTargetIdentifiers();
    }

    public int[] getTargetLength() {
        return delegate.getTargetLength();
    }

    public int getNumberOfTargets() {
        return delegate.getNumberOfTargets();
    }

    public int getNumberOfQueries() {
        return delegate.getNumberOfQueries();
    }

    public boolean isConstantQueryLengths() {
        return delegate.isConstantQueryLengths();
    }

    public int getConstantQueryLength() {
        return delegate.getConstantQueryLength();
    }

    public IndexedIdentifier getQueryIdentifiers() {
        return delegate.getQueryIdentifiers();
    }

    @Override
    public long getStartByteOffset(int startReferenceIndex, int startPosition) {
        return delegate.getStartByteOffset(startReferenceIndex, startPosition);
    }

    @Override
    public boolean getQueryIndicesWerePermuted() {
        return delegate.getQueryIndicesWerePermuted();
    }


    @Override
    public boolean getHasQueryIndexOccurrences() {
        return delegate.getHasQueryIndexOccurrences();
    }

    @Override
    public long getEndByteOffset(int startReferenceIndex, int startPosition, int endReferenceIndex, int endPosition) {
        return delegate.getEndByteOffset(startReferenceIndex, startPosition, endReferenceIndex, endPosition);
    }

}
