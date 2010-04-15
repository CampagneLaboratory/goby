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

package edu.cornell.med.icb.goby.alignments.filters;

import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.util.Arrays;

/**
 * @author Kevin Dorff
 * @author Fabien Campagne
 */
public class BestScoreAmbiguityAlignmentFilter extends AbstractAlignmentEntryFilter {
    /**
     * Used to log debug and informational messages.
     */
    public static final Log LOG = LogFactory.getLog(BestScoreAmbiguityAlignmentFilter.class);

    /**
     * An array of of read-name-index to the number of fewest mismatches.
     */
    private final float[] indexToBestScore;

    /**
     * A array of read-name-index to the count of reads with the same fewest mismatches.
     */
    private final short[] indexToCountAtBestScore;

    /**
     * The k value for the filter.
     */
    private final int k;

    /**
     * Number of entries that will pass the filter. Populated after postProcessing() has been called.
     *
     * @return The number of times shouldRetainEntry will return true.
     */
    public int getShouldWrite() {
        return shouldWrite;
    }

    /**
     * The number of entries that should be written.
     */
    private int shouldWrite;

    /**
     * The start time.
     */
    private final long startTime;

    private int willSkip;

    /**
     * Constructor.
     *
     * @param k                The maximum number of alignment at best score.
     * @param maxNumberOfReads the maximum number of POSSIBLE reads we could encounter
     */
    public BestScoreAmbiguityAlignmentFilter(final int k, final int maxNumberOfReads) {

        indexToBestScore = new float[maxNumberOfReads];
        Arrays.fill(indexToBestScore, Float.MIN_VALUE);
        indexToCountAtBestScore = new short[maxNumberOfReads];
        Arrays.fill(indexToCountAtBestScore, (short) 0);

        this.k = k;
        this.startTime = System.currentTimeMillis();
    }

    public void setHeader(final AlignmentReader reader) {
        // do nothing
    }

    @Override
    public void setTargetIdentifiers(final IndexedIdentifier targets) {

    }

    /**
     * Inspect an entry (will be called during a first pass of reading the entries).
     * The lowest number of mismatches for a read-name-index is noted and the number of
     * entries at that lowest number of mismatches is noted for that read-name-index.
     *
     * @param entry the entry
     */
    @Override
    public void inspectEntry(final Alignments.AlignmentEntry entry) {
        final int index = entry.getQueryIndex();
        final float score = entry.getScore();


        final float previousScore = indexToBestScore[index];
        if (previousScore == score) {
            // Increment the count for this quality
            indexToCountAtBestScore[index] += 1;

        } else if (previousScore < score) {
            // We have a new best score: start over with this as the new best score
            indexToBestScore[index] = score;
            indexToCountAtBestScore[index] = 1;
        }
    }

    /**
     * Determine if this entry should be retained (will be called during a second
     * pass of reading the entries).
     *
     * @param entry the entry
     * @return true if it should be retained
     */
    @Override
    public boolean shouldRetainEntry(final Alignments.AlignmentEntry entry) {
        final int index = entry.getQueryIndex();
        final float score = entry.getScore();
        final float keepHighestScore = indexToBestScore[index];
        if (keepHighestScore == Float.MIN_VALUE) {
            return false;
        }
        if (keepHighestScore == score) {
            final short count = indexToCountAtBestScore[index];
            if (count > k) {
                // the entry matches too many reference locations. We do not keep it.
                return false;
            } else {
                return true;
            }
        } else {
            // This entry isn't at the appropriate quality.
            return false;
        }
    }

    /**
     * Post processing (after all inspectEntry's have been called).
     * Filter the indexToFewestMismatchesMap / indexToCountAtFewestMismatchesMap maps
     * to help control the number of items that are written.  This will also populate
     * the shouldWrite field.
     */
    @Override
    public synchronized void postProcessing() {
    }

    public int getWillSkip() {
        return willSkip;
    }
}
