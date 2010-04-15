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

import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.identifier.IndexedIdentifier;

import java.util.Arrays;

/**
 * Filter out reference matches if that do not have the best score of each alignment with the
 * same query.
 * @author Fabien Campagne
 */
public class BestScoreOnlyAlignmentFilter extends AbstractAlignmentEntryFilter {
    /**
     * A array of of read-name-index to the number of fewest mismatches.
     */
    private final float[] indexToBestScore;

    /**
     * Constructor.
     *
     * @param maxNumberOfReads the maximum number of POSSIBLE reads we could encounter
     */
    public BestScoreOnlyAlignmentFilter(final int maxNumberOfReads) {
        super();
        indexToBestScore = new float[maxNumberOfReads];
        Arrays.fill(indexToBestScore, Float.MIN_VALUE);
    }

    @Override
    public void setTargetIdentifiers(final IndexedIdentifier targets) {
        // do nothing
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
        if (previousScore < score) {
            // update the best score for this query
            indexToBestScore[index] = score;
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
            return true;
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
}
