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

package edu.cornell.med.icb.alignments.filters;

import edu.cornell.med.icb.alignments.Alignments;
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
     * @param k                The maximum number of alignment at best score.
     * @param maxNumberOfReads the maximum number of POSSIBLE reads we could encounter
     */
    public BestScoreOnlyAlignmentFilter(final int k, final int maxNumberOfReads) {
        indexToBestScore = new float[maxNumberOfReads];
        Arrays.fill(indexToBestScore, Float.MIN_VALUE);

    }

    @Override
    public void setHeader(final IndexedIdentifier targets) {
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
