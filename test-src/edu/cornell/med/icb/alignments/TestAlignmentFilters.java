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

package edu.cornell.med.icb.alignments;

import edu.cornell.med.icb.alignments.filters.PercentMismatchesQualityFilter;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import org.junit.Test;

/**
 * @author Fabien Campagne
 *         Date: May 6, 2009
 *         Time: 11:44:09 AM
 */
public class TestAlignmentFilters {
    @Test
    public void testFewMismatchesFilter() {

        Alignments.AlignmentEntry.Builder builder;
        Alignments.AlignmentEntry entry;

        builder = buildMinimalEntry();

        builder.setNumberOfIndels(0);
        builder.setNumberOfMismatches(0);
        entry = builder.build();

        PercentMismatchesQualityFilter filter = new PercentMismatchesQualityFilter();
        assertTrue("entry must be kept", filter.keepEntry(100, entry));

        builder = buildMinimalEntry();
        builder.setNumberOfIndels(2);
        builder.setNumberOfMismatches(3);
        builder.setQueryAlignedLength(1);
        entry = builder.build();

        filter = new PercentMismatchesQualityFilter();
        assertFalse("entry must be removed", filter.keepEntry(1, entry));
        assertFalse("entry must be removed", filter.keepEntry(10, entry));
        assertTrue("entry must be kept", filter.keepEntry(100, entry));
        assertTrue("entry must be kept", filter.keepEntry(101, entry));

        builder = buildMinimalEntry();
        builder.setNumberOfIndels(2);
        builder.setNumberOfMismatches(3);

        entry = builder.build();

        filter = new PercentMismatchesQualityFilter();
        filter.setParameters("threshold=0.02");         // 2% differences allowed.
        assertFalse("entry must be removed", filter.keepEntry(1, entry));
        assertFalse("entry must be removed", filter.keepEntry(10, entry));
        assertFalse("entry must be removed", filter.keepEntry(100, entry));
        assertFalse("entry must be removed", filter.keepEntry(101, entry));
        assertTrue("entry must be kept", filter.keepEntry(250, entry));
        assertTrue("entry must be kept", filter.keepEntry(251, entry));
    }
    
    private Alignments.AlignmentEntry.Builder buildMinimalEntry() {
        final Alignments.AlignmentEntry.Builder builder;
        builder = Alignments.AlignmentEntry.newBuilder();

        builder.setQueryIndex(0);
        builder.setTargetIndex(0);
        builder.setPosition(10);
        builder.setQueryAlignedLength(100);
        builder.setMatchingReverseStrand(false);

        return builder;
    }
}
