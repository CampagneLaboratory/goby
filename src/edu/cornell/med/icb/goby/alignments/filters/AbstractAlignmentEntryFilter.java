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

package edu.cornell.med.icb.goby.alignments.filters;

import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.identifier.IndexedIdentifier;

/**
 * @author Fabien Campagne
 */


/**
 * Abstract class for AlignmentEntryFilter. These assist with merging alignment files.
 *
 * @author Kevin Dorff
 */
public abstract class AbstractAlignmentEntryFilter {
    /**
     * Give the filter access to targets of the merged alignment.
     *
     * @param targets
     */
    public abstract void setHeader(final IndexedIdentifier targets);

    /**
     * Called during first pass processing of every entry.
     *
     * @param entry the entry to inspect for candidacy in the merge
     */
    public abstract void inspectEntry(final Alignments.AlignmentEntry entry);

    /**
     * Called after first pass (inspecEntry) processing is complete.
     */
    public abstract void postProcessing();

    /**
     * Returns true if the entry should be retained.
     *
     * @param entry the entry to be inspected.
     */
    public abstract boolean shouldRetainEntry(final Alignments.AlignmentEntry entry);



}
