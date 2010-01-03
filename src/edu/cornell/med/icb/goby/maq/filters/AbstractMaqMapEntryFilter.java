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

package edu.cornell.med.icb.goby.maq.filters;

import edu.cornell.med.icb.goby.maq.MaqMapEntry;
import edu.cornell.med.icb.goby.maq.MaqMapHeader;

/**
 * Abstract class for MaqMapentryFilter. These assist with merging Maq Map files.
 * @author Kevin Dorff
 */
public abstract class AbstractMaqMapEntryFilter {
    /**
     * Called when there is a new or updated MaqMapHeader, ie, a new MapMapReader was created.
     * @param header the new header or updated header
     */
    public abstract void setHeader(final MaqMapHeader header);
    /**
     * Called during first pass processing of every entry.
     * @param entry the entry to inspect for candidacy in the merge
     */
    public abstract void inspectEntry(final MaqMapEntry entry);
    /**
     * Called after first pass processing in is complete.
     */
    public abstract void postProcessing();
    /**
     * Called during second pass processing of every entry to find reads to keep / omit.
     * @param entry the entry to check if it should be retained in a merge
     */
    public abstract boolean shouldRetainEntry(final MaqMapEntry entry);
    /**
     * A sanity check, how many records should be written. Generally calculated in
     * postProcessing().
     * @return the number of records that should be written.
     */
    public abstract int getShouldWrite();

    /** One megabyte constant. */
    private static final long ONE_MEGABYTE = 1024 * 1024;

    /**
     * Return the heap size in in use / total as a string.
     * @return heap size.
     */
    public final String getHeapSize() {
        final long totalMemory = Runtime.getRuntime().totalMemory();
        final long freeMemory = Runtime.getRuntime().freeMemory();
        final Long memoryUseMB = (totalMemory - freeMemory) / ONE_MEGABYTE;
        final Long memoryAvailMB = totalMemory / ONE_MEGABYTE;

        return String.format("%,d / %,d MB", memoryUseMB, memoryAvailMB);
    }
}
