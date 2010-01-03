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
 * This class assists with filtering MaqMapEntry values when merging non-transcript-based
 * alignments.
 * @author Kevin Dorff
 */
public class FewestMismatchesFilter extends AbstractMaqMapEntryFilter {
    /** The registered indexes - which read-name-index values have been already inspected. */
    private int maxRegisteredIndex;

    /** A array of of read-name-index to the number of fewest mismatches. */
    private final short[] indexToFewestMismatches;

    /** A array of read-name-index to the count of reads with the same fewest mismatches. */
    private final int[] indexToCountAtFewestMismatches;

    /** The number times inspectEntry as called. */
    private int numInspected;

    /** The k value for the filter. */
    private final int k;

    /** The number of entries that should be written. */
    private int shouldWrite;

    /** The start time. */
    private long startTime;

    /**
     * Constructor.
     * @param k the k value
     * @param maxNumberOfReads the maximum number of POSSIBLE reads we could encounter
     */
    public FewestMismatchesFilter(final int k, final int maxNumberOfReads) {
        maxRegisteredIndex = -1;
        indexToFewestMismatches = new short[maxNumberOfReads];
        for (int i = 0; i < maxNumberOfReads; i++) {
            indexToFewestMismatches[i] = -1;
        }
        indexToCountAtFewestMismatches = new int[maxNumberOfReads];
        for (int i = 0; i < maxNumberOfReads; i++) {
            indexToCountAtFewestMismatches[i] = -1;
        }
        numInspected = 0;
        this.k = k;
        this.startTime = System.currentTimeMillis();
    }

    /**
     * Set the new / updated header that is being used when filtering these entries.
     * @param header the header.
     */
    @Override
    public void setHeader(final MaqMapHeader header) {
        // Nothing to do
    }

    /**
     * Inspect an entry (will be called during a first pass of reading the entries).
     * The lowest number of mismatches for a read-name-index is noted and the number of
     * entries at that lowest number of mismatches is noted for that read-name-index.
     * @param entry the entry
     */
    @Override
    public synchronized void inspectEntry(final MaqMapEntry entry) {
        final int index = entry.getReadNameIndex();
        final short numMisMatches = entry.getNumMisMatches();
        numInspected++;
        if (index > maxRegisteredIndex) {
            // new index we haven't observed before
            maxRegisteredIndex = index;
            indexToFewestMismatches[index] = numMisMatches;
            indexToCountAtFewestMismatches[index] = 1;
        } else {
            final short prevQuality = indexToFewestMismatches[index];
            if (prevQuality == numMisMatches) {
                // Increment the count for this quality
                indexToCountAtFewestMismatches[index] += 1;
            } else  if (prevQuality > numMisMatches) {
                // We have a new LOWER quality score, start over with this as the
                // new best quality
                indexToFewestMismatches[index] = numMisMatches;
                indexToCountAtFewestMismatches[index] = 1;
            }
        }
        if ((numInspected % 100000) == 0) {
            // Track how long it took to run the last 100,000 records (# per second)
            final long curTime = System.currentTimeMillis();
            final int msPer100k = (int) (curTime - startTime);
            startTime = curTime;
            System.out.printf(
                    "#Inspected = %,d, #Reg'd = %,d, ms/100k=%,d, memory used / total = %s%n",
                    numInspected, maxRegisteredIndex, msPer100k, getHeapSize());
        }
    }

    /**
     * Determine if this entry should be retained (will be called during a second
     * pass of reading the entries).
     * @param entry the entry
     * @return true if it should be retained
     */
    @Override
    public synchronized boolean shouldRetainEntry(final MaqMapEntry entry) {
        final int index = entry.getReadNameIndex();
        final short numMisMatches = entry.getNumMisMatches();
        final short keepNumMismatches = indexToFewestMismatches[index];
        if (keepNumMismatches == -1) {
            return false;
        }
        if (keepNumMismatches == numMisMatches) {
            int count = indexToCountAtFewestMismatches[index];
            count--;
            if (count == 0) {
                // This is the last entry at this quality
                indexToFewestMismatches[index] = -1;
                indexToCountAtFewestMismatches[index] = -1;
            } else {
                // At least one more entry at this quality
                indexToCountAtFewestMismatches[index] = count;
            }
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
        assert k >= 1 : "k value must be at least 1";
        shouldWrite = 0;
        for (int index = 0; index <= maxRegisteredIndex; index++) {
            final int numDupes = indexToCountAtFewestMismatches[index];
            if (numDupes > k) {
                indexToFewestMismatches[index] = -1;
                indexToCountAtFewestMismatches[index] = -1;
            } else {
                shouldWrite += numDupes;
            }
        }
    }

    /**
     * Get the number of entries that should be written (sanity check).
     * @return the number of entries that should be written
     */
    @Override
    public int getShouldWrite() {
        return shouldWrite;
    }
}
