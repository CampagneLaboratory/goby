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

import edu.cornell.med.icb.identifier.IndexedIdentifier;
import edu.cornell.med.icb.goby.maq.MaqMapEntry;
import it.unimi.dsi.lang.MutableString;
import static org.junit.Assert.assertEquals;
import org.junit.Test;

/**
 * Test the FewestMismatchesFilter class.
 * @author Kevin Dorff
 */
public class TestFewestMismatchesFilter {
    /** Index identifiers for making nameIndex for MaqMapEntry constructor. */
    private IndexedIdentifier indexedIdentifiers = new IndexedIdentifier();

    /**
     * General tests for class.
     */
    @Test
    public void testFewestMismatchesFilter() {
        final FewestMismatchesFilter store = new FewestMismatchesFilter(2, 10);

        // 1 is inspected 3 times
        MaqMapEntry entry = makeEntry("zero", 0);
        store.inspectEntry(entry);
        assertEquals(0, entry.getReadNameIndex());

        // 2 is inspected 2 times
        entry = makeEntry("one", 0);
        store.inspectEntry(entry);
        assertEquals(1, entry.getReadNameIndex());

        entry = makeEntry("zero", 0);
        store.inspectEntry(entry);
        assertEquals(0, entry.getReadNameIndex());

        // 3 is inspected 1 times
        entry = makeEntry("two", 0);
        store.inspectEntry(entry);
        assertEquals(2, entry.getReadNameIndex());

        entry = makeEntry("one", 0);
        store.inspectEntry(entry);
        assertEquals(1, entry.getReadNameIndex());

        entry = makeEntry("zero", 0);
        store.inspectEntry(entry);
        assertEquals(0, entry.getReadNameIndex());

        store.postProcessing();

        assertEquals(false, store.shouldRetainEntry(makeEntry("zero", 0)));
        assertEquals(true, store.shouldRetainEntry(makeEntry("one", 0)));
        assertEquals(true, store.shouldRetainEntry(makeEntry("two", 0)));
        assertEquals(false, store.shouldRetainEntry(makeEntry("not_there", 0)));

        assertEquals(3, store.getShouldWrite());
    }

    /**
     * Make a simple MaqMapEntry.
     * @param name read name
     * @param numMisMatches number of mismatches
     * @return the new entry
     */
    private MaqMapEntry makeEntry(final String name, final int numMisMatches) {
        final MaqMapEntry entry = new MaqMapEntry(128, indexedIdentifiers);
        entry.setInfo1((short) numMisMatches);
        entry.setReadName(new MutableString(name));
        return entry;
    }
}
