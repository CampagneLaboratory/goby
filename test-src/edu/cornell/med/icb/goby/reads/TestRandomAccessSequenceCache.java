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

package edu.cornell.med.icb.goby.reads;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;
import org.junit.Test;

import java.io.IOException;
import java.io.StringReader;

/**
 * @author Fabien Campagne
 *         Date: May 19, 2009
 *         Time: 4:01:04 PM
 */
public class TestRandomAccessSequenceCache {
    @Test
    public void testEncodeDecodeOneSequence() throws IOException {

        final String seqs = ">1\n" +
                "ACTGNNNACTG\n";
        final RandomAccessSequenceCache cache = new RandomAccessSequenceCache();
        cache.loadFasta(new StringReader(seqs));
        int i = 0;
        assertEquals('A', cache.get(0, i++));
        assertEquals('C', cache.get(0, i++));
        assertEquals('T', cache.get(0, i++));
        assertEquals('G', cache.get(0, i++));
        assertEquals('N', cache.get(0, i++));
        assertEquals('N', cache.get(0, i++));
        assertEquals('N', cache.get(0, i++));
        assertEquals('A', cache.get(0, i++));
        assertEquals('C', cache.get(0, i++));
        assertEquals('T', cache.get(0, i++));
        assertEquals('G', cache.get(0, i++));
        try {
            // past the end of the sequence.
            cache.get(0, i++);
            fail();
        } catch (AssertionError e) {
            // OK, an assertion was triggered.
        }
    }

    @Test
    public void testEncodeDecodeMoreSequences() throws IOException {

        final String seqs = ">1\n" +
                "ACTGNNNACTG\n" +
                ">2\n" +
                "NNNNNANNNNN";
        final RandomAccessSequenceCache cache = new RandomAccessSequenceCache();
        cache.loadFasta(new StringReader(seqs));
        int i = 0;
        assertEquals('A', cache.get(0, i++));
        assertEquals('C', cache.get(0, i++));
        assertEquals('T', cache.get(0, i++));
        assertEquals('G', cache.get(0, i++));
        assertEquals('N', cache.get(0, i++));
        assertEquals('N', cache.get(0, i++));
        assertEquals('N', cache.get(0, i++));
        assertEquals('A', cache.get(0, i++));
        assertEquals('C', cache.get(0, i++));
        assertEquals('T', cache.get(0, i++));
        assertEquals('G', cache.get(0, i++));


        i = 0;
        assertEquals('N', cache.get(1, i++));
        assertEquals('N', cache.get(1, i++));
        assertEquals('N', cache.get(1, i++));
        assertEquals('N', cache.get(1, i++));
        assertEquals('N', cache.get(1, i++));
        assertEquals('A', cache.get(1, i++));
        assertEquals('N', cache.get(1, i++));
        assertEquals('N', cache.get(1, i++));
        assertEquals('N', cache.get(1, i++));
        assertEquals('N', cache.get(1, i++));
        assertEquals('N', cache.get(1, i++));
    }
}
