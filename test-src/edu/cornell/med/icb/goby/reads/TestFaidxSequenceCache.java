/*
 * Copyright (C) 2009-2012 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.reads;

import org.junit.Test;

import java.io.IOException;

import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

/**
 * @author Fabien Campagne
 *         Date: May 19, 2009
 *         Time: 4:01:04 PM
 */
public class TestFaidxSequenceCache {
    @Test
    public void testLoadDual() {
        try {
            final DualRandomAccessSequenceCache cache = new DualRandomAccessSequenceCache();
            cache.load("test-data/faidx/file1.fasta");
            assertEquals(2, cache.size());

            final DualRandomAccessSequenceCache cache2 = new DualRandomAccessSequenceCache();
            cache2.load("test-data/faidx/file1");
            assertEquals(2, cache2.size());

        } catch (ClassNotFoundException e) {
            assertTrue(false);
        } catch (IOException e) {
            assertTrue(false);
        }
    }

    @Test
    public void testEncodeDecodeOneSequence() throws IOException {


        final RandomAccessSequenceInterface cache = new PicardFastaIndexedSequence("test-data/faidx/file1.fasta");

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


        final RandomAccessSequenceInterface cache = new PicardFastaIndexedSequence("test-data/faidx/file1.fasta");
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

    @Test
    public void testRemoveNewLines() throws IOException {


        final PicardFastaIndexedSequence cache = new PicardFastaIndexedSequence("test-data/faidx/file2.fasta");
        cache.print(0);
        int i = 0;
        do {
            i = testPattern(cache, i);
        } while (i < cache.getLength(0));
    }

    private int j = 0;

    private int testPattern(PicardFastaIndexedSequence cache, int i) {
        char[] bases = {'A', 'C', 'T', 'G'};

        if (j >= bases.length) {
            j = 0;
        }
        final char expected = bases[j++];
        assertEquals('A', cache.get(0, i++));

        final char actual = cache.get(0, i++);
        assertEquals(String.format("middle base (actual=%c) differs from expected (%c) at position %d", actual, expected, i), expected, actual);

        assertEquals('C', cache.get(0, i++));
        return i;
    }
}
