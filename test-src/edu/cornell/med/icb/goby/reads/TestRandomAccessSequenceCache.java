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
                "ACTGNNNACTGMARARRAQA\n";
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

        // Non GATC bases become N's
        assertEquals('N', cache.get(0, i++));
        assertEquals('A', cache.get(0, i++));
        assertEquals('N', cache.get(0, i++));
        assertEquals('A', cache.get(0, i++));
        assertEquals('N', cache.get(0, i++));
        assertEquals('N', cache.get(0, i++));
        assertEquals('A', cache.get(0, i++));
        assertEquals('N', cache.get(0, i++));
        assertEquals('A', cache.get(0, i++));

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
