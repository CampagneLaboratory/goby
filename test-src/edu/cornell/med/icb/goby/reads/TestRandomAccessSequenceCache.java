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

import it.unimi.dsi.lang.MutableString;
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
                "aCtGNNNACTGMARARRAQA\n";
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

    /**
     * Test creation of a sequence
     * @throws IOException
     */
    @Test
    public void testMakeGenome() throws IOException {
        final int firstPos = 2;
        final String genomeStr = "TGAATGAGACCTA";
        final RandomAccessSequenceCache genome = makeGenome(firstPos, genomeStr);

        assertEquals(firstPos + genomeStr.length(), genome.getLength(0));
        assertEquals('N', genome.get(0, 0));
        assertEquals('N', genome.get(0, 1));
        assertEquals('T', genome.get(0, 2));
        assertEquals('G', genome.get(0, 3));
        assertEquals('A', genome.get(0, 4));
    }

    /**
     * The goal of this is to create a reference where sequence
     * first appears at 0-based position firstPos. Positions prior to
     * this will be padded with N's.
     * This seems to work but the length of targetIndex 0 seems to be
     * off by 1. ??
     * @param firstPos the 0-based position where sequence should start
     * @param sequence the reference sequence to create
     * @return the sequence cache
     * @throws IOException io exception thrown
     */
    private RandomAccessSequenceCache makeGenome(final int firstPos, final CharSequence sequence) throws IOException {
        final MutableString genomeSeq = new MutableString();
        genomeSeq.append(">1\n");
        for (int i = 0; i < firstPos; i++) {
            genomeSeq.append('N');
        }
        genomeSeq.append(sequence).append('\n');
        System.out.println(genomeSeq);
        final RandomAccessSequenceCache cache = new RandomAccessSequenceCache();
        cache.loadFasta(new StringReader(genomeSeq.toString()));
        return cache;
    }

}
