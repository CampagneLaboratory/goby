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

package edu.cornell.med.icb.goby.alignments;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import org.junit.Test;

import java.io.IOException;
import java.io.StringReader;
import java.util.NoSuchElementException;

public class TestLastReader {
    private static final double EPSILON = 1.0E-10;

    @Test(expected = NoSuchElementException.class)
    public void testRead() throws IOException {
        final String mafInput = getMafInput();
        final LastParser reader = new LastParser(new StringReader(mafInput));
        assertTrue(reader.hasNext());
        reader.next();
        assertEquals(42.0f, reader.getScore(), EPSILON);
        assertEquals(2, reader.getAlignedSequences().size());
        assertEquals("Y", reader.getAlignedSequences().get(0).sequenceIdentifier.toString());
        assertEquals(2709736, reader.getAlignedSequences().get(0).alignedStart);
        assertEquals(42, reader.getAlignedSequences().get(0).alignedLength);
        assertEquals('+', reader.getAlignedSequences().get(0).strand);
        assertEquals(57772954, reader.getAlignedSequences().get(0).sequenceLength);
        assertEquals("CGCTTGCAGTGAGCTGAGATTGCGCCACTGCACTCCAGCCTG", reader.getAlignedSequences().get(0).alignment.toString());

        assertEquals("1", reader.getAlignedSequences().get(1).sequenceIdentifier.toString());
        assertEquals(0, reader.getAlignedSequences().get(1).alignedStart);
        assertEquals(42, reader.getAlignedSequences().get(1).alignedLength);
        assertEquals('+', reader.getAlignedSequences().get(1).strand);
        assertEquals(42, reader.getAlignedSequences().get(1).sequenceLength);
        assertEquals("CGCTTGCAGTGAGCTGAGATTGCGCCACTGCACTCCAGCCTG", reader.getAlignedSequences().get(1).alignment.toString());


        assertTrue(reader.hasNext());
        reader.next();
        assertEquals(41.0f, reader.getScore(), EPSILON);
        assertEquals("Y", reader.getAlignedSequences().get(0).sequenceIdentifier.toString());
        assertEquals("1", reader.getAlignedSequences().get(1).sequenceIdentifier.toString());
        assertEquals(2, reader.getAlignedSequences().size());

        assertTrue(reader.hasNext());
        reader.next();
        assertEquals(41.0f, reader.getScore(), EPSILON);
        assertEquals(2, reader.getAlignedSequences().size());

        assertTrue(reader.hasNext());
        reader.next();
        assertEquals("Y", reader.getAlignedSequences().get(0).sequenceIdentifier.toString());
        assertEquals("1", reader.getAlignedSequences().get(1).sequenceIdentifier.toString());

        assertEquals(40.0f, reader.getScore(), EPSILON);
        assertEquals(2, reader.getAlignedSequences().size());

        assertFalse(reader.hasNext());

        // next should throw an exception if hasNext is false
        reader.next();
    }


    private String getMafInput() {
        return "# LAST version 50\n" +
                "#\n" +
                "# a=7 b=1 c=100000 e=40 d=24 x=27 y=10\n" +
                "# u=0 s=2 m=10 l=1 k=1 i=134217728 w=1000 t=-1 g=1 j=3\n" +
                "# humandb\n" +
                "#\n" +
                "#    A  C  G  T\n" +
                "# A  1 -1 -1 -1\n" +
                "# C -1  1 -1 -1\n" +
                "# G -1 -1  1 -1\n" +
                "# T -1 -1 -1  1\n" +
                "#\n" +
                "# Coordinates are 0-based.  For - strand matches, coordinates\n" +
                "# in the reverse complement of the 2nd sequence are used.\n" +
                "#\n" +
                "# name start alnSize strand seqSize alignment\n" +
                "#\n" +
                "a score=42\n" +
                "s Y 2709736 42 + 57772954 CGCTTGCAGTGAGCTGAGATTGCGCCACTGCACTCCAGCCTG\n" +
                "s 1       0 42 +       42 CGCTTGCAGTGAGCTGAGATTGCGCCACTGCACTCCAGCCTG\n" +
                "\n" +
                "a score=41\n" +
                "s Y 2905156 41 + 57772954 GCTTGCAGTGAGCTGAGATTGCGCCACTGCACTCCAGCCTG\n" +
                "s 1       1 41 +       42 GCTTGCAGTGAGCTGAGATTGCGCCACTGCACTCCAGCCTG\n" +
                "\n" +
                "a score=41\n" +
                "s Y 3923174 41 + 57772954 GCTTGCAGTGAGCTGAGATTGCGCCACTGCACTCCAGCCTG\n" +
                "s 1       1 41 +       42 GCTTGCAGTGAGCTGAGATTGCGCCACTGCACTCCAGCCTG\n" +
                "\n" +
                "a score=40\n" +
                "s Y 5043226 42 + 57772954 CAGGCTGGAGTGCAGTGGCGCAATCTCGGCTCACTGCAAGCG\n" +
                "s 1       0 42 -       42 CAGGCTGGAGTGCAGTGGCGCAATCTCAGCTCACTGCAAGCG\n" +
                "\n" +
                "# CPU time: 0.28 seconds";
    }

    private String getMafInputWithQuality() {
        return "# fraglen=111.0 sdev=65.9758 disjoint=0.01 genome=2864785220\n" +
                "a score=540 mismap=0\n" +
                "s 4         17885592 90 + 191154276 CCACGTTTTTTCCTGGGCTGCTTACTGTCTTTTCGATCTAATCCATCTTCAGTATTTTCAGAGGTTCCATCAACTGTTCCATTTTTAGAA\n" +
                "s 5130000/1        7 90 -        97 CCACGTTTTTTCCTGGGCTGCTTACTGTCTTTTCGATCTAATCCATCTTCAGTATTTTCAGAGGTTCCATCAACTGTTCCATTTTTAGAA\n" +
                "q 5130000/1                         ^^^^^^^```_`ccccccccbc_`hhhchcchhhddd`hhdbc`bdddhhhdhchdhdhhhhhhdhdchdde`hhhdadhcccccccccc\n" +
                "\n" +
                "a score=558 mismap=0\n" +
                "s 4         17885595 93 + 191154276 CGTTTTTTCCTGGGCTGCTTACTGTCTTTTCGATCTAATCCATCTTCAGTATTTTCAGAGGTTCCATCAACTGTTCCATTTTTAGAATTAAGA\n" +
                "s 5130000/2        0 93 +        97 CGTTTTTTCCTGGGCTGCTTACTGTCTTTTCGATCTAATCCATCTTCAGTATTTTCAGAGGTTCCATCAACTGTTCCATTTTTAGAATTAAGA\n" +
                "q 5130000/2                         cccccccccchdadhhh`eddhcdhdhhhhhhdhdhchdhhhdddb`cbdhh`dddhhhcchchhh`_cbcccccccc`_```^^^^^^^^^`";
    }

    @Test(expected = NoSuchElementException.class)
    public void testReadQualities() throws IOException {
        final String mafInput = getMafInputWithQuality();
        final LastParser reader = new LastParser(new StringReader(mafInput));
        assertTrue(reader.hasNext());
        reader.next();
        assertEquals(540f, reader.getScore(), EPSILON);
        int sequenceLength = reader.getAlignedSequences().get(1).alignment.length();
        ByteArrayList qualityScores = reader.getAlignedSequences().get(1).getQualityScores();
        assertEquals(61, qualityScores.getByte(0));
        assertEquals(66, qualityScores.getByte(sequenceLength-1));

        assertTrue(reader.hasNext());
        reader.next();
        assertEquals(558f, reader.getScore(), EPSILON);
        qualityScores = reader.getAlignedSequences().get(1).getQualityScores();
        assertEquals(66, qualityScores.getByte(0));
        assertEquals(61, qualityScores.getByte(sequenceLength-1));

        assertFalse(reader.hasNext());

        // next should throw an exception if hasNext is false
        reader.next();

    }


}
