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

import static org.junit.Assert.*;
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

}
