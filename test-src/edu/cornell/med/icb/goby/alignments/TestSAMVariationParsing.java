/*
 * Copyright (C) 2010 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.alignments;

import org.junit.Test;
import static org.junit.Assert.assertEquals;
import it.unimi.dsi.lang.MutableString;

import java.util.regex.Pattern;
import java.util.regex.Matcher;

import edu.cornell.med.icb.goby.modes.SAMToCompactMode;

/**
 * Test the import of sequence variations from SAM bwa format.
 *
 * @author Fabien Campagne
 *         Date: Mar 10, 2010
 *         Time: 1:59:53 PM
 */
public class TestSAMVariationParsing {
    @Test
    // point mutations only
    public void testMDStringParsing1() {

        MutableString referenceSequence = new MutableString();

        String stringReadSequence = "TTTCCCACATTTCCCATCACCACTACTACGGATACAGAACGGGG";
        MutableString expectedReferenceSequence = new MutableString("TTTCCCAAATTTCACATCACTACTACTACGGATACAGAACGGGG");
        testMdAttribute("44M", "7A5A6T23", stringReadSequence, referenceSequence);
        display(referenceSequence, stringReadSequence, expectedReferenceSequence);
        assertEquals("the reconstructed reference sequence should match expected.", expectedReferenceSequence, referenceSequence);


    }

    @Test
    //deletion in the read
    public void testMDStringParsing2() {

        MutableString referenceSequence = new MutableString();

        String stringReadSequence =                                 "TTTCCCAAATTTCACATCACTACTACACGGATACAGAACGGGG";
        MutableString expectedReferenceSequence = new MutableString("TTTCCCAAATTTCACATCACTACTACTACGGATACAGAACGGGG");
        testMdAttribute("26M1D17M", "26^T17", stringReadSequence, referenceSequence);
        display(referenceSequence, stringReadSequence, expectedReferenceSequence);
        assertEquals("the reconstructed reference sequence should match expected.", expectedReferenceSequence, referenceSequence);

    }

    private void display(MutableString referenceSequence, String stringReadSequence, MutableString expectedReferenceSequence) {
        System.out.println(String.format("read       =%s\nexpected   =%s\nrecons. ref=%s", stringReadSequence,
                expectedReferenceSequence, referenceSequence));
        System.out.flush();
    }

    @Test
    // insertions in the read only
    public void testMDStringParsing3() {

        MutableString referenceSequence = new MutableString();

        String stringReadSequence = "TAAAACCTAAAAAAAAAAAAAAACCCC";
        MutableString expectedReferenceSequence = new MutableString("TAAAA--TAAAAAAAAAAAAAAACCCC");
        testMdAttribute("5M2I20M", "25", stringReadSequence, referenceSequence);
        display(referenceSequence, stringReadSequence, expectedReferenceSequence);
        assertEquals("the reconstructed reference sequence should match expected.", expectedReferenceSequence, referenceSequence);

    }

    @Test
    // insertions with mutations
    public void testMDStringParsing4() {

        MutableString referenceSequence = new MutableString();

        String stringReadSequence = "TTTTGATGAAGTCTCTGTGTCCTGGGGCATCAATGATGGTCACA";
        MutableString expectedReferenceSequence = new MutableString("TTTTGACGAAGTCTCTATGTCCT-GGGCATCAATGATGGTCACA");
        testMdAttribute("23M1I20M", "6C9A26", stringReadSequence, referenceSequence);
        display(referenceSequence, stringReadSequence, expectedReferenceSequence);
        assertEquals("the reconstructed reference sequence should match expected.", expectedReferenceSequence, referenceSequence);

    }

    @Test
       // insertions with mutations
       public void testMDStringParsing5() {

           MutableString referenceSequence = new MutableString();

           String stringReadSequence = "TTTCCCAAATTTCACATCACTACACTACGGATACAGAACGGGG";
           MutableString expectedReferenceSequence = new MutableString( "TTTCCCAAATTTCACATCACTACTACTACGGATACAGAACGGGG");
           testMdAttribute("23M1D20M", "23^T20", stringReadSequence, referenceSequence);
           display(referenceSequence, stringReadSequence, expectedReferenceSequence);
           assertEquals("the reconstructed reference sequence should match expected.", expectedReferenceSequence, referenceSequence);

       }

    private void testMdAttribute(String CIGAR, String mdAttribute, String stringReadSequence, MutableString referenceSequence) {
        MutableString readSequence = new MutableString(stringReadSequence);
        MutableString readPostInsertions = new MutableString(stringReadSequence);

        SAMToCompactMode.produceReferenceSequence(CIGAR, mdAttribute, readSequence, readPostInsertions, referenceSequence);


    }


}
