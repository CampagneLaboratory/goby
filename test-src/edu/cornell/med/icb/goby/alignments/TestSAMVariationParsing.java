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

import edu.cornell.med.icb.goby.modes.SAMToCompactMode;
import it.unimi.dsi.lang.MutableString;
import static org.junit.Assert.assertEquals;
import org.junit.Test;

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

        final MutableString referenceSequence = new MutableString();

        final String stringReadSequence = "TTTCCCACATTTCCCATCACCACTACTACGGATACAGAACGGGG";
        final MutableString expectedReferenceSequence = new MutableString("TTTCCCAAATTTCACATCACTACTACTACGGATACAGAACGGGG");
        testMdAttribute("44M", "7A5A6T23", stringReadSequence, referenceSequence);
        display(referenceSequence, stringReadSequence, expectedReferenceSequence);
        assertEquals("the reconstructed reference sequence should match expected.", expectedReferenceSequence, referenceSequence);


    }

    @Test
    //deletion in the read
    public void testMDStringParsing2() {

        final MutableString referenceSequence = new MutableString();

        final String stringReadSequence = "TTTCCCAAATTTCACATCACTACTACACGGATACAGAACGGGG";
        final MutableString expectedReferenceSequence = new MutableString("TTTCCCAAATTTCACATCACTACTACTACGGATACAGAACGGGG");
        testMdAttribute("26M1D17M", "26^T17", stringReadSequence, referenceSequence);
        display(referenceSequence, stringReadSequence, expectedReferenceSequence);
        assertEquals("the reconstructed reference sequence should match expected.", expectedReferenceSequence, referenceSequence);

    }

    private void display(final MutableString referenceSequence, final String stringReadSequence, final MutableString expectedReferenceSequence) {
        System.out.println(String.format("read       =%s\nexpected   =%s\nrecons. ref=%s", stringReadSequence,
                expectedReferenceSequence, referenceSequence));
        System.out.flush();
    }

    @Test
    // insertions in the read only
    public void testMDStringParsing3() {

        final MutableString referenceSequence = new MutableString();

        final String stringReadSequence = "TAAAACCTAAAAAAAAAAAAAAACCCC";
        final MutableString expectedReferenceSequence = new MutableString("TAAAA--TAAAAAAAAAAAAAAACCCC");
        testMdAttribute("5M2I20M", "25", stringReadSequence, referenceSequence);
        display(referenceSequence, stringReadSequence, expectedReferenceSequence);
        assertEquals("the reconstructed reference sequence should match expected.", expectedReferenceSequence, referenceSequence);

    }

    @Test
    // insertions with mutations
    public void testMDStringParsing4() {

        final MutableString referenceSequence = new MutableString();

        final String stringReadSequence = "TTTTGATGAAGTCTCTGTGTCCTGGGGCATCAATGATGGTCACA";
        final MutableString expectedReferenceSequence = new MutableString("TTTTGACGAAGTCTCTATGTCCT-GGGCATCAATGATGGTCACA");
        testMdAttribute("23M1I20M", "6C9A26", stringReadSequence, referenceSequence);
        display(referenceSequence, stringReadSequence, expectedReferenceSequence);
        assertEquals("the reconstructed reference sequence should match expected.", expectedReferenceSequence, referenceSequence);

    }

    @Test
    // insertions with mutations
    public void testMDStringParsing5() {

        final MutableString referenceSequence = new MutableString();

        final String stringReadSequence = "TTTCCCAAATTTCACATCACTACACTACGGATACAGAACGGGG";
        final MutableString expectedReferenceSequence = new MutableString("TTTCCCAAATTTCACATCACTACTACTACGGATACAGAACGGGG");
        testMdAttribute("23M1D20M", "23^T20", stringReadSequence, referenceSequence);
        display(referenceSequence, stringReadSequence, expectedReferenceSequence);
        assertEquals("the reconstructed reference sequence should match expected.", expectedReferenceSequence, referenceSequence);

    }

    @Test
    //cigar: 38M1I4M2D1M mdAttribute: 42^AA1
    public void testMDStringParsing6() {

        final MutableString referenceSequence = new MutableString();

        final String stringReadSequence =                                 "CCATGACCAACATAACTGTGGTGTCATGCATTTGGTATCTTTTT";
        final MutableString expectedReferenceSequence = new MutableString("CCATGACCAACATAACTGTGGTGTCATGCATTTGGTAT-TTTTAAT");
        testMdAttribute("38M1I4M2D1M", "42^AA1", stringReadSequence, referenceSequence);
        display(referenceSequence, stringReadSequence, expectedReferenceSequence);
        assertEquals("the reconstructed reference sequence should match expected.", expectedReferenceSequence, referenceSequence);

    }

    private void testMdAttribute(final String CIGAR, final String mdAttribute, final String stringReadSequence, final MutableString referenceSequence) {
        final MutableString readSequence = new MutableString(stringReadSequence);
        final MutableString readPostInsertions = new MutableString(stringReadSequence);

        SAMToCompactMode.produceReferenceSequence(CIGAR, mdAttribute, readSequence, readPostInsertions, referenceSequence);


    }


}
