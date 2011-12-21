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

package edu.cornell.med.icb.goby.readers;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import org.junit.Test;

import java.io.IOException;

/**
 * Describe class here.
 *
 * @author Kevin Dorff
 */
public class TestFastXReader {
    @Test
    public void testFastA1() throws IOException {
        int entryNum = 0;
        for (final FastXEntry entry :
                new FastXReader("test-data/fastx-test-data/test-fasta-1.fa")) {
            if (entryNum == 0) {
                assertEquals("GATACCA", entry.getSequence().toString());
                assertEquals("GATACCA".length(), entry.getReadLength());
                assertEquals(0, entry.getQuality().length());
                assertEquals("GATACCA", entry.getEntrySansHeader().toString());
                assertEquals(">the header line\nGATACCA", entry.getEntry().toString());
                assertEquals('>', entry.getHeaderSymbol());
            }
            if (entryNum == 1) {
                assertEquals("TAGACCA", entry.getSequence().toString());
                assertEquals("TAGACCA".length(), entry.getReadLength());
                assertEquals(0, entry.getQuality().length());
                assertEquals("TAGACCA", entry.getEntrySansHeader().toString());
                assertEquals(">second entry\nTAGACCA", entry.getEntry().toString());
                assertEquals('>', entry.getHeaderSymbol());
            }

            assertTrue("Entry " + entryNum + " is not complete", entry.isEntryComplete());
            entryNum++;
        }
        assertEquals(2, entryNum);
    }

    @Test
    public void testFastA2Multiline() throws IOException {
        int entryNum = 0;
        for (final FastXEntry entry :
                new FastXReader("test-data/fastx-test-data/test-fasta-2.fa")) {
            if (entryNum == 0) {
                assertEquals("GATACCACATCA", entry.getSequence().toString());
                assertEquals("GATACCACATCA".length(), entry.getReadLength());
                assertEquals(0, entry.getQuality().length());
                assertEquals(">the header line\nGATACCA\nCATCA", entry.getEntry().toString());
                assertEquals("GATACCA\nCATCA", entry.getEntrySansHeader().toString());
                assertEquals('>', entry.getHeaderSymbol());
            }
            if (entryNum == 1) {
                assertEquals("TAGACCATAGG", entry.getSequence().toString());
                assertEquals("TAGACCATAGG".length(), entry.getReadLength());
                assertEquals(0, entry.getQuality().length());
                assertEquals(">second entry\nTAGACCA\nTAGG", entry.getEntry().toString());
                assertEquals("TAGACCA\nTAGG", entry.getEntrySansHeader().toString());
                assertEquals('>', entry.getHeaderSymbol());
            }
            assertTrue("Entry " + entryNum + " is not complete", entry.isEntryComplete());
            entryNum++;
        }
        assertEquals(2, entryNum);
    }

    @Test
    public void testFastq1() throws IOException {
        int entryNum = 0;
        for (final FastXEntry entry :
                new FastXReader("test-data/fastx-test-data/test-fastq-1.fq")) {
            if (entryNum == 0) {
                assertEquals("GATACCACATCA", entry.getSequence().toString());
                assertEquals("GATACCACATCA".length(), entry.getReadLength());
                assertEquals("123456789012", entry.getQuality().toString());
                assertEquals("123456789012".length(), entry.getQuality().length());
                assertEquals("@the header line\nGATACCA\nCATCA\n"
                        + "+the quality header\n1234\n567890\n12",
                        entry.getEntry().toString());
                assertEquals(
                        "GATACCA\nCATCA\n+the quality header\n1234\n567890\n12",
                        entry.getEntrySansHeader().toString());
                assertEquals('@', entry.getHeaderSymbol());
            }
            if (entryNum == 1) {
                assertEquals("TAGACCATAGG", entry.getSequence().toString());
                assertEquals("TAGACCATAGG".length(), entry.getReadLength());
                assertEquals("12345678902", entry.getQuality().toString());
                assertEquals("12345678902".length(), entry.getQuality().length());
                assertEquals("@second entry\nTAGACCA\nTAGG\n"
                        + "+the quality header2\n12345678902", entry.getEntry().toString());
                assertEquals("TAGACCA\nTAGG\n+the quality header2\n12345678902",
                        entry.getEntrySansHeader().toString());
                assertEquals('@', entry.getHeaderSymbol());
            }
            assertTrue("Entry " + entryNum + " is not complete", entry.isEntryComplete());
            entryNum++;
        }
        assertEquals(2, entryNum);
    }

    @Test
    public void testFastq2SymbolClash() throws IOException {
        int entryNum = 0;
        for (final FastXEntry entry :
                new FastXReader("test-data/fastx-test-data/test-fastq-2.fq")) {
            if (entryNum == 0) {
                assertEquals("GATACCACATCA", entry.getSequence().toString());
                assertEquals("GATACCACATCA".length(), entry.getReadLength());
                assertEquals("@234>67890+2", entry.getQuality().toString());
                assertEquals("the header line", entry.getEntryHeader().toString());
                assertEquals("@234>67890+2".length(), entry.getQuality().length());
                assertEquals("@the header line\nGATACCA\nCATCA\n"
                        + "+the quality header\n@234\n>67890\n+2",
                        entry.getEntry().toString());
                assertEquals(
                        "GATACCA\nCATCA\n+the quality header\n@234\n>67890\n+2",
                        entry.getEntrySansHeader().toString());
                assertEquals('@', entry.getHeaderSymbol());
            }
            if (entryNum == 1) {
                assertEquals("TAGACCATAGG", entry.getSequence().toString());
                assertEquals("TAGACCATAGG".length(), entry.getReadLength());
                assertEquals(">2345678902", entry.getQuality().toString());
                assertEquals(">2345678902".length(), entry.getQuality().length());
                assertEquals("@second entry\nTAGACCA\nTAGG\n"
                        + "+the quality header2\n>2345678902", entry.getEntry().toString());
                assertEquals("TAGACCA\nTAGG\n+the quality header2\n>2345678902",
                        entry.getEntrySansHeader().toString());
                assertEquals('@', entry.getHeaderSymbol());
            }
            assertTrue("Entry " + entryNum + " is not complete", entry.isEntryComplete());
            entryNum++;
        }
        assertEquals(2, entryNum);
    }

    @Test
    public void testCasavaFilter() throws IOException {
        final FastXReader reader = new FastXReader("test-data/fastx-test-data/sample_casava18.fq.gz");
        reader.setUseCasavaQualityFilter(true);
        int numRead = 0;
        for (final FastXEntry entry : reader) {
            numRead++;
        }
        assertEquals("Incorrect number of records read", 48062, numRead);
    }

    @Test
    public void testCasavaFilterDisabled() throws IOException {
        final FastXReader reader = new FastXReader("test-data/fastx-test-data/sample_casava18.fq.gz");
        int numRead = 0;
        for (final FastXEntry entry : reader) {
            numRead++;
        }
        assertEquals("Incorrect number of records read", 50000, numRead);
    }
}
