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

package edu.cornell.med.icb.goby.fastxReaders;

import static org.junit.Assert.assertEquals;
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
            entryNum++;
        }
        assertEquals(2, entryNum);
    }

}
