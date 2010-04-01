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

package edu.cornell.med.icb.goby.modes;

import edu.cornell.med.icb.goby.readers.FastXEntry;
import edu.cornell.med.icb.goby.readers.FastXReader;
import edu.cornell.med.icb.goby.reads.QualityEncoding;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.FileUtils;
import org.junit.After;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * Validates the functionality of {@link edu.cornell.med.icb.goby.modes.CompactToFastaMode}.
 */
public class TestCompactToFastaMode {
    /**
     * Expected sequence strings.
     */
    private final String[] expectedSequence = {
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            "TTTTGGAAGAGGTTCCACTGTAGGGTGTGGAAGGCAAG",
            "CTCTCTCGCCATTTAGATACCGTTTAACTTCACTGTTC",
            "CATTGAAGCCATGGATCCCAGCATCCTGAAGGGAGAGC",
            "CCGTGGGCAGCCTGGACGCACACACTGCGCTGGGGCTG"
    };

    /**
     * Expected quality scores. (encoded Illumina strings shown along with byte values)
     */
    private final byte[][] expectedQualityScores = {
            // XXS[SS_X[S\XRJZXXS_ZSSSMOM^QXVNMMZ_R_T
            {24, 24, 19, 27, 19, 19, 31, 24, 27, 19, 28, 24, 18, 10, 26, 24, 24, 19, 31,
             26, 19, 19, 19, 13, 15, 13, 30, 17, 24, 22, 14, 13, 13, 26, 31, 18, 31, 20},
            // S__bb`_bba^`W\VRGXU\X]a`_FXOZYFT]TZRPX
            {19, 31, 31, 34, 34, 32, 31, 34, 34, 33, 30, 32, 23, 28, 22, 18, 7, 24, 21,
             28, 24, 29, 33, 32, 31, 6, 24, 15, 26, 25, 6, 20, 29, 20, 26, 18, 16, 24},
            // aa`aaa]a`W^baa^aaaa_abab`aY]_``[D[b_\]
            {33, 33, 32, 33, 33, 33, 29, 33, 32, 23, 30, 34, 33, 33, 30, 33, 33, 33, 33,
             31, 33, 34, 33, 34, 32, 33, 25, 29, 31, 32, 32, 27, 4, 27, 34, 31, 28, 29},
            // abaaaaa`XP`_a_T^^T_aaXa`Z``a]]_ZZY_``^
            {33,  34,  33,  33,  33,  33,  33, 32, 24, 16, 32, 31, 33, 31, 20, 30, 30, 20,
             31, 33, 33, 24, 33, 32, 26, 32, 32, 33, 29, 29, 31, 26, 26, 25, 31, 32, 32, 30},
            // a`_P_`\]]`a[ZZ]VU_XXXOKTMSYBBBBBBBBBBB
            {33, 32, 31, 16, 31, 32, 28, 29, 29, 32, 33, 27, 26, 26, 29, 22, 21, 31,
             24, 24, 24, 15, 11, 20, 13, 19, 25, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}

    };

    /**
     * Temporary files for the tests so they can be deleted after the fact.
     */
    private final List<File> testFiles = new LinkedList<File>();

    /**
     * Remove any test files.
     * @throws IOException if any of the test files cannot be deleted
     */
    @After
    public void cleanup() throws IOException {
        final Iterator<File> files = testFiles.iterator();
        while (files.hasNext()) {
            final File file = files.next();
            if (file.exists()) {
                FileUtils.forceDeleteOnExit(file);
            }
            files.remove();
        }
    }

    /**
     * Create a temporary file for testing.  The file will be scheduled for deletion
     * when the test completes.
     * @param prefix The prefix string to be used in generating the file's name.
     * @param suffix The suffix string to be used in generating the file's name.
     * @return The newly-created empty file
     * @throws IOException if the file cannot be created
     */
    private File createTempFile(final String prefix, final String suffix) throws IOException {
        final File file = File.createTempFile(prefix, suffix);
        testFiles.add(file);
        return file;
    }

    /**
     * Test conversion of a compact file with no quality scores to FASTA format.
     * @throws IOException if the files cannot be read/written properly
     */
    @Test
    public void toFastaNoQuality() throws IOException {
        final String inputFilename = "test-data/compact-reads/s_1_sequence_short.compact-reads";
        final File fastaFile = createTempFile("toFasta", ".fasta");
        final String outputFilename = fastaFile.toString();
        final CompactToFastaMode compactToFastaMode = new CompactToFastaMode();
        compactToFastaMode.setInputFilename(inputFilename);
        compactToFastaMode.setOutputFilename(outputFilename);
        compactToFastaMode.setOutputFormat(CompactToFastaMode.OutputFormat.FASTA);
        compactToFastaMode.execute();

        final FastXReader reader = new FastXReader(outputFilename);
        assertEquals("File should be in FASTA format", "fa", reader.getFileType());

        int index = 0;
        for (final FastXEntry entry : reader) {
            assertEquals("Entry " + index + "symbol is not correct", '>', entry.getHeaderSymbol());
            assertEquals("Read length is not correct", 35, entry.getReadLength());
            final MutableString quality = entry.getQuality();
            assertNotNull("Quality string should never be null", quality);
            assertEquals("There should be no quality values", 0, quality.length());
            assertTrue("Entry " + index + " is not complete", entry.isEntryComplete());
            index++;
        }

        assertEquals(73, index);
        reader.close();
    }

    /**
     * Test conversion of a compact file with no quality scores to Illumina/FASTQ format.
     * @throws IOException if the files cannot be read/written properly
     */
    @Test
    public void toFastqIlluminaNoQuality() throws IOException {
        final String inputFilename = "test-data/compact-reads/s_1_sequence_short.compact-reads";
        final File fastqFile = createTempFile("llumina", ".fastq");
        final String outputFilename = fastqFile.toString();
        final CompactToFastaMode compactToFastaMode = new CompactToFastaMode();
        compactToFastaMode.setInputFilename(inputFilename);
        compactToFastaMode.setOutputFilename(outputFilename);
        compactToFastaMode.setOutputFormat(CompactToFastaMode.OutputFormat.FASTQ);
        compactToFastaMode.setQualityEncoding(QualityEncoding.ILLUMINA);
        compactToFastaMode.execute();

        final FastXReader reader = new FastXReader(outputFilename);
        assertEquals("File should be in FASTQ format", "fq", reader.getFileType());

        int index = 0;
        for (final FastXEntry entry : reader) {
            assertEquals("Entry " + index + "symbol is not correct", '@', entry.getHeaderSymbol());
            assertEquals("Read length is not correct", 35, entry.getReadLength());
            final MutableString quality = entry.getQuality();
            assertNotNull("Quality string should never be null", quality);
            assertEquals("There should be some quality values", 35, quality.length());
            // check quality scores
            for (int i = 0; i < 35; i++) {
                final char qualityCharacter = quality.charAt(i);
                assertEquals("Entry " + index + " has incorrect quality score at index " + i, 40,
                        QualityEncoding.ILLUMINA.asciiEncodingToQualityScore(qualityCharacter));
            }
            assertTrue("Entry " + index + " is not complete", entry.isEntryComplete());
            index++;
        }

        assertEquals(73, index);
        reader.close();
    }

    /**
     * Test conversion of a compact file with no quality scores to Sanger/FASTQ format.
     * @throws IOException if the files cannot be read/written properly
     */
    @Test
    public void toFastqSangerNoQuality() throws IOException {
        final String inputFilename = "test-data/compact-reads/s_1_sequence_short.compact-reads";
        final File fastqFile = createTempFile("Sanger", ".fastq");
        final String outputFilename = fastqFile.toString();
        final CompactToFastaMode compactToFastaMode = new CompactToFastaMode();
        compactToFastaMode.setInputFilename(inputFilename);
        compactToFastaMode.setOutputFilename(outputFilename);
        compactToFastaMode.setOutputFormat(CompactToFastaMode.OutputFormat.FASTQ);
        compactToFastaMode.setQualityEncoding(QualityEncoding.SANGER);
        compactToFastaMode.execute();

        final FastXReader reader = new FastXReader(outputFilename);
        assertEquals("File should be in FASTQ format", "fq", reader.getFileType());

        int index = 0;
        for (final FastXEntry entry : reader) {
            assertEquals("Entry " + index + "symbol is not correct", '@', entry.getHeaderSymbol());
            assertEquals("Read length is not correct", 35, entry.getReadLength());
            final MutableString quality = entry.getQuality();
            assertNotNull("Quality string should never be null", quality);
            assertEquals("There should be some quality values", 35, quality.length());
            // check quality scores
            for (int i = 0; i < 35; i++) {
                final char qualityCharacter = quality.charAt(i);
                assertEquals("Entry " + index + " has incorrect quality score at index " + i, 40,
                        QualityEncoding.SANGER.asciiEncodingToQualityScore(qualityCharacter));
            }
            assertTrue("Entry " + index + " is not complete", entry.isEntryComplete());
            index++;
        }

        assertEquals(73, index);
        reader.close();
    }

    /**
     * Test conversion of a compact file with quality scores to FASTA format.
     * @throws IOException if the files cannot be read/written properly
     */
    @Test
    public void toFastaWithQuality() throws IOException {
        final String inputFilename = "test-data/compact-reads/five-with-quality.compact-reads";
        final File fastaFile = createTempFile("toFasta", ".fasta");
        final String outputFilename = fastaFile.toString();
        final CompactToFastaMode compactToFastaMode = new CompactToFastaMode();
        compactToFastaMode.setInputFilename(inputFilename);
        compactToFastaMode.setOutputFilename(outputFilename);
        compactToFastaMode.setOutputFormat(CompactToFastaMode.OutputFormat.FASTA);
        compactToFastaMode.execute();

        final FastXReader reader = new FastXReader(outputFilename);
        assertEquals("File should be in FASTA format", "fa", reader.getFileType());

        int index = 0;
        for (final FastXEntry entry : reader) {
            assertEquals("Entry " + index + "symbol is not correct", '>', entry.getHeaderSymbol());
            assertEquals("Seqence for entry " + index + " is not correct",
                    expectedSequence[index], entry.getSequence().toString());
            final MutableString quality = entry.getQuality();
            assertNotNull("Quality string should never be null", quality);
            assertEquals("There should be no quality values", 0, quality.length());
            assertTrue("Entry " + index + " is not complete", entry.isEntryComplete());
            index++;
        }

        assertEquals(5, index);
        reader.close();
    }

    /**
     * Test conversion of a compact file with no quality scores to Illumina/FASTQ format.
     * @throws IOException if the files cannot be read/written properly
     */
    @Test
    public void toFastqIlluminaWithQuality() throws IOException {
        final String inputFilename = "test-data/compact-reads/five-with-quality.compact-reads";
        final File fastqFile = createTempFile("llumina", ".fastq");
        final String outputFilename = fastqFile.toString();
        final CompactToFastaMode compactToFastaMode = new CompactToFastaMode();
        compactToFastaMode.setInputFilename(inputFilename);
        compactToFastaMode.setOutputFilename(outputFilename);
        compactToFastaMode.setOutputFormat(CompactToFastaMode.OutputFormat.FASTQ);
        compactToFastaMode.setQualityEncoding(QualityEncoding.ILLUMINA);
        compactToFastaMode.execute();

        final FastXReader reader = new FastXReader(outputFilename);
        assertEquals("File should be in FASTQ format", "fq", reader.getFileType());

        int index = 0;
        for (final FastXEntry entry : reader) {
            assertEquals("Entry " + index + "symbol is not correct", '@', entry.getHeaderSymbol());
            assertEquals("Seqence for entry " + index + " is not correct",
                    expectedSequence[index], entry.getSequence().toString());
            final MutableString quality = entry.getQuality();
            assertNotNull("Quality string should never be null", quality);
            assertEquals("There should be some quality values", 38, quality.length());
            // check quality scores
            for (int i = 0; i < 38; i++) {
                final char qualityCharacter = quality.charAt(i);
                assertEquals("Entry " + index + " has incorrect quality score at index " + i,
                        expectedQualityScores[index][i],
                        QualityEncoding.ILLUMINA.asciiEncodingToQualityScore(qualityCharacter));
            }
            assertTrue("Entry " + index + " is not complete", entry.isEntryComplete());
            index++;
        }

        assertEquals(5, index);
        reader.close();
    }

    /**
     * Test conversion of a compact file with no quality scores to Sanger/FASTQ format.
     * @throws IOException if the files cannot be read/written properly
     */
    @Test
    public void toFastqSangerWithQuality() throws IOException {
        final String inputFilename = "test-data/compact-reads/five-with-quality.compact-reads";
        final File fastqFile = createTempFile("llumina", ".fastq");
        final String outputFilename = fastqFile.toString();
        final CompactToFastaMode compactToFastaMode = new CompactToFastaMode();
        compactToFastaMode.setInputFilename(inputFilename);
        compactToFastaMode.setOutputFilename(outputFilename);
        compactToFastaMode.setOutputFormat(CompactToFastaMode.OutputFormat.FASTQ);
        compactToFastaMode.setQualityEncoding(QualityEncoding.SANGER);
        compactToFastaMode.execute();

        final FastXReader reader = new FastXReader(outputFilename);
        assertEquals("File should be in FASTQ format", "fq", reader.getFileType());

        int index = 0;
        for (final FastXEntry entry : reader) {
            assertEquals("Entry " + index + "symbol is not correct", '@', entry.getHeaderSymbol());
            assertEquals("Seqence for entry " + index + " is not correct",
                    expectedSequence[index], entry.getSequence().toString());
            final MutableString quality = entry.getQuality();
            assertNotNull("Quality string should never be null", quality);
            assertEquals("There should be some quality values", 38, quality.length());
            // check quality scores
            for (int i = 0; i < 38; i++) {
                final char qualityCharacter = quality.charAt(i);
                assertEquals("Entry " + index + " has incorrect quality score at index " + i,
                        expectedQualityScores[index][i],
                        QualityEncoding.SANGER.asciiEncodingToQualityScore(qualityCharacter));
            }
            assertTrue("Entry " + index + " is not complete", entry.isEntryComplete());
            index++;
        }

        assertEquals(5, index);
        reader.close();
    }

    /**
     * Solexa encoding is not supported at this time.
     * @throws IOException if the files cannot be read/written properly
     */
    @Test(expected = UnsupportedOperationException.class)
    public void toSolexa() throws IOException {
        final CompactToFastaMode compactToFastaMode = new CompactToFastaMode();
        compactToFastaMode.setQualityEncoding(QualityEncoding.SOLEXA);
        compactToFastaMode.execute();
    }
}
