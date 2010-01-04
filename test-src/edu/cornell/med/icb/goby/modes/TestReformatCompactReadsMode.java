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

package edu.cornell.med.icb.goby.modes;

import edu.cornell.med.icb.goby.reads.Reads;
import edu.cornell.med.icb.goby.reads.ReadsReader;
import org.apache.commons.io.FileUtils;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.List;

/**
 * Validates the functionality of {@link edu.cornell.med.icb.goby.modes.ReformatCompactReadsMode}.
 */
public class TestReformatCompactReadsMode {
    /**
     * Validates that the reformat compact reads mode is capable of writing the same contents.
     * @throws IOException if there is a problem reading or writing to the files
     */
    @Test
    public void noChange() throws IOException {
        final ReformatCompactReadsMode reformat = new ReformatCompactReadsMode();

        final String inputFilename = "test-data/compact-reads/s_1_sequence_short.compact-reads";
        reformat.setInputFilenames(inputFilename);

        final String outputFilename = "test-results/reformat-test.compact-reads";
        reformat.setOutputFile(outputFilename);
        reformat.execute();

        assertTrue(FileUtils.contentEquals(new File(inputFilename), new File(outputFilename)));
    }

    /**
     * Validates that the reformat compact reads mode is capable of reformatting when
     * given positions at the extreme minimum and maximum values.
     * @throws IOException if there is a problem reading or writing to the files
     */
    @Test
    public void startAndEndAtExtremes() throws IOException {
        final ReformatCompactReadsMode reformat = new ReformatCompactReadsMode();

        final String inputFilename = "test-data/compact-reads/s_1_sequence_short.compact-reads";
        reformat.setInputFilenames(inputFilename);
        reformat.setStartPosition(0L);
        reformat.setEndPosition(Long.MAX_VALUE);
        final String outputFilename = "test-results/reformat-test-extremes.compact-reads";
        reformat.setOutputFile(outputFilename);
        reformat.execute();

        assertTrue(FileUtils.contentEquals(new File(inputFilename), new File(outputFilename)));
    }

    /**
     * Validates that the reformat compact reads mode is capable of reformatting when
     * given positions that exactly match the length of the file.
     * @throws IOException if there is a problem reading or writing to the files
     */
    @Test
    public void startAndEndAtExactLength() throws IOException {
        final ReformatCompactReadsMode reformat = new ReformatCompactReadsMode();

        final String inputFilename = "test-data/compact-reads/s_1_sequence_short.compact-reads";
        reformat.setInputFilenames(inputFilename);
        reformat.setStartPosition(0L);
        reformat.setEndPosition(new File(inputFilename).length() - 1);
        final String outputFilename = "test-results/reformat-test-start-end.compact-reads";
        reformat.setOutputFile(outputFilename);
        reformat.execute();

        assertTrue(FileUtils.contentEquals(new File(inputFilename), new File(outputFilename)));
    }

    /**
     * Validates that the reformat compact reads mode is capable of writing the same contents.
     * @throws IOException if there is a problem reading or writing to the files
     */
    @Test
    public void zeroEntries() throws IOException {
        final ReformatCompactReadsMode reformat = new ReformatCompactReadsMode();

        final String inputFilename = "test-data/compact-reads/s_1_sequence_short.compact-reads";
        reformat.setInputFilenames(inputFilename);
        reformat.setStartPosition(0L);
        reformat.setEndPosition(0L);
        final String outputFilename = "test-results/reformat-test-zero-zero.compact-reads";
        reformat.setOutputFile(outputFilename);
        reformat.execute();

        final File inputFile = new File(inputFilename);
        final File outputFile = new File(outputFilename);
        assertFalse("The reformatted file should not be the same as the original",
                FileUtils.contentEquals(inputFile, outputFile));

        final ReadsReader reader = new ReadsReader(FileUtils.openInputStream(outputFile));
        assertFalse("There should be no reads in this file", reader.hasNext());
    }

    /**
     * Validates that a file can be reformatted to change the chunk size.
     * @throws IOException if there is a problem reading or writing to the files
     */
    @Test
    public void reformatChunkSize() throws IOException {
        final ReformatCompactReadsMode reformat = new ReformatCompactReadsMode();

        final String inputFilename = "test-data/compact-reads/s_1_sequence_short.compact-reads";
        reformat.setInputFilenames(inputFilename);
        reformat.setSequencePerChunk(1);
        final String outputFilename = "test-results/reformat-test-chunk-size.compact-reads";
        reformat.setOutputFile(outputFilename);
        reformat.execute();

        final File inputFile = new File(inputFilename);
        final File outputFile = new File(outputFilename);
        assertFalse("The reformatted file should not be the same as the original",
                FileUtils.contentEquals(inputFile, outputFile));

        final ReadsReader inputReader = new ReadsReader(FileUtils.openInputStream(inputFile));
        assertTrue("There should be reads in this file", inputReader.hasNext());

        final List<Reads.ReadEntry> inputEntries = new ArrayList<Reads.ReadEntry>(73);
        for (final Reads.ReadEntry entry : inputReader) {
            inputEntries.add(entry);
        }

        final ReadsReader outputReader = new ReadsReader(FileUtils.openInputStream(outputFile));
        assertTrue("There should be reads in this file", outputReader.hasNext());

        final List<Reads.ReadEntry> outputEntries = new ArrayList<Reads.ReadEntry>(73);
        for (final Reads.ReadEntry entry : outputReader) {
            outputEntries.add(entry);
        }

        assertEquals("The entries of both files should be equal", inputEntries, outputEntries);
    }

    /**
     * Validates that a subset of a compact reads file can be written.
     * @throws IOException if there is a problem reading or writing to the files
     */
    @Test
    public void reformatStartOfCompactFile() throws IOException {
        final ReformatCompactReadsMode reformat = new ReformatCompactReadsMode();

        final String inputFilename =
                "test-data/compact-reads/s_1_sequence_short_1_per_chunk.compact-reads";
        reformat.setInputFilenames(inputFilename);
        reformat.setStartPosition(0);
        reformat.setEndPosition(10);
        final String outputFilename = "test-results/reformat-test-start.compact-reads";
        reformat.setOutputFile(outputFilename);
        reformat.execute();

        final File inputFile = new File(inputFilename);
        final File outputFile = new File(outputFilename);
        assertFalse("The reformatted file should not be the same as the original",
                FileUtils.contentEquals(inputFile, outputFile));

        final ReadsReader reader = new ReadsReader(FileUtils.openInputStream(outputFile));
        assertTrue("There should be reads in this file", reader.hasNext());
        final Reads.ReadEntry entry = reader.next();
        assertNotNull("Entry should not be null", entry);
        assertEquals("Reader returned the wrong sequence string",
                "CTCATGTTCATACACCTNTCCCCCATTCTCCTCCT",
                entry.getSequence().toString(Charset.defaultCharset().name()));
        assertFalse("There should be no other reads in this file", reader.hasNext());
    }

    /**
     * Validates that a subset of a compact reads file can be written.
     * @throws IOException if there is a problem reading or writing to the files
     */
    @Test
    public void reformatSubsetOfCompactFile() throws IOException {
        final ReformatCompactReadsMode reformat = new ReformatCompactReadsMode();

        final String inputFilename =
                "test-data/compact-reads/s_1_sequence_short_1_per_chunk.compact-reads";
        reformat.setInputFilenames(inputFilename);
        reformat.setStartPosition(11);
        final String outputFilename = "test-results/reformat-test-subset.compact-reads";
        reformat.setOutputFile(outputFilename);
        reformat.execute();

        final File inputFile = new File(inputFilename);
        final File outputFile = new File(outputFilename);
        assertFalse("The reformatted file should not be the same as the original",
                FileUtils.contentEquals(inputFile, outputFile));

        final ReadsReader reader = new ReadsReader(FileUtils.openInputStream(outputFile));
        assertTrue("There should be reads in this file", reader.hasNext());

        int numberOfEntries = 0;
        for (final Reads.ReadEntry entry : reader) {
            assertNotNull("Entry should not be null: " + numberOfEntries, entry);
            numberOfEntries++;
        }

        // we should have skipped the first entry
        assertEquals("There should be 72 entries in the test file", 72, numberOfEntries);
    }

    /**
     * Validates that setting a maximum read length will propertly exclude reads from
     * being written to the output.
     * @throws IOException if there is a problem reading or writing to the files
     */
    @Test
    public void excludeReadLengthsOf23() throws IOException {
        final ReformatCompactReadsMode reformat = new ReformatCompactReadsMode();

        final String inputFilename =
                "test-data/compact-reads/s_1_sequence_short_1_per_chunk.compact-reads";
        reformat.setInputFilenames(inputFilename);

        // there are no reads in the input file longer than 23
        reformat.setMaxReadLength(23);
        final String outputFilename =
                "test-results/reformat-test-exclude-read-lengths.compact-reads";
        reformat.setOutputFile(outputFilename);
        reformat.execute();

        final File inputFile = new File(inputFilename);
        final File outputFile = new File(outputFilename);
        assertFalse("The reformatted file should not be the same as the original",
                FileUtils.contentEquals(inputFile, outputFile));

        final ReadsReader reader = new ReadsReader(FileUtils.openInputStream(outputFile));
        assertFalse("There should be no reads in this file", reader.hasNext());
    }

    /**
     * Validates that setting a maximum read length will not exclude reads that are
     * within the limit from being written to the output.
     * @throws IOException if there is a problem reading or writing to the files
     */
    @Test
    public void excludeReadLengthsAt35() throws IOException {
        final ReformatCompactReadsMode reformat = new ReformatCompactReadsMode();

        final String inputFilename = "test-data/compact-reads/s_1_sequence_short.compact-reads";
        reformat.setInputFilenames(inputFilename);
        reformat.setMaxReadLength(35);

        final String outputFilename =
                "test-results/reformat-test-exclude-read-lengths-at-extrreme.compact-reads";
        reformat.setOutputFile(outputFilename);
        reformat.execute();

        assertTrue(FileUtils.contentEquals(new File(inputFilename), new File(outputFilename)));
    }

    /**
     * Validates that setting a read length trim value will not exclude reads that are
     * within the limit from being written to the output.
     * @throws IOException if there is a problem reading or writing to the files
     */
    @Test
    public void trimReadLengthsAt35() throws IOException {
        final ReformatCompactReadsMode reformat = new ReformatCompactReadsMode();

        final String inputFilename = "test-data/compact-reads/s_1_sequence_short.compact-reads";
        reformat.setInputFilenames(inputFilename);
        reformat.setTrimReadLength(35);

        final String outputFilename =
                "test-results/reformat-test-exclude-read-lengths-at-extrreme.compact-reads";
        reformat.setOutputFile(outputFilename);
        reformat.execute();

        assertTrue(FileUtils.contentEquals(new File(inputFilename), new File(outputFilename)));
    }

    /**
     * Validates that setting a read length trim value will write the trimmed reads to
     * the output properly.
     * @throws IOException if there is a problem reading or writing to the files
     */
    public void trimReadLengthsAt23() throws IOException {
        final ReformatCompactReadsMode reformat = new ReformatCompactReadsMode();

        final String inputFilename ="test-data/compact-reads/s_1_sequence_short.compact-reads";
        reformat.setInputFilenames(inputFilename);
        final String outputFilename = "test-results/reformat-test-start.compact-reads";
        reformat.setOutputFile(outputFilename);
        reformat.execute();
        reformat.setTrimReadLength(23);

        final File inputFile = new File(inputFilename);
        final File outputFile = new File(outputFilename);
        assertFalse("The reformatted file should not be the same as the original",
                FileUtils.contentEquals(inputFile, outputFile));

        final ReadsReader reader = new ReadsReader(FileUtils.openInputStream(outputFile));
        assertTrue("There should be reads in this file", reader.hasNext());

        int readCount = 0;
        for (final Reads.ReadEntry entry : reader) {
            if (readCount == 0) {
                assertEquals("Reader returned the wrong sequence string",
                        "CTCATGTTCATACACCTNTCCCCCATTCTCCTCCT".subSequence(0, 22),
                        entry.getSequence().toString(Charset.defaultCharset().name()));
            }
            readCount++;
            assertEquals("Entry ", readCount + " was not trimmed", entry.getReadLength());
        }
        assertEquals("There should have been 73 entries in the reformatted file", 73, readCount);
    }
}
