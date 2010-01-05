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

package edu.cornell.med.icb.goby.reads;

import org.apache.commons.io.input.NullInputStream;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;

/**
 * Validates the functionality of the {@link edu.cornell.med.icb.goby.reads.ReadsReader}.
 */
public class TestReadsReader {
    /**
     * Make sure that an empty file is handled properly.
     */
    @Test
    public void emptyFile() {
        // A zero length stream
        final ReadsReader reader = new ReadsReader(new NullInputStream(0));
        assertFalse("There should be no reads in a zero-length file", reader.hasNext());
    }

    /**
     * Make sure that a file with no entries is handled properly.
     */
    @Test
    public void noEntriesInFile() {
        // A simulated empty compact reads file
        final ReadsReader reader =
                new ReadsReader(new NullInputStream(MessageChunksWriter.DELIMITER_LENGTH + 4));
        assertFalse("There should be no reads in this file", reader.hasNext());
    }

    /**
     * Test that uses just the first byte of a reads file.  It should still return at least
     * one read.
     * @throws IOException if the file cannot be read
     */
    @Test
    public void wholeFile() throws IOException {
        final ReadsReader reader = new ReadsReader(0,
                new File("test-data/compact-reads/s_1_sequence_short.compact-reads").length(),
                "test-data/compact-reads/s_1_sequence_short.compact-reads");
        assertTrue("There should be at least one read in this file", reader.hasNext());
        int numberOfEntries = 0;
        for (final Reads.ReadEntry entry : reader) {
            assertNotNull("Entry should not be null: " + numberOfEntries, entry);
            numberOfEntries++;
        }
        assertEquals("There should be 73 entries in the test file", 73, numberOfEntries);
    }

    /**
     * Test that uses just the first byte of a reads file.  It should still return at least
     * one read.
     * @throws IOException if the file cannot be read
     */
    @Test
    public void firstByte() throws IOException {
        final ReadsReader reader =
                new ReadsReader(0, 1, "test-data/compact-reads/s_1_sequence_short.compact-reads");
        assertTrue("There should be at least one read in this file", reader.hasNext());
        int numberOfEntries = 0;
        for (final Reads.ReadEntry entry : reader) {
            assertNotNull("Entry should not be null: " + numberOfEntries, entry);
            numberOfEntries++;
        }
        assertEquals("There should be 73 entries in the test file", 73, numberOfEntries);
    }

    /**
     * Test that starts reading from the second byte of a reads file.  In this case it should
     * report zero entries since we started from past the initial header.  The file
     * used for testing has all the reads in a single chunk.
     * @throws IOException if the file cannot be read
     */
    @Test
    public void secondByte() throws IOException {
        final ReadsReader reader = new ReadsReader(1, 2,
                "test-data/compact-reads/s_1_sequence_short.compact-reads");
        assertFalse("There should be no reads in this file", reader.hasNext());
    }

    /**
     * Test that a range outside of the total file length doesn't throw an error but
     * more importantly, that it doesn't return any data.
     * @throws IOException if the file cannot be read
     */
    @Test
    public void outOfBounds() throws IOException {
        final ReadsReader reader = new ReadsReader(
                new File("test-data/compact-reads/s_1_sequence_short.compact-reads").length() + 10,
                Long.MAX_VALUE, "test-data/compact-reads/s_1_sequence_short.compact-reads");
        assertFalse("There should be no reads in this file", reader.hasNext());
    }

    /**
     * Test that a range with zero length returns no results.
     * @throws IOException if the file cannot be read
     */
    @Test
    public void zeroStartZeroEnd() throws IOException {
        final ReadsReader reader =
                new ReadsReader(0, 0, "test-data/compact-reads/s_1_sequence_short.compact-reads");
        assertFalse("There should be no reads in this file", reader.hasNext());
    }

    /**
     * Test that a range with zero length returns no results.
     * @throws IOException if the file cannot be read
     */
    @Test
    public void zeroLengthRange() throws IOException {
        final ReadsReader reader =
                new ReadsReader(0, 0, "test-data/compact-reads/s_1_sequence_short.compact-reads");
        assertFalse("There should be no reads in this file", reader.hasNext());
    }

    /**
     * Read from a file that only has a single sequence per chunk.
     * @throws IOException if the file cannot be read
     */
    @Test
    public void oneSequencePerChunk() throws IOException {
        final ReadsReader reader = new ReadsReader(0, 1,
                "test-data/compact-reads/s_1_sequence_short_1_per_chunk.compact-reads");
        assertTrue("There should be at least one read in this file", reader.hasNext());
        final Reads.ReadEntry entry = reader.next();
        assertNotNull("Entry should not be null", entry);
        assertEquals("Reader returned the wrong sequence string",
                "CTCATGTTCATACACCTNTCCCCCATTCTCCTCCT",
                entry.getSequence().toString(Charset.defaultCharset().name()));
        assertFalse("There should be no other reads in this file", reader.hasNext());
    }

    /**
     * Validate functionality of {@link ReadsReader#next()}.
     * @throws IOException if the file cannot be read
     */
    @Test
    public void firstAndSecondEntry() throws IOException {
        final ReadsReader reader = new ReadsReader(
                "test-data/compact-reads/s_1_sequence_short.compact-reads");
        assertTrue("There should be at least one read in this file", reader.hasNext());
        final Reads.ReadEntry entry = reader.next();
        assertNotNull("Entry should not be null", entry);
        assertEquals("Reader returned the wrong sequence string",
                "CTCATGTTCATACACCTNTCCCCCATTCTCCTCCT",
                entry.getSequence().toString(Charset.defaultCharset().name()));

        assertTrue("There should be at least two reads in this file", reader.hasNext());
        final Reads.ReadEntry entry2 = reader.next();
        assertNotNull("Entry should not be null", entry2);
        assertEquals("Reader returned the wrong sequence string",
                "GTAGGGGCCCCCTTTGCNTCCCTGGTGGCAACTGG",
                entry2.getSequence().toString(Charset.defaultCharset().name()));

        assertTrue("There should be at least three reads in this file", reader.hasNext());
    }

    /**
     * Check to make sure that calling {@link ReadsReader#next()} without first calling
     * {@link ReadsReader#hasNext()} produces the correct result.
     * @throws IOException if the file cannot be read
     */
    // @Test TODO: calling ReadsReader#next() before ReadsReader#hasNext() throws an exception
    public void nextWithNoHasFirst() throws IOException {
        final ReadsReader reader = new ReadsReader(
                "test-data/compact-reads/s_1_sequence_short_1_per_chunk.compact-reads");
        final Reads.ReadEntry entry = reader.next();
        assertNotNull("Entry should not be null", entry);
        assertEquals("Reader returned the wrong sequence string",
                "CTCATGTTCATACACCTNTCCCCCATTCTCCTCCT",
                entry.getSequence().toString(Charset.defaultCharset().name()));

        final Reads.ReadEntry entry2 = reader.next();
        assertNotNull("Entry should not be null", entry2);
        assertEquals("Reader returned the wrong sequence string",
                "GTAGGGGCCCCCTTTGCNTCCCTGGTGGCAACTGG",
                entry2.getSequence().toString(Charset.defaultCharset().name()));
    }

    /**
     * Check to make sure that calling {@link ReadsReader#next()} without first calling
     * {@link ReadsReader#hasNext()} produces the correct result.
     * @throws IOException if the file cannot be read
     */
    // @Test TODO: calling ReadsReader#next() before ReadsReader#hasNext() throws an exception
    public void nextWithNoHasFirstAndRange() throws IOException {
        final ReadsReader reader = new ReadsReader(0, 1,
                "test-data/compact-reads/s_1_sequence_short_1_per_chunk.compact-reads");
        final Reads.ReadEntry entry = reader.next();
        assertNotNull("Entry should not be null", entry);
        assertEquals("Reader returned the wrong sequence string",
                "CTCATGTTCATACACCTNTCCCCCATTCTCCTCCT",
                entry.getSequence().toString(Charset.defaultCharset().name()));
        assertFalse("There should be at least one read in this file", reader.hasNext());
    }
}
