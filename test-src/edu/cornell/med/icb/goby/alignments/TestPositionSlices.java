/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 *  This file is part of the Goby IO API.
 *
 *     The Goby IO API is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     The Goby IO API is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with the Goby IO API.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.cornell.med.icb.goby.alignments;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.junit.BeforeClass;
import org.junit.AfterClass;
import org.junit.Test;
import static org.junit.Assert.*;
import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.io.File;

/**
 * @author Fabien Campagne
 *         Date: Dec 10, 2010
 *         Time: 2:08:43 PM
 */
public class TestPositionSlices {

    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(TestSkipTo.class);
    private static final String BASE_TEST_DIR = "test-results/alignments-position-slices";

    @BeforeClass
    public static void initializeTestDirectory() throws IOException {
        if (LOG.isDebugEnabled()) {
            LOG.debug("Creating base test directory: " + BASE_TEST_DIR);
        }
        FileUtils.forceMkdir(new File(BASE_TEST_DIR));
    }

    @AfterClass
    public static void cleanupTestDirectory() throws IOException {
        if (LOG.isDebugEnabled()) {
            LOG.debug("Deleting base test directory: " + BASE_TEST_DIR);
        }
        FileUtils.forceDeleteOnExit(new File(BASE_TEST_DIR));
    }

    @Test
    public void testAllEntries() throws IOException {
        final String basename = "align-position-slices-0";
        buildAlignment(basename);

        {  // check that we can read everything when the boundaries include everything:
            final AlignmentReader reader =
                    new AlignmentReaderImpl(FilenameUtils.concat(BASE_TEST_DIR, basename),
                            0, 400, 4, 500);

            check(reader, 1, 12);
            check(reader, 1, 13);
            check(reader, 1, 13);
            check(reader, 1, 13);
            check(reader, 2, 123);
            check(reader, 2, 300);
            check(reader, 2, 300);
            reader.close();
        }
    }

    @Test
    public void testOneSlice() throws IOException {
        final String basename = "align-position-slices-1";
        buildAlignment(basename);


        {// check that we can read everything only 1 12-13:
            final AlignmentReader reader =
                    new AlignmentReaderImpl(FilenameUtils.concat(BASE_TEST_DIR, basename),
                            0, 0, 1, 13);

            check(reader, 1, 12);
            check(reader, 1, 13);
            check(reader, 1, 13);
            check(reader, 1, 13);
            assertFalse(reader.hasNext());
            reader.close();
        }
    }

    @Test
    public void testTwoSlice() throws IOException {
        final String basename = "align-position-slices-2";
        buildAlignment(basename);


        {// check that we can read everything only 1 12-13:
            final AlignmentReader reader =
                    new AlignmentReaderImpl(FilenameUtils.concat(BASE_TEST_DIR, basename),
                            0, 0, 1, 12);

            check(reader, 1, 12);


            assertFalse(reader.hasNext());
            reader.close();
        }
    }

    @Test
    public void testThreeSlice() throws IOException {
        final String basename = "align-position-slices-3";
        buildAlignment(basename);


        {// check that we can read everything only between 1 13 and 2 123:
            final AlignmentReader reader =
                    new AlignmentReaderImpl(FilenameUtils.concat(BASE_TEST_DIR, basename),
                            1, 13, 2, 123);

            check(reader, 1, 13);
            check(reader, 1, 13);
            check(reader, 1, 13);
            check(reader, 2, 123);
            assertFalse(reader.hasNext());
            reader.close();
        }
    }

    @Test
    public void testFourSlice() throws IOException {
        final String basename = "align-position-slices-4";
        buildAlignment(basename);


        {// check that we can read everything only between 1 13 and 1 13:
            final AlignmentReader reader =
                    new AlignmentReaderImpl(FilenameUtils.concat(BASE_TEST_DIR, basename),
                            1, 13, 1, 13);

            check(reader, 1, 13);
            check(reader, 1, 13);
            check(reader, 1, 13);

            assertFalse(reader.hasNext());
            reader.close();
        }
    }

    @Test
    public void testBeforeSlice() throws IOException {
        final String basename = "align-position-slices-5";
        buildAlignment(basename);


        {// check that we can read skipTo before the start of the slice and only get entries that are within the slice
            // boundaries.
            final AlignmentReader reader =
                    new AlignmentReaderImpl(FilenameUtils.concat(BASE_TEST_DIR, basename),
                            2, 300, 2, 300);
            // (1,0) occurs before the slice boundaries (1,13)-(1,13)
            final Alignments.AlignmentEntry alignmentEntry = reader.skipTo(2, 0);

            assertNotNull(alignmentEntry);
            assertEquals(2, alignmentEntry.getTargetIndex());
            assertEquals(300, alignmentEntry.getPosition());
            check(reader, 2, 300);
            

            assertFalse(reader.hasNext());
            reader.close();
        }
    }

    private void buildAlignment(String basename) throws IOException {
        final AlignmentWriter writer =
                new AlignmentWriter(FilenameUtils.concat(BASE_TEST_DIR, basename));
        writer.setNumAlignmentEntriesPerChunk(2);

        final int numTargets = 3;
        final int[] targetLengths = new int[numTargets];

        for (int referenceIndex = 0; referenceIndex < numTargets; referenceIndex++) {
            targetLengths[referenceIndex] = 1000;
        }
        writer.setTargetLengths(targetLengths);
        // we write this alignment sorted:

        writer.setSorted(true);
        // chunk 1:
        writer.setAlignmentEntry(0, 1, 12, 30, false);
        writer.appendEntry();

        writer.setAlignmentEntry(0, 1, 13, 1, false);
        writer.appendEntry();
        // chunk 2:
        writer.setAlignmentEntry(0, 1, 13, 2, false);
        writer.appendEntry();

        writer.setAlignmentEntry(0, 1, 13, 3, false);
        writer.appendEntry();
        // chunk 3:
        writer.setAlignmentEntry(0, 2, 123, 30, false);
        writer.appendEntry();
        writer.setAlignmentEntry(0, 2, 300, 1, false);
        writer.appendEntry();
        // chunk 4:
        writer.setAlignmentEntry(0, 2, 300, 2, false);
        writer.appendEntry();
       // TODO remove the last entry and figure out why TestPositionSlices.testBeforeSlice cannot find the second/last (2,300) entry!
        writer.setAlignmentEntry(0, 3, 300, 2, false);
              writer.appendEntry();

        writer.close();
        writer.printStats(System.out);
    }

    private void check(AlignmentReader reader, int referenceIndex, int referencePosition) {
        assertTrue("Reader must have another entry", reader.hasNext());
        final Alignments.AlignmentEntry a = reader.next();
        assertNotNull(a);
        assertEquals(referenceIndex, a.getTargetIndex());
        assertEquals(referencePosition, a.getPosition());
    }
}
