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

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.junit.AfterClass;

import static org.junit.Assert.*;

import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: May 20, 2009
 *         Time: 6:33:41 PM
 */
public class TestConcatAlignmentReader {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(TestConcatAlignmentReader.class);

    private static final String BASE_TEST_DIR = "test-results/alignments/concat/";

    private int numEntriesIn101;
    private int numQueries101;
    private int numEntriesIn102;
    private int numQueries102;
    private final int numTargets = 5;
    private String outputBasename1;
    private String outputBasename2;
    private int count102;
    private int count101;
    private int constantQueryLength = 40;

    @Test
    public void testLoadTwo() throws IOException {
        final int count;

        final ConcatAlignmentReader concatReader = new ConcatAlignmentReader(outputBasename1, outputBasename2);
        count = countAlignmentEntries(concatReader);
        assertEquals(count101 + count102, count);
        concatReader.readHeader();

        assertEquals(numQueries101 + numQueries102, concatReader.getNumberOfQueries());
        assertEquals(numTargets, concatReader.getNumberOfTargets());

    }

    @Test
    public void testLoadTwoAlignerName() throws IOException {


        final ConcatAlignmentReader concatReader = new ConcatAlignmentReader(outputBasename1, outputBasename2);

        concatReader.readHeader();

        assertNotNull(concatReader.getAlignerName());
        assertNotNull(concatReader.getAlignerVersion());
        assertEquals("[first-aligner, second-aligner]", concatReader.getAlignerName());
        assertEquals("[version-first-aligner, version-second-aligner]", concatReader.getAlignerVersion());
    }

    @Test
    public void testLoadTwoAdjustFalse() throws IOException {
        final int count;

        final ConcatAlignmentReader concatReader = new ConcatAlignmentReader(false, outputBasename2, outputBasename1);
        count = countAlignmentEntries(concatReader);
        assertEquals(count101 + count102, count);
        concatReader.readHeader();
        System.out.println(concatReader.getSmallestSplitQueryIndex());
        System.out.println(concatReader.getLargestSplitQueryIndex());
        assertEquals(numQueries101 + numQueries102, concatReader.getNumberOfQueries());
        assertEquals(numTargets, concatReader.getNumberOfTargets());

    }

    @Test
    public void testQueryIndices() throws IOException {
        final ConcatAlignmentReader concatReader = new ConcatAlignmentReader(true, outputBasename1, outputBasename2);
        while (concatReader.hasNext()) {
            final Alignments.AlignmentEntry alignmentEntry = concatReader.next();

            if (alignmentEntry.getScore() == 50) {
                assertTrue("query index is too small: " + alignmentEntry.getQueryIndex(),
                        alignmentEntry.getQueryIndex() >= numQueries101);
            } else if (alignmentEntry.getScore() == 30) {
                assertTrue("query index is too big: " + alignmentEntry.getQueryIndex(),
                        alignmentEntry.getQueryIndex() < numQueries101);
            } else {
                fail("only scores possible are 30 and 50.");
            }

        }
    }

    @Test
    public void testQueryIndicesNoAdjustment() throws IOException {
        final ConcatAlignmentReader concatReader = new ConcatAlignmentReader(false, outputBasename1, outputBasename2);

        while (concatReader.hasNext()) {
            final Alignments.AlignmentEntry alignmentEntry = concatReader.next();

            if (alignmentEntry.getScore() == 50) {
                assertTrue(alignmentEntry.getQueryIndex() <= Math.max(numQueries101, numQueries102));
            } else if (alignmentEntry.getScore() == 30) {
                assertTrue(alignmentEntry.getQueryIndex() <= Math.max(numQueries101, numQueries102));
            } else {
                fail("only scores possible are 30 and 50.");
            }

        }
    }

    @Test
    public void testLoadNonAmbiguousOnly() throws IOException {

        // we now install a factory that removes entries whose reads match ambiguously:
        final ConcatAlignmentReader concatReader = new ConcatAlignmentReader(
                new NonAmbiguousAlignmentReaderFactory(),
                false, outputBasename1, outputBasename2);

        // the concat reads should now return only entries from basename2:
        int count = countAlignmentEntries(concatReader);
        assertEquals(count102, count);

    }


    private int countAlignmentEntries(final AbstractAlignmentReader reader) {
        int count = 0;

        while (reader.hasNext()) {
            final Alignments.AlignmentEntry alignmentEntry = reader.next();
            //System.out.println("found entry: " + alignmentEntry);
            assert alignmentEntry.hasPosition();
            count++;
        }
        return count;
    }

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

    @Before
    public void setUp() throws IOException {
        {
            outputBasename1 = FilenameUtils.concat(BASE_TEST_DIR, "concat-align-101");
            final AlignmentWriter writer = new AlignmentWriter(outputBasename1);
            writer.setNumAlignmentEntriesPerChunk(1000);
            writer.setAlignerName("first-aligner");
            writer.setAlignerVersion("version-first-aligner");
            final int numQuery = 10;
            int position = 100;
            final int score = 30;


            for (int targetIndex = 0; targetIndex < numTargets; targetIndex++) {
                for (int queryIndex = 0; queryIndex < numQuery; queryIndex++) {
                    writer.setAlignmentEntry(queryIndex, targetIndex, position++, score, false, constantQueryLength);
                    writer.appendEntry();
                    numEntriesIn101++;
                    count101++;
                }
            }
            numQueries101 = numQuery;

            writer.close();
            // reads in basename1 hit in 10 different places in the genome. They should be filtered when
            // removing ambiguous reads.
            AlignmentTooManyHitsWriter tmhWriter = new AlignmentTooManyHitsWriter(outputBasename1, 2);
            for (int queryIndex = 0; queryIndex < numQuery; queryIndex++) {
                tmhWriter.append(queryIndex, 10, 30);
            }
            tmhWriter.close();
        }
        {
            outputBasename2 = FilenameUtils.concat(BASE_TEST_DIR, "concat-align-102");
            final AlignmentWriter writer = new AlignmentWriter(outputBasename2);
            writer.setAlignerName("second-aligner");
            writer.setAlignerVersion("version-second-aligner");
            writer.setNumAlignmentEntriesPerChunk(1000);

            final int numQuery = 13;

            int position = 1;
            final int score = 50;
            for (int targetIndex = 0; targetIndex < numTargets; targetIndex++) {
                for (int queryIndex = 0; queryIndex < numQuery; queryIndex++) {
                    writer.setAlignmentEntry(queryIndex, targetIndex, position++, score, false, constantQueryLength);

                    writer.appendEntry();
                    numEntriesIn102++;
                    count102++;
                }
            }
            numQueries102 = numQuery;

            writer.close();
            // reads in basename 2 have just one hit, they are never ambiguous
            AlignmentTooManyHitsWriter tmhWriter = new AlignmentTooManyHitsWriter(outputBasename2, 2);
            for (int queryIndex = 0; queryIndex < numQuery; queryIndex++) {
                tmhWriter.append(queryIndex, 1, 30);
            }
            tmhWriter.close();
        }


    }

    @Test
    public void testLoadTwoFromFileURLs() throws IOException {
        final int count;

        final ConcatAlignmentReader concatReader = new ConcatAlignmentReader("file://./" + outputBasename1, "file://./" + outputBasename2);
        count = countAlignmentEntries(concatReader);
        assertEquals(count101 + count102, count);
        concatReader.readHeader();

        assertEquals(numQueries101 + numQueries102, concatReader.getNumberOfQueries());
        assertEquals(numTargets, concatReader.getNumberOfTargets());

    }

    // this test does URL connections to DropBox, so we give it up to 60 seconds before failing:
    @Test(timeout = 60000)
    public void testLoadTwoFromHttpURLs() throws IOException {
        final int count;
        /*
        final AlignmentReaderFactory alignmentReaderFactory,
                                 final boolean adjustQueryIndices,
                                 final int startReferenceIndex,
                                 final int startPosition,
                                 final int endReferenceIndex,
                                 final int endPosition,
                                 final String... basenames
         */
        // There are exactly 12 entries between position 33031693 and 33031798
        final ConcatAlignmentReader concatReader = new ConcatAlignmentReader(new DefaultAlignmentReaderFactory(),
                false, 21, 33031693, 21, 33031798,
                "http://dl.dropbox.com/u/357497/KHTFWNT-419-bis6-chr22-simulated-flat.entries",
                "http://dl.dropbox.com/u/357497/MCQPRWA-419-bis6-chr22-simulated-spikes.entries");
        count = countAlignmentEntries(concatReader);
        assertEquals(12, count);
        concatReader.readHeader();

        assertEquals(976, concatReader.getNumberOfQueries());
        assertEquals(84, concatReader.getNumberOfTargets());

    }
}
