/*
 * Copyright (C) 2009-2012 Institute for Computational Biomedicine,
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

import edu.cornell.med.icb.goby.alignments.perms.PermutationReader;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

import static org.junit.Assert.*;

/**
 * Check  that Goby can concatenate alignments that have query index permutations.
 * @author Fabien Campagne
 *         Date: May 20, 2009
 *         Time: 6:33:41 PM
 */
public class TestConcatAlignmentReaderWithPerms {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(TestConcatAlignmentReaderWithPerms.class);

    private static final String BASE_TEST_DIR = "test-results/alignments/concat-perm/";

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

    @Test
    public void testLoadTwoAdjustFalse() throws IOException {
        final int count;

        final ConcatAlignmentReader concatReader = new ConcatAlignmentReader(false, outputBasename2, outputBasename1);
        count = countAlignmentEntries(concatReader);
        assertEquals(count101 + count102, count);
        concatReader.readHeader();
        System.out.println(concatReader.getSmallestSplitQueryIndex());
        System.out.println(concatReader.getLargestSplitQueryIndex());
        // since both input alignments have permutations and distinct query indices, the total number of query index is
        // the sum of the two:
        assertEquals(numQueries101 + numQueries102, concatReader.getNumberOfQueries());

        concatReader.close();
        final ConcatAlignmentReader reader = new ConcatAlignmentReader(false, outputBasename2, outputBasename1);

        final String globalPermBasename = FilenameUtils.concat(BASE_TEST_DIR, "concatenated-permutation");

        concatReader.getConcatPerm().concatenate(globalPermBasename);
        PermutationReader permReader = new PermutationReader(globalPermBasename);

        while (concatReader.hasNext()) {
            final Alignments.AlignmentEntry alignmentEntry = reader.next();
            final int smallIndex = alignmentEntry.getQueryIndex();
            int queryIndex = permReader.getQueryIndex(smallIndex);
            final boolean score = alignmentEntry.getScore() == 30;
            if (score) {

                assertTrue(queryIndex >= 2000 && queryIndex <= 2010);
            }
            if (alignmentEntry.getScore() == 50) {

                assertTrue(queryIndex >= 1000 && queryIndex <= 1013);
            }

        }

        assertEquals(numTargets, concatReader.getNumberOfTargets());

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
            final AlignmentWriterImpl writer = new AlignmentWriterImpl(outputBasename1);
            writer.setNumAlignmentEntriesPerChunk(1000);
            writer.setAlignerName("first-aligner");
            writer.setAlignerVersion("version-first-aligner");
            final int numQuery = 10;
            int position = 100;
            final int score = 30;

            writer.setPermutation(true);
            for (int targetIndex = 0; targetIndex < numTargets; targetIndex++) {
                for (int queryIndex = 0; queryIndex < numQuery; queryIndex++) {
                    Alignments.AlignmentEntry.Builder newEntry = Alignments.AlignmentEntry.newBuilder();
                    newEntry.setQueryIndex(2000 + queryIndex);
                    newEntry.setTargetIndex(targetIndex);
                    newEntry.setScore((float) score);
                    newEntry.setPosition(position++);
                    newEntry.setMatchingReverseStrand(false);
                    newEntry.setMultiplicity(1);
                    newEntry.setQueryLength(constantQueryLength);
                    newEntry.setQueryIndexOccurrences(numTargets);
                    writer.appendEntry(newEntry.build());
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
            final AlignmentWriterImpl writer = new AlignmentWriterImpl(outputBasename2);
            writer.setAlignerName("second-aligner");
            writer.setAlignerVersion("version-second-aligner");
            writer.setNumAlignmentEntriesPerChunk(1000);
            writer.setPermutation(true);
            final int numQuery = 13;

            int position = 1;
            final int score = 50;
            for (int targetIndex = 0; targetIndex < numTargets; targetIndex++) {
                for (int queryIndex = 0; queryIndex < numQuery; queryIndex++) {

                    Alignments.AlignmentEntry.Builder newEntry = Alignments.AlignmentEntry.newBuilder();
                    newEntry.setQueryIndex(1000 + queryIndex);
                    newEntry.setTargetIndex(targetIndex);
                    newEntry.setScore((float) score);
                    newEntry.setPosition(position++);
                    newEntry.setMatchingReverseStrand(false);
                    newEntry.setMultiplicity(1);
                    newEntry.setQueryLength(constantQueryLength);
                    newEntry.setQueryIndexOccurrences(numTargets);
                    writer.appendEntry(newEntry.build());
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

}
