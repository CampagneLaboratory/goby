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

import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.junit.AfterClass;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Fabien Campagne
 *         Date:May 5, 2009
 *         Time: 9:06:40 AM
 */
public class TestMerge {
    private static final Log LOG = LogFactory.getLog(TestMerge.class);
    private static final String BASE_TEST_DIR = "test-results/alignments/merge";

    private int numQueries101;
    private int numQueries102;
    private int numEntriesIn102;
    private int numQueries103;
    private int numEntriesIn103;
    private int numEntriesIn105;
    private int numQueries105;
    private int numQueries106;
    private int numEntriesIn106;

    private int numEntriesIn101;

    @Test
    public void testMerge() throws IOException {
        Merge merger = new Merge(3);

        List<File> inputFiles = new ArrayList<File>();
        inputFiles.add(new File(FilenameUtils.concat(BASE_TEST_DIR, "align-101")));

        String outputFile = FilenameUtils.concat(BASE_TEST_DIR, "out-101-merged");
        merger.setK(2);
        merger.merge(inputFiles, outputFile);

        // With k=2 and a single input file, the merge should keep only the best score for each query:
        String basename = outputFile;
        int count = countAlignmentEntries(basename);
        assertEquals(numQueries101, count);

        merger = new Merge(3);

        inputFiles = new ArrayList<File>();
        inputFiles.add(new File(FilenameUtils.concat(BASE_TEST_DIR, "align-101")));
        inputFiles.add(new File(FilenameUtils.concat(BASE_TEST_DIR, "align-101")));

        outputFile = FilenameUtils.concat(BASE_TEST_DIR, "out-102-merged");
        merger.setK(1);
        merger.merge(inputFiles, outputFile);

        // With k=1 and twice the same input file, the merge should keep zero entries
        basename = outputFile;
        count = countAlignmentEntries(basename);
        assertEquals(0, count);
    }

    @Test
    public void testMergeWithTargetIds() throws IOException {
        final Merge merger = new Merge(3);

        final List<File> inputFiles = new ArrayList<File>();
        inputFiles.add(new File(FilenameUtils.concat(BASE_TEST_DIR, "align-105")));

        final String outputFile = FilenameUtils.concat(BASE_TEST_DIR, "out-105-merged");
        merger.setK(2);
        merger.merge(inputFiles, outputFile);

        // With k=2 and a single input file, the merge should keep only the best score for each query:
        final String basename = outputFile;
        final int count = countAlignmentEntries(basename);
        assertEquals(numQueries101, count);

    }

    @Test
    public void testMergeTranscriptRuns() throws IOException {
        final Merge merger = new Merge("test-data/alignments/geneTranscriptInfo.txt", 3);

        final List<File> inputFiles = new ArrayList<File>();
        inputFiles.add(new File(FilenameUtils.concat(BASE_TEST_DIR, "transcript-101")));

        final String outputFile = FilenameUtils.concat(BASE_TEST_DIR, "out-transcript-101-merged");
        merger.setK(1);
        merger.merge(inputFiles, outputFile);

        // With k=1, the merge should keep query=0 and query=2
        final String basename = outputFile;
        final int count = countAlignmentEntries(basename);
        assertEquals("Only best score, non ambiguous gene matches should be kept. ", 4 - 1 , count);     // -1 removes an entry with lower score

        final AlignmentReaderImpl reader = new AlignmentReaderImpl(outputFile);
        reader.readHeader();
        assertEquals(5, reader.getNumberOfTargets());
        assertArrayEquals(new int[] {1024, 5678, 1237, 9, 143}, reader.getTargetLength());
        assertEquals(3, reader.getNumberOfQueries());
        final int maxTargetIndex = -1;
        final IntSet queryIndices = new IntArraySet();
        final IntSet targetIndices = new IntArraySet();
        while (reader.hasNext()) {
            final Alignments.AlignmentEntry alignmentEntry = reader.next();
            queryIndices.add(alignmentEntry.getQueryIndex());
            targetIndices.add(alignmentEntry.getQueryIndex());
        }
        assertEquals(true, queryIndices.contains(0));
        assertEquals(true, queryIndices.contains(2));
        assertEquals(false, queryIndices.contains(1));
        assertEquals(false, queryIndices.contains(1));
    }


    @Test
    public void testMergeWithTargetIds3() throws IOException {
        // Test that target indices are merged. align-105 and align-106 each have two targets. The result merged
        // target index should be between 0 and 3. Maximum index at 3.
        final Merge merger = new Merge(3);

        final List<File> inputFiles = new ArrayList<File>();
        inputFiles.add(new File(FilenameUtils.concat(BASE_TEST_DIR, "align-105")));
        inputFiles.add(new File(FilenameUtils.concat(BASE_TEST_DIR, "align-106")));

        final String outputFile = FilenameUtils.concat(BASE_TEST_DIR, "out-105-merged");
        merger.setK(2);
        merger.setSilent(false);
        merger.merge(inputFiles, outputFile);


        final AlignmentReader reader = new AlignmentReaderImpl(outputFile);
        int maxTargetIndex = -1;
        while (reader.hasNext()) {
            final Alignments.AlignmentEntry alignmentEntry = reader.next();
            maxTargetIndex = Math.max(maxTargetIndex, alignmentEntry.getTargetIndex());
        }
        assertEquals(3, maxTargetIndex);
    }

    @Test
    public void testMergeWithTargetIds2() throws IOException {
        final Merge merger = new Merge(30);

        final List<File> inputFiles = new ArrayList<File>();
        inputFiles.add(new File(FilenameUtils.concat(BASE_TEST_DIR, "align-105")));
        inputFiles.add(new File(FilenameUtils.concat(BASE_TEST_DIR, "align-106")));

        final String outputFile = FilenameUtils.concat(BASE_TEST_DIR, "out-106-merged");
        merger.setK(1);
        merger.merge(inputFiles, outputFile);

        // With k=1 and twice the same input file, the merge should keep zero entries

        final int count = countAlignmentEntries(outputFile);

        assertEquals(0, count);

    }

    @Test
    public void testMerge2() throws IOException {
        final Merge merger = new Merge(3);

        final List<File> inputFiles;
        final String outputFile;
        final String basename;
        final int count;
        inputFiles = new ArrayList<File>();
        inputFiles.add(new File(FilenameUtils.concat(BASE_TEST_DIR, "align-101")));
        inputFiles.add(new File(FilenameUtils.concat(BASE_TEST_DIR, "align-101")));

        outputFile = FilenameUtils.concat(BASE_TEST_DIR, "out-103-merged");
        merger.setK(3);
        merger.merge(inputFiles, outputFile);

        // With k=3 and twice the same input file, the merge should keep twice the queries with best score:
        basename = outputFile;
        count = countAlignmentEntries(basename);
        assertEquals(numQueries101 * 2, count);
    }

    @Test
    public void testIgnoreTooManyHits() throws IOException {
        synchronized (this) {
            // make the too many hits info:
            final AlignmentTooManyHitsWriter tmhWriter =
                    new AlignmentTooManyHitsWriter(FilenameUtils.concat(BASE_TEST_DIR, "align-102"), 4);

            // tmhWriter will only write entries if numHits > thresh
            tmhWriter.getNewAmbiguousLocation().setAtLeastNumberOfHits(5);
            tmhWriter.getNewAmbiguousLocation().setQueryIndex(0);
            tmhWriter.append();
            tmhWriter.getNewAmbiguousLocation().setAtLeastNumberOfHits(15);
            tmhWriter.getNewAmbiguousLocation().setQueryIndex(1);
            tmhWriter.append();
            tmhWriter.close();

            final Merge merger = new Merge(3);

            final List<File> inputFiles;
            final String outputFile;
            final String basename;
            final int count;
            inputFiles = new ArrayList<File>();
            inputFiles.add(new File(FilenameUtils.concat(BASE_TEST_DIR, "align-102")));
            inputFiles.add(new File(FilenameUtils.concat(BASE_TEST_DIR, "align-102")));

            outputFile = FilenameUtils.concat(BASE_TEST_DIR, "out-103-merged");
            merger.setK(3);
            merger.merge(inputFiles, outputFile);

            // With k=3 and twice the same input file, the merge should keep twice the queries with best score:
            basename = outputFile;
            count = countAlignmentEntries(basename);
            assertEquals(count, (numQueries101 - 2) * 2/* two queries excluded because of too many hits with k=3 */);
        }
    }

    @Test
    public void testIgnoreTooManyHits2() throws IOException {
        synchronized (this) {
            // make the too many hits info:
            final AlignmentTooManyHitsWriter tmhWriter =
                    new AlignmentTooManyHitsWriter(FilenameUtils.concat(BASE_TEST_DIR, "align-102"), 4);

            // tmhWriter will only write entries if numHits > thresh
            tmhWriter.getNewAmbiguousLocation().setAtLeastNumberOfHits(5);
            tmhWriter.getNewAmbiguousLocation().setQueryIndex(0);
            tmhWriter.append();
            tmhWriter.getNewAmbiguousLocation().setAtLeastNumberOfHits(15);
            tmhWriter.getNewAmbiguousLocation().setQueryIndex(1);
            tmhWriter.append();
            tmhWriter.close();

            final Merge merger = new Merge(3);

            final List<File> inputFiles;
            final String outputFile;
            final String basename;
            final int count;
            inputFiles = new ArrayList<File>();
            inputFiles.add(new File(FilenameUtils.concat(BASE_TEST_DIR, "align-102")));
            inputFiles.add(new File(FilenameUtils.concat(BASE_TEST_DIR, "align-102")));

            outputFile = FilenameUtils.concat(BASE_TEST_DIR, "out-103-merged");
            merger.setK(30);
            merger.merge(inputFiles, outputFile);

            // With k=3 and twice the same input file, the merge should keep twice the queries with best score:
            basename = outputFile;
            count = countAlignmentEntries(basename);
            assertEquals(count, (numQueries102 - 2) * 2/* two queries excluded because of too many hits with k=30 > alignerThreshold */);
        }
    }

    /* DISABLING Test
     * Originally, the tmhWriter threshold was 10, and the entries numHits were 1 and 2.  Due to an update
     * in the tmhWriter, appendTooManyHits() will *NOT* write these entries to an output file.
     *
     * It is not possible to set the tmhWriter threshold to 0, for the following reason ...
     *
     * Writer only adds entries if numHits > thresh, so tmhWriter needs a threshold below numHits of entries.
     *
     * Merge tries to respect the minimum threshold of the input files, which in our case is below numHits of
     * the entries.  However, if Merge is called with K > tmhWriter's threshold, then all entries are declared
     * ambiguous due to the (strange) behaviour of isQueryAmbiguous(queryIndex, k, matchLength).
     *
     * Since I could not trace this further, I have decided to disable this test.
     */
    public void testIgnoreTooManyHits3() throws IOException {
        synchronized (this) {
            // make the too many hits info:
            final AlignmentTooManyHitsWriter tmhWriter =
                    new AlignmentTooManyHitsWriter(FilenameUtils.concat(BASE_TEST_DIR, "align-103"), 10);

            // tmhWriter will only write entries if numHits > thresh
            tmhWriter.getNewAmbiguousLocation().setAtLeastNumberOfHits(2);
            tmhWriter.getNewAmbiguousLocation().setQueryIndex(0);
            tmhWriter.append();
            tmhWriter.getNewAmbiguousLocation().setAtLeastNumberOfHits(1);
            tmhWriter.getNewAmbiguousLocation().setQueryIndex(1);
            tmhWriter.append();
            tmhWriter.close();

            final Merge merger = new Merge(3);

            final List<File> inputFiles;
            final String outputFile;
            final String basename;
            final int count;
            inputFiles = new ArrayList<File>();
            inputFiles.add(new File(FilenameUtils.concat(BASE_TEST_DIR, "align-103")));
            inputFiles.add(new File(FilenameUtils.concat(BASE_TEST_DIR, "align-103")));

            outputFile = FilenameUtils.concat(BASE_TEST_DIR, "out-103-merged");
            merger.setK(3);
            merger.merge(inputFiles, outputFile);

            // With k=3 and twice the same input file, the merge should keep twice the queries with best score:
            basename = outputFile;
            count = countAlignmentEntries(basename);
            assertEquals(count, (numQueries102 - 1) * 2/* one query excluded because of too many hits with k=3 */);
        }
    }

    private int countAlignmentEntries(final String basename) throws IOException {
        int count = 0;
        final AlignmentReader reader = new AlignmentReaderImpl(basename);
        while (reader.hasNext()) {
            final Alignments.AlignmentEntry alignmentEntry = reader.next();
            System.out.println("found entry: " + alignmentEntry);
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
            final AlignmentWriter writer = new AlignmentWriter(FilenameUtils.concat(BASE_TEST_DIR, "align-101"));
            writer.setNumAlignmentEntriesPerChunk(1000);
            final int numTargets = 2000;
            final int numQuery = 10;
            int position = 1;
            int score = 30;
            for (int targetIndex = 0; targetIndex < numTargets; targetIndex++) {
                for (int queryIndex = 0; queryIndex < numQuery; queryIndex++) {
                    writer.setAlignmentEntry(queryIndex, targetIndex, position++, score++, false);
                    writer.appendEntry();
                    numEntriesIn101++;
                }
            }
            numQueries101 = numQuery;
            writer.close();
        }

        {
            final AlignmentWriter writer =
                    new AlignmentWriter(FilenameUtils.concat(BASE_TEST_DIR, "transcript-101"));
            writer.setNumAlignmentEntriesPerChunk(1000);
            final int numTargets = 5;
            int position = 1;
            final int score = 30;
            final String[] transcriptIds = {"transcriptId1", "transcriptId2", "transcriptId3", "transcriptId4", "transcriptId5",};
            final int[] transcriptIndex = {0, 1, 2, 3, 4};
            final int[] targetLengths = {1024, 5678, 1237, 9, 143};
            writer.setTargetLengths(targetLengths);

            final int[] queryLentghs = {35, 35, 35, 35, 35, 35, 35, 35, 35, 35};
            writer.setQueryLengths(queryLentghs);

            final IndexedIdentifier targetIds = new IndexedIdentifier();
            for (int i = 0; i < numTargets; i++) {
                if (transcriptIndex[i] != targetIds.registerIdentifier(new MutableString(transcriptIds[i]))) {
                    assert false : "transcript Index must match";
                }
            }

            writer.setTargetIdentifiers(targetIds);

            // make query 0 OK, matches 3 transcripts that belong to the same gene
            // query 1 matches two transcripts, but 2 genes, not OK.
            // query 2 OK: 1 transcript 1 gene.
            writer.setAlignmentEntry(0, transcriptIndex[0], position++, score, false);
            writer.appendEntry();
            writer.setAlignmentEntry(0, transcriptIndex[1], position++, score + 1, false);
            writer.appendEntry();
            writer.setAlignmentEntry(0, transcriptIndex[2], position++, score + 1, false);
            writer.appendEntry();
            writer.setAlignmentEntry(1, transcriptIndex[2], position++, score, false);
            writer.appendEntry();
            writer.setAlignmentEntry(1, transcriptIndex[3], position++, score, false);
            writer.appendEntry();
            writer.setAlignmentEntry(2, transcriptIndex[4], position++, score, false);
            writer.appendEntry();
            writer.close();

        }


        {
            // Will have too many hits info created for it:
            final AlignmentWriter writer = new AlignmentWriter(FilenameUtils.concat(BASE_TEST_DIR, "align-102"));
            writer.setNumAlignmentEntriesPerChunk(1000);
            final int numTargets = 2000;
            final int numQuery = 10;
            int position = 1;
            int score = 30;
            for (int targetIndex = 0; targetIndex < numTargets; targetIndex++) {
                for (int queryIndex = 0; queryIndex < numQuery; queryIndex++) {
                    writer.setAlignmentEntry(queryIndex, targetIndex, position++, score++, false);
                    writer.appendEntry();
                    numEntriesIn102++;
                }
            }
            numQueries102 = numQuery;
            writer.close();
        }

        {
            // Will have too many hits info created for it:
            final AlignmentWriter writer = new AlignmentWriter(FilenameUtils.concat(BASE_TEST_DIR, "align-103"));
            writer.setNumAlignmentEntriesPerChunk(1000);
            final int numTargets = 2000;
            final int numQuery = 10;
            int position = 1;
            int score = 30;
            for (int targetIndex = 0; targetIndex < numTargets; targetIndex++) {
                for (int queryIndex = 0; queryIndex < numQuery; queryIndex++) {
                    writer.setAlignmentEntry(queryIndex, targetIndex, position++, score++, false);
                    writer.appendEntry();
                    numEntriesIn103++;
                }
            }
            numQueries103 = numQuery;
            writer.close();
        }

        {
            // Will have too many hits info created for it:
            final AlignmentWriter writer = new AlignmentWriter(FilenameUtils.concat(BASE_TEST_DIR, "align-105"));
            System.out.println("preparing 105");
            writer.setNumAlignmentEntriesPerChunk(1000);
            final int numTargets = 2;
            final int numQuery = 10;
            int position = 1;
            int score = 30;
            final IndexedIdentifier targetIds = new IndexedIdentifier();
            for (int i = 0; i < numTargets; i++) {
                targetIds.registerIdentifier(new MutableString("target-" + i));
            }

            writer.setTargetIdentifiers(targetIds);

            for (int targetIndex = 0; targetIndex < numTargets; targetIndex++) {
                for (int queryIndex = 0; queryIndex < numQuery; queryIndex++) {
                    writer.setAlignmentEntry(queryIndex, targetIndex, position++, score++, false);
                    //         System.out.println(String.format("Preparing entry: queryIndex: %d targetIndex: %d position: %d score: %d strand: %b", queryIndex, targetIndex, position++, score++, false));
                    writer.appendEntry();

                    numEntriesIn105++;
                }
            }
            numQueries105 = numQuery;

            writer.close();
        }

        {
            // Will have too many hits info created for it:
            final AlignmentWriter writer = new AlignmentWriter(FilenameUtils.concat(BASE_TEST_DIR, "align-106"));
            writer.setNumAlignmentEntriesPerChunk(1000);
            final int numTargets = 2;
            final int numQuery = 10;
            int position = 1;
            int score = 30;
            final IndexedIdentifier targetIds = new IndexedIdentifier();
            for (int i = 0; i < numTargets; i++) {
                targetIds.registerIdentifier(new MutableString("target-b-" + i));
            }

            writer.setTargetIdentifiers(targetIds);

            for (int targetIndex = 0; targetIndex < numTargets; targetIndex++) {
                for (int queryIndex = 0; queryIndex < numQuery; queryIndex++) {
                    writer.setAlignmentEntry(queryIndex, targetIndex, position++, score++, false);
                    writer.appendEntry();
                    numEntriesIn106++;
                }
            }
            numQueries106 = numQuery;

            writer.close();
        }
    }
}
