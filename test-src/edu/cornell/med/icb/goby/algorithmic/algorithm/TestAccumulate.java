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

package edu.cornell.med.icb.goby.algorithmic.algorithm;

import edu.cornell.med.icb.goby.algorithmic.data.Read;
import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.AlignmentWriter;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.counts.CountsReader;
import edu.cornell.med.icb.goby.counts.CountsWriter;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.junit.AfterClass;
import static org.junit.Assert.assertEquals;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

public class TestAccumulate {
    private static final Log LOG = LogFactory.getLog(TestAccumulate.class);
    private static final String BASE_TEST_DIR = "test-results/accumulate";
    private String testDir;
    private ComputeCount computeCount;

    @Test
    public void testBaseCount() {
        final ObjectList<Read> reads = new ObjectArrayList<Read>();
        final Read read1 = new Read(3, 3 + 5);
        final Read read2 = new Read(5, 5 + 5);
        final Read read3 = new Read(3, 3 + 4);
        reads.add(read1);
        reads.add(read2);
        reads.add(read3);

        computeCount.populate(reads);
        computeCount.accumulate();
//        assertEquals(0, computeCount.getValue(1, computeCount.startKeys, computeCount.starts));
//        assertEquals(2, computeCount.getValue(3, computeCount.startKeys, computeCount.starts));
        assertEquals(3, computeCount.totalCountOnReference());
    }

    @Test
    public void testBaseCount2() throws IOException {
        computeCount.startPopulating();
        computeCount.populate(3, 8);
        computeCount.populate(9, 10);
        computeCount.populate(5, 10);
        computeCount.populate(3, 7);
        computeCount.populate(8, 12);
        computeCount.populate(15, 20);
        computeCount.accumulate();

        computeCount.baseCount(); // final algorithm for base count without writer
    }

    @Test
    public void testBaseCountPipeLine() throws IOException {
        final String basename = FilenameUtils.concat(testDir, "align-reads");
        AlignmentWriter alignmentWriter = null;
        try {
            alignmentWriter = new AlignmentWriter(basename);
            Alignments.AlignmentEntry.Builder currentEntry = alignmentWriter.getAlignmentEntry();
            currentEntry.setScore(30);
            currentEntry.setPosition(5);
            currentEntry.setQueryIndex(1);
            currentEntry.setTargetIndex(1);
            currentEntry.setMatchingReverseStrand(false);
            currentEntry.setQueryAlignedLength(5);
            alignmentWriter.appendEntry();

            currentEntry = alignmentWriter.getAlignmentEntry();
            currentEntry.setScore(30);
            currentEntry.setPosition(3);
            currentEntry.setQueryIndex(1);
            currentEntry.setTargetIndex(1);
            currentEntry.setMatchingReverseStrand(false);
            currentEntry.setQueryAlignedLength(5);
            alignmentWriter.appendEntry();

            currentEntry = alignmentWriter.getAlignmentEntry();
            currentEntry.setScore(30);
            currentEntry.setPosition(3);
            currentEntry.setQueryIndex(1);
            currentEntry.setTargetIndex(1);
            currentEntry.setMatchingReverseStrand(false);
            currentEntry.setQueryAlignedLength(4);
            alignmentWriter.appendEntry();

            currentEntry = alignmentWriter.getAlignmentEntry();
            currentEntry.setScore(30);
            currentEntry.setPosition(8);
            currentEntry.setQueryIndex(1);
            currentEntry.setTargetIndex(1);
            currentEntry.setMatchingReverseStrand(false);
            currentEntry.setQueryAlignedLength(4);
            alignmentWriter.appendEntry();

            alignmentWriter.printStats(System.out);
        } finally {
            if (alignmentWriter != null) {
                try {
                    alignmentWriter.close();
                } catch (IOException e) { // NOPMD
                    // nothing to do - ignore
                }
            }
        }

        alignmentWriter.printStats(System.out);

        AlignmentReader alignmentReader = null;
        try {
            alignmentReader = new AlignmentReader(basename);
            computeCount.startPopulating();
            while (alignmentReader.hasNextAligmentEntry()) {
                final Alignments.AlignmentEntry alignmentEntry = alignmentReader.nextAlignmentEntry();
                final int startPosition = alignmentEntry.getPosition();
                final int alignmentLength = alignmentEntry.getQueryAlignedLength();
                System.out.println("start " + startPosition + " length " + alignmentLength);
                //shifted the ends populating by 1
                computeCount.populate(startPosition, startPosition + alignmentLength);
            }
        } finally {
            if (alignmentReader != null) {
                alignmentReader.close();
            }
        }
        final String countsFile = FilenameUtils.concat(testDir, "align-count");
        CountsWriter countsWriter = null;
        CountsReader countsReader = null;
        try {
            countsWriter = new CountsWriter(new FileOutputStream(countsFile));
            computeCount.accumulate();
            computeCount.baseCount(countsWriter);
            countsReader = new CountsReader(new FileInputStream(countsFile));
            final int[] exp = {0, 0, 0, 2, 2, 3, 3, 3, 3, 2, 2, 1, 1};
            int i = 0;
            while (countsReader.hasNextPosition()) {
                assertEquals(exp[i], countsReader.nextCountAtPosition());
                i++;
            }
        } finally {
            if (countsWriter != null) {
                try {
                    countsWriter.close();
                } catch (IOException e) { // NOPMD
                    // nothing to do - ignore
                }
            }

            if (countsReader != null) {
                try {
                    countsReader.close();
                } catch (IOException e) { // NOPMD
                    // nothing to do - ignore
                }
            }
        }
    }

     @Test
    public void testStartEnd() {
        final ComputeCount counter = new ComputeCount();
        counter.startPopulating();
        counter.populate(1, 10);
        counter.populate(1, 9);
        counter.populate(1, 8);
        counter.populate(1, 7);   // second read ends at 7
        counter.populate(1, 7);

        assertEquals(5, counter.getNumberOfReadsWithStartAt(1));
        assertEquals(0, counter.getNumberOfReadsWithStartAt(0));
        assertEquals(0, counter.getNumberOfReadsWithStartAt(2));

        assertEquals(1, counter.getNumberOfReadsWithEndAt(9));
        assertEquals(1, counter.getNumberOfReadsWithEndAt(8));
        assertEquals("two reads must end at position 7",2, counter.getNumberOfReadsWithEndAt(7));
        assertEquals(0, counter.getNumberOfReadsWithEndAt(6));
        assertEquals(0, counter.getNumberOfReadsWithEndAt(5));
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
        computeCount = new ComputeCount();
        final File tempDir = File.createTempFile("accumulate", "test", new File(BASE_TEST_DIR));
        FileUtils.deleteQuietly(tempDir);
        FileUtils.forceMkdir(tempDir);
        testDir = FilenameUtils.concat(BASE_TEST_DIR, tempDir.getName());

        if (LOG.isDebugEnabled()) {
            LOG.debug("Using test directory: " + testDir);
        }
    }
}
