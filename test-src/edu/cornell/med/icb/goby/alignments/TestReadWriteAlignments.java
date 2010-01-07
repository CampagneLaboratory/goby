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
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.junit.AfterClass;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: Apr 30, 2009
 *         Time: 6:15:52 PM
 */
public class TestReadWriteAlignments  {
    private static final Log LOG = LogFactory.getLog(TestReadWriteAlignments.class);
    private static final String BASE_TEST_DIR = "test-results/alignments";

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
    public void testReadWriteEntries() throws IOException {
        final AlignmentWriter writer = new AlignmentWriter(FilenameUtils.concat(BASE_TEST_DIR, "align-101"));
        writer.setNumAlignmentEntriesPerChunk(1000);
        final int numReads = 2000;
        final int numTargets = 10;
        int numExpected = 0;
        int position = 1;
        for (int referenceIndex = 0; referenceIndex < numTargets; referenceIndex++) {
            for (int queryIndex = 0; queryIndex < numReads; queryIndex++) {
                writer.setAlignmentEntry(queryIndex, referenceIndex, position++, 30, false);
                writer.appendEntry();
                numExpected++;
            }
        }

        writer.close();
        writer.printStats(System.out);

        int count = 0;
        final AlignmentReader reader = new AlignmentReader(FilenameUtils.concat(BASE_TEST_DIR, "align-101"));
        int maxQueryIndex = -1;
        int maxTargetIndex = -1;
        while (reader.hasNextAligmentEntry()) {
            final Alignments.AlignmentEntry alignmentEntry = reader.nextAlignmentEntry();
            assert alignmentEntry.hasPosition();
            maxQueryIndex = Math.max(maxQueryIndex, alignmentEntry.getQueryIndex());
            maxTargetIndex = Math.max(maxTargetIndex, alignmentEntry.getTargetIndex());
            count++;
        }
        assertEquals(numExpected, count);
        assertEquals("number of queries must match", numReads - 1, maxQueryIndex);
        assertEquals("number of targets must match", numTargets - 1, maxTargetIndex);
    }

    @Test
    public void testReadWriteHeader() throws IOException {
        AlignmentWriter writer = new AlignmentWriter(FilenameUtils.concat(BASE_TEST_DIR, "align-102"));
        final IndexedIdentifier queryIds = new IndexedIdentifier();
        assertNotNull(queryIds.keySet());
        writer.setQueryIdentifiers(queryIds);
        writer.setTargetIdentifiers(queryIds);
        writer.close();

        writer = new AlignmentWriter(FilenameUtils.concat(BASE_TEST_DIR, "align-102"));
        final IndexedIdentifier queryIds2 = new IndexedIdentifier();
        final IndexedIdentifier targetIds2 = new IndexedIdentifier();

        queryIds2.put(new MutableString("query:1"), 1);
        queryIds2.put(new MutableString("query:2"), 2);

        targetIds2.put(new MutableString("target:1"), 1);

        assertNotNull(queryIds.keySet());
        writer.setQueryIdentifiers(queryIds2);
        writer.setTargetIdentifiers(targetIds2);
        writer.close();


        AlignmentReader reader = new AlignmentReader(FilenameUtils.concat(BASE_TEST_DIR, "align-102"));
        reader.readHeader();
        assertEquals(1, reader.getQueryIdentifiers().getInt(new MutableString("query:1")));
        assertEquals(-1, reader.getTargetIdentifiers().getInt(new MutableString("query:1")));
        assertEquals(1, reader.getTargetIdentifiers().getInt(new MutableString("target:1")));
        assertEquals(-1, reader.getTargetIdentifiers().getInt(new MutableString("target:2")));

        writer = new AlignmentWriter(FilenameUtils.concat(BASE_TEST_DIR, "align-103"));

        assertNotNull(queryIds.keySet());
        writer.setQueryIdentifiers(queryIds2);
        writer.setTargetIdentifiers(targetIds2);
        writer.setQueryLengths(new int[]{0, 34, 84});
        writer.close();
        reader = new AlignmentReader(FilenameUtils.concat(BASE_TEST_DIR, "align-103"));
        reader.readHeader();
        final int[] lengths = reader.getQueryLength();
        assertEquals("at least three lengths must be available, since we entered a queryIndex=2",
                3, lengths.length);
        assertEquals(34, lengths[1]);
        assertEquals(84, lengths[2]);

        writer = new AlignmentWriter(FilenameUtils.concat(BASE_TEST_DIR, "align-103"));

        assertNotNull(queryIds.keySet());
        writer.setQueryIdentifiers(queryIds2);
        writer.setTargetIdentifiers(targetIds2);
        try {
            writer.setQueryLengths(new int[]{34, 84});
        } catch (AssertionError e) {  // TODO - change to a real exception -
            // success, not enough elements in length array to account for queryIndex=2
            return;
        }
        assertFalse("the length array is too short. This case must be detected.", true);
    }
}
