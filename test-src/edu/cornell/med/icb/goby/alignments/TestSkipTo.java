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

package edu.cornell.med.icb.goby.alignments;

import edu.cornell.med.icb.goby.compression.HybridChunkCodec1;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.*;

import static org.junit.Assert.*;

/**
 * @author Fabien Campagne
 *         Date: Apr 30, 2009
 *         Time: 6:15:52 PM
 */
public class TestSkipTo {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(TestSkipTo.class);
    private static final String BASE_TEST_DIR = "test-results/alignments-skip-to";
    private int numEntriesPerChunk = 2;
    private int constantQueryLength = 40;


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
    public void testFewSkips1() throws IOException {
        final String basename = "align-skip-to-1";
        final AlignmentWriterImpl writer =
                new AlignmentWriterImpl(FilenameUtils.concat(BASE_TEST_DIR, basename));
        writer.setNumAlignmentEntriesPerChunk(numEntriesPerChunk);

        final int numTargets = 3;
        final int[] targetLengths = new int[numTargets];

        for (int referenceIndex = 0; referenceIndex < numTargets; referenceIndex++) {
            targetLengths[referenceIndex] = 1000;
        }
        writer.setTargetLengths(targetLengths);
        // we write this alignment sorted:

        writer.setSorted(true);

        writer.setAlignmentEntry(0, 1, 12, 30, false, constantQueryLength);
        writer.appendEntry();

        writer.setAlignmentEntry(0, 1, 13, 30, false, constantQueryLength);
        writer.appendEntry();

        writer.setAlignmentEntry(0, 1, 13, 30, false, constantQueryLength);
        writer.appendEntry();

        writer.setAlignmentEntry(0, 2, 123, 30, false, constantQueryLength);
        writer.appendEntry();
        writer.setAlignmentEntry(0, 2, 300, 30, false, constantQueryLength);
        writer.appendEntry();
        writer.setAlignmentEntry(0, 2, 300, 30, false, constantQueryLength);
        writer.appendEntry();
        writer.close();
        writer.printStats(System.out);

        final AlignmentReader reader =
                new AlignmentReaderImpl(FilenameUtils.concat(BASE_TEST_DIR, basename));

        final Alignments.AlignmentEntry a = reader.skipTo(0, 0);
        assertNotNull(a);
        assertEquals(1, a.getTargetIndex());
        assertEquals(12, a.getPosition());
        final Alignments.AlignmentEntry b = reader.skipTo(0, 0);
        assertEquals(1, b.getTargetIndex());
        assertEquals(13, b.getPosition());

        final Alignments.AlignmentEntry c = reader.skipTo(0, 0);
        assertEquals(1, c.getTargetIndex());
        assertEquals(13, c.getPosition());

        final Alignments.AlignmentEntry d = reader.skipTo(2, 300);
        assertEquals(2, d.getTargetIndex());
        assertEquals(300, d.getPosition());

        final Alignments.AlignmentEntry e = reader.skipTo(2, 300);
        assertEquals(2, e.getTargetIndex());
        assertEquals(300, e.getPosition());
        assertFalse(reader.hasNext());


    }

    @Test
    public void testFewSkips1_WithConcat() throws IOException {
        final String basename = "align-skip-to-1-concat";
        final AlignmentWriterImpl writer =
                new AlignmentWriterImpl(FilenameUtils.concat(BASE_TEST_DIR, basename));
        writer.setNumAlignmentEntriesPerChunk(numEntriesPerChunk);

        final int numTargets = 3;
        final int[] targetLengths = new int[numTargets];

        for (int referenceIndex = 0; referenceIndex < numTargets; referenceIndex++) {
            targetLengths[referenceIndex] = 1000;
        }
        writer.setTargetLengths(targetLengths);
        // we write this alignment sorted:

        writer.setSorted(true);

        writer.setAlignmentEntry(0, 1, 12, 30, false, constantQueryLength);
        writer.appendEntry();

        writer.setAlignmentEntry(0, 1, 13, 30, false, constantQueryLength);
        writer.appendEntry();

        writer.setAlignmentEntry(0, 1, 13, 30, false, constantQueryLength);
        writer.appendEntry();

        writer.setAlignmentEntry(0, 2, 123, 30, false, constantQueryLength);
        writer.appendEntry();
        writer.setAlignmentEntry(0, 2, 300, 30, false, constantQueryLength);
        writer.appendEntry();
        writer.setAlignmentEntry(0, 2, 300, 30, false, constantQueryLength);
        writer.appendEntry();
        writer.close();
        writer.printStats(System.out);

        final ConcatSortedAlignmentReader reader = new ConcatSortedAlignmentReader(
                FilenameUtils.concat(BASE_TEST_DIR, basename),
                FilenameUtils.concat(BASE_TEST_DIR, basename));


        Alignments.AlignmentEntry a = reader.skipTo(0, 0);
        assertNotNull(a);
        assertEquals(1, a.getTargetIndex());
        assertEquals(12, a.getPosition());
        a = reader.skipTo(0, 0);
        assertNotNull(a);
        assertEquals(1, a.getTargetIndex());
        assertEquals(12, a.getPosition());

        final Alignments.AlignmentEntry b = reader.skipTo(0, 0);
        assertEquals(1, b.getTargetIndex());
        assertEquals(13, b.getPosition());

        Alignments.AlignmentEntry c = reader.next();
        assertEquals(1, c.getTargetIndex());
        assertEquals(13, c.getPosition());

        c = reader.next();
        assertEquals(1, c.getTargetIndex());
        assertEquals(13, c.getPosition());
        c = reader.next();
        assertEquals(1, c.getTargetIndex());
        assertEquals(13, c.getPosition());

        Alignments.AlignmentEntry d = reader.skipTo(2, 300);
        assertEquals(2, d.getTargetIndex());
        assertEquals(300, d.getPosition());

        d = reader.skipTo(2, 300);
        assertEquals(2, d.getTargetIndex());
        assertEquals(300, d.getPosition());

        Alignments.AlignmentEntry e = reader.skipTo(2, 300);
        assertNotNull(e);
        assertEquals(2, e.getTargetIndex());
        assertEquals(300, e.getPosition());
        e = reader.skipTo(2, 300);
        assertEquals(2, e.getTargetIndex());
        assertEquals(300, e.getPosition());
        assertFalse(reader.hasNext());


    }

    @Test
    public void testFewSkips2() throws IOException {
        final String basename = "align-skip-to-2";
        final AlignmentWriterImpl writer =
                new AlignmentWriterImpl(FilenameUtils.concat(BASE_TEST_DIR, basename));
        writer.setNumAlignmentEntriesPerChunk(numEntriesPerChunk);

        final int numTargets = 3;
        final int[] targetLengths = new int[numTargets];

        for (int referenceIndex = 0; referenceIndex < numTargets; referenceIndex++) {
            targetLengths[referenceIndex] = 1000;
        }
        writer.setTargetLengths(targetLengths);
        // we write this alignment sorted:

        writer.setSorted(true);

        writer.setAlignmentEntry(0, 1, 12, 30, false, constantQueryLength);
        writer.appendEntry();

        writer.setAlignmentEntry(0, 1, 13, 30, false, constantQueryLength);
        writer.appendEntry();

        writer.setAlignmentEntry(0, 1, 13, 30, false, constantQueryLength);
        writer.appendEntry();

        writer.setAlignmentEntry(0, 2, 123, 30, false, constantQueryLength);
        writer.appendEntry();
        writer.setAlignmentEntry(0, 2, 300, 30, false, constantQueryLength);
        writer.appendEntry();
        writer.setAlignmentEntry(0, 2, 300, 30, false, constantQueryLength);
        writer.appendEntry();
        writer.close();
        writer.printStats(System.out);

        final AlignmentReader reader =
                new AlignmentReaderImpl(FilenameUtils.concat(BASE_TEST_DIR, basename));


        final Alignments.AlignmentEntry c = reader.skipTo(2, 0);
        assertEquals(2, c.getTargetIndex());
        assertEquals(123, c.getPosition());


        final Alignments.AlignmentEntry d = reader.skipTo(2, 300);
        assertEquals(2, d.getTargetIndex());
        assertEquals(300, d.getPosition());

        final Alignments.AlignmentEntry e = reader.skipTo(2, 300);
        assertEquals(2, e.getTargetIndex());
        assertEquals(300, e.getPosition());
        assertFalse(reader.hasNext());


    }

    @Test
    public void testEmptyAlignment() throws IOException {
        final String basename = "align-skip-to-3";
        final AlignmentWriterImpl writer =
                new AlignmentWriterImpl(FilenameUtils.concat(BASE_TEST_DIR, basename));
        writer.setNumAlignmentEntriesPerChunk(numEntriesPerChunk);

        final int numTargets = 3;
        final int[] targetLengths = new int[numTargets];

        for (int referenceIndex = 0; referenceIndex < numTargets; referenceIndex++) {
            targetLengths[referenceIndex] = 1000;
        }
        writer.setTargetLengths(targetLengths);
        // we write this alignment sorted:

        writer.setSorted(true);

        writer.close();
        writer.printStats(System.out);

        final AlignmentReader reader =
                new AlignmentReaderImpl(FilenameUtils.concat(BASE_TEST_DIR, basename));

        final Alignments.AlignmentEntry c = reader.skipTo(2, 0);
        assertNull(c);

        assertFalse(reader.hasNext());


    }


    @Test
    public void testFewSkips4() throws IOException {
        final String basename = "align-skip-to-2";
        final AlignmentWriterImpl writer =
                new AlignmentWriterImpl(FilenameUtils.concat(BASE_TEST_DIR, basename));
        writer.setNumAlignmentEntriesPerChunk(numEntriesPerChunk);

        final int numTargets = 3;
        final int[] targetLengths = new int[numTargets];

        for (int referenceIndex = 0; referenceIndex < numTargets; referenceIndex++) {
            targetLengths[referenceIndex] = 1000;
        }
        writer.setTargetLengths(targetLengths);
        // we write this alignment sorted:

        writer.setSorted(true);

        writer.setAlignmentEntry(0, 1, 12, 30, false, constantQueryLength);
        writer.appendEntry();

        writer.setAlignmentEntry(0, 1, 13, 30, false, constantQueryLength);
        writer.appendEntry();

        writer.setAlignmentEntry(0, 1, 13, 30, false, constantQueryLength);
        writer.appendEntry();

        writer.setAlignmentEntry(0, 2, 123, 30, false, constantQueryLength);
        writer.appendEntry();
        writer.setAlignmentEntry(0, 2, 300, 30, false, constantQueryLength);
        writer.appendEntry();
        writer.setAlignmentEntry(0, 2, 300, 30, false, constantQueryLength);
        writer.appendEntry();
        writer.close();
        writer.printStats(System.out);

        final AlignmentReader reader =
                new AlignmentReaderImpl(FilenameUtils.concat(BASE_TEST_DIR, basename));


        final Alignments.AlignmentEntry c = reader.skipTo(1, 13);
        assertEquals(1, c.getTargetIndex());
        assertEquals(13, c.getPosition());


        final Alignments.AlignmentEntry d = reader.next();
        assertEquals(1, d.getTargetIndex());
        assertEquals(13, d.getPosition());

        final Alignments.AlignmentEntry e = reader.next();
        assertEquals(2, e.getTargetIndex());
        assertEquals(123, e.getPosition());


    }

    @Test
    public void testSkipWithUrl() throws IOException {

        // AlignmentReader reader = new AlignmentReaderImpl("http://dl.dropbox.com/u/357497/UANMNXR-hybrid-domain.header");
        //  AlignmentReader reader = new AlignmentReaderImpl("/data/igv-test/UANMNXR-hybrid-domain-reindexed.entries");
        AlignmentReader reader = new AlignmentReaderImpl("http://dl.dropbox.com/u/357497/EJOYQAZ-small-hybrid.entries");
        reader.readHeader();
        reader.reposition(0, 1256375);
        Alignments.AlignmentEntry entry = reader.skipTo(0, 1256375);
        assertNotNull(entry);
        assertEquals(0, entry.getTargetIndex());
        assertEquals(16676, entry.getQueryIndex());
        assertTrue(1256375 <= entry.getPosition());

    }

     @Test
    public void testOldHybrid() throws IOException {

        // AlignmentReader reader = new AlignmentReaderImpl("http://dl.dropbox.com/u/357497/UANMNXR-hybrid-domain.header");
        //  AlignmentReader reader = new AlignmentReaderImpl("/data/igv-test/UANMNXR-hybrid-domain-reindexed.entries");
        AlignmentReader reader = new AlignmentReaderImpl("http://gobyweb.apps.campagnelab.org/data/H_T_D/MYHZZJH/MYHZZJH-hybrid-domain.entries");
        reader.readHeader();
        reader.reposition(0, 1256375);
        Alignments.AlignmentEntry entry = reader.skipTo(0, 1256375);
        assertNotNull(entry);
        assertEquals(0, entry.getTargetIndex());
        assertEquals(1160266, entry.getQueryIndex());
        assertTrue(1256375 <= entry.getPosition());

    }
    /*
    chr11:67,501,982-67,505,747
     */

       @Test
    public void testOldHybridHZ() throws IOException {

        // AlignmentReader reader = new AlignmentReaderImpl("http://dl.dropbox.com/u/357497/UANMNXR-hybrid-domain.header");
        //  AlignmentReader reader = new AlignmentReaderImpl("/data/igv-test/UANMNXR-hybrid-domain-reindexed.entries");
        AlignmentReader reader = new AlignmentReaderImpl("http://gobyweb.apps.campagnelab.org/data/H_T_D/HZFWPTI/HZFWPTI-hybrid-domain.entries");
        reader.readHeader();
        reader.reposition(10, 67501982);
        Alignments.AlignmentEntry entry = reader.skipTo(10, 67501982);
        assertNotNull(entry);
        assertEquals(10, entry.getTargetIndex());
        assertEquals(1921857, entry.getQueryIndex());
        assertTrue(67501982 <= entry.getPosition());

    }
    @Test
    public void testSkipToHybrid() throws IOException {

        AlignmentReader reader = new AlignmentReaderImpl("test-data/alignment-hybrid-codec/EJOYQAZ-small-hybrid.entries");
        reader.readHeader();
        reader.reposition(0, 1014810);
        Alignments.AlignmentEntry entry = reader.skipTo(0, 1014810);
        assertNotNull(entry);
        assertEquals(0, entry.getTargetIndex());
        assertTrue(1014810 <= entry.getPosition());

    }




    @Test
    public void testHybridWindow() throws IOException {

        AlignmentReader reader = new AlignmentReaderImpl(0x304,0x1999,"test-data/alignment-hybrid-codec/EJOYQAZ-small-hybrid.entries");
        reader.readHeader();

        Alignments.AlignmentEntry entry = reader.next();
        assertNotNull(entry);
        System.out.println(entry.getQueryIndex());
        System.out.flush();
        assertEquals(97, entry.getQueryIndex());


    }

    @Test
    public void testUrlRange() throws IOException {

        long startOffset = 0;
        RepositionableInputStream is = new RepositionableInputStream("http://dl.dropbox.com/u/357497/UANMNXR-hybrid-domain-reindexed.entries");
        readFromOffset(is, startOffset);
        readFromOffset(is, 2);
        readFromOffset(is, 9);

    }

    private void readFromOffset(RepositionableInputStream is, long startOffset) throws IOException {
        System.out.println("Reading from " + startOffset);
        is.position(startOffset);
        read10(is, startOffset);
    }


    private void read10(InputStream is, long startOffset) throws IOException {
        DataInputStream dis = new DataInputStream(is);
        for (int i = 0; i < 10; i++) {
            final byte b = dis.readByte();
            System.out.println(Byte.toString(b));
            if (i == 9) {
                switch ((int) startOffset) {
                    case 0:
                        assertEquals(2, b);
                        break;
                    case 2:
                        assertEquals(-59, b);
                        break;
                    case 9:
                        assertEquals(0, b);
                        break;

                }
            }
        }
    }

}
