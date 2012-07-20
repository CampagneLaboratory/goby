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

import edu.cornell.med.icb.goby.util.FileExtensionHelper;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;
import junit.framework.Assert;
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.lang.RandomStringUtils;

import static junitx.framework.Assert.assertNotEquals;
import static org.junit.Assert.*;

import org.junit.Test;

import java.io.IOException;

/**
 * Basic tests for the {@link AlignmentReaderImpl}.
 */
public class TestAlignmentReader {
    String basename = "test-data/bam/Example.entries";

    @Test
    public void testWithSlice() throws IOException {
        long end;
        long start;

        final AlignmentReaderImpl readerFromStart = new AlignmentReaderImpl(0L, 1L, basename);
        assertTrue(readerFromStart.hasNext());
        final Alignments.AlignmentEntry firstEntry = readerFromStart.next();
        final AlignmentReaderImpl readerAt10K = new AlignmentReaderImpl(146387L, 146387L + 1, basename);
        assertTrue(readerAt10K.hasNext());
        final Alignments.AlignmentEntry firstEntry10KSlice = readerAt10K.next();

        assertNotEquals(firstEntry10KSlice, firstEntry);
    }

    @Test
    public void anotherTest() throws IOException {
        final AlignmentReaderImpl readerFromStart = new AlignmentReaderImpl(0L, 1L, basename);
        assertTrue(readerFromStart.hasNext());
        final Alignments.AlignmentEntry firstEntry = readerFromStart.next();
        readerFromStart.close();

        final IntSet entrySet0 = new IntOpenHashSet();
        String basename = "test-data/bam/Example.entries";
        IterateAlignments it = new IterateAlignments() {
            @Override
            public void processAlignmentEntry(AlignmentReader alignmentReader, Alignments.AlignmentEntry alignmentEntry) {
                if (entrySet0.isEmpty()) {
                    System.out.println(alignmentEntry);
                    entrySet0.add(alignmentEntry.getQueryIndex());
                }
            }
        };
        //    it.iterate(0L, 1L, basename);
        //   System.out.println("----------------");
        it.iterate(146387L, 146387L + 1, basename);
        assertTrue(!entrySet0.isEmpty());

        assertNotEquals("query indices of first entry must differ", entrySet0.iterator().nextInt(), firstEntry.getQueryIndex());

    }

    /**
     * Validate that the method
     * {@link AlignmentReaderImpl#getBasename(String)}
     * produces the proper results.
     */
    @Test
    public void basename() {
        assertNull("Basename should be null", AlignmentReaderImpl.getBasename(null));
        assertEquals("Basename should be unchanged", "", AlignmentReaderImpl.getBasename(""));
        assertEquals("Basename should be unchanged",
                "foobar", AlignmentReaderImpl.getBasename("foobar"));
        assertEquals("Basename should be unchanged",
                "foobar.txt", AlignmentReaderImpl.getBasename("foobar.txt"));

        for (final String extension : FileExtensionHelper.COMPACT_ALIGNMENT_FILE_EXTS) {
            final String basename = RandomStringUtils.randomAlphabetic(8);
            final String filename = basename + extension;
            assertEquals("Basename not stripped properly from '" + filename + "'",
                    basename, AlignmentReaderImpl.getBasename(basename));
        }

        assertEquals("Only the extension should have been removed",
                "foo.entries.bar", AlignmentReaderImpl.getBasename("foo.entries.bar.entries"));
        assertEquals("Only the extension should have been removed",
                "entries.foo.bar", AlignmentReaderImpl.getBasename("entries.foo.bar.entries"));
    }

    /**
     * Validate that the method
     * {@link AlignmentReaderImpl#getBasenames(String[])}
     * produces the proper results.
     */
    @Test
    public void basenames() {
        assertTrue("Basename array should be empty",
                ArrayUtils.isEmpty(AlignmentReaderImpl.getBasenames()));

        final String[] nullArray = {null};
        assertArrayEquals("Basename array should contain a single null element",
                nullArray, AlignmentReaderImpl.getBasenames((String) null));

        final String[] emptyStringArray = {""};
        assertArrayEquals("Basename array should contain a single empty string element",
                emptyStringArray, AlignmentReaderImpl.getBasenames(""));

        final String[] foobarArray = {"foobar"};
        assertArrayEquals("Basenames should be unchanged",
                foobarArray, AlignmentReaderImpl.getBasenames("foobar"));

        final String[] foobarTxtArray = {"foobar.txt"};
        assertArrayEquals("Basenames should be unchanged",
                foobarTxtArray, AlignmentReaderImpl.getBasenames("foobar.txt"));

        assertArrayEquals("Basenames should be unchanged",
                ArrayUtils.addAll(foobarArray, foobarTxtArray),
                AlignmentReaderImpl.getBasenames("foobar", "foobar.txt"));

        final String basename = "mybasename";
        final String[] basenameArray = {basename};
        final String[] filenames =
                new String[FileExtensionHelper.COMPACT_ALIGNMENT_FILE_EXTS.length];
        for (int i = 0; i < FileExtensionHelper.COMPACT_ALIGNMENT_FILE_EXTS.length; i++) {
            filenames[i] = basename + FileExtensionHelper.COMPACT_ALIGNMENT_FILE_EXTS[i];

        }
        assertArrayEquals("Basename not stripped properly from " + ArrayUtils.toString(filenames),
                basenameArray, AlignmentReaderImpl.getBasenames(filenames));
    }


    @Test
    public void canRead() {

        assertFalse(AlignmentReaderImpl.canRead("https://dm.genomespace.org/datamanager/file/Home/igvtest/breasttumor.acgh.info.txt"));

    }
}
