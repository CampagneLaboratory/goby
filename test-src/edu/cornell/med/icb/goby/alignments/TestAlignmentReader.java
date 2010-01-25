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
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.lang.RandomStringUtils;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;
import org.junit.Test;

/**
 * Basic tests for the {@link edu.cornell.med.icb.goby.alignments.AlignmentReader}.
 */
public class TestAlignmentReader {
    /**
     * Validate that the method
     * {@link edu.cornell.med.icb.goby.alignments.AlignmentReader#getBasename(String)}
     * produces the proper results.
     */
    @Test
    public void basename() {
        assertNull("Basename should be null", AlignmentReader.getBasename(null));
        assertEquals("Basename should be unchanged", "", AlignmentReader.getBasename(""));
        assertEquals("Basename should be unchanged",
                "foobar", AlignmentReader.getBasename("foobar"));
        assertEquals("Basename should be unchanged",
                "foobar.txt", AlignmentReader.getBasename("foobar.txt"));

        for (final String extension : FileExtensionHelper.COMPACT_ALIGNMENT_FILE_EXTS) {
            final String basename = RandomStringUtils.randomAlphabetic(8);
            final String filename = basename + extension;
            assertEquals("Basename not stripped properly from '" + filename + "'",
                    basename, AlignmentReader.getBasename(basename));
        }

        assertEquals("Only the extension should have been removed",
                "foo.entries.bar", AlignmentReader.getBasename("foo.entries.bar.entries"));
        assertEquals("Only the extension should have been removed",
                "entries.foo.bar", AlignmentReader.getBasename("entries.foo.bar.entries"));
    }

    /**
     * Validate that the method
     * {@link edu.cornell.med.icb.goby.alignments.AlignmentReader#getBasenames(String[])}
     * produces the proper results.
     */
    @Test
    public void basenames() {
        assertTrue("Basename array should be empty",
                ArrayUtils.isEmpty(AlignmentReader.getBasenames(null)));

        final String[] nullArray = { null };
        assertArrayEquals("Basename array should contain a single null element",
                nullArray, AlignmentReader.getBasenames((String) null));

        final String[] emptyStringArray = { "" };
        assertArrayEquals("Basename array should contain a single empty string element",
                emptyStringArray, AlignmentReader.getBasenames(""));

        final String[] foobarArray = { "foobar" };
        assertArrayEquals("Basenames should be unchanged",
                foobarArray, AlignmentReader.getBasenames("foobar"));

        final String[] foobarTxtArray = { "foobar.txt" };
        assertArrayEquals("Basenames should be unchanged",
                foobarTxtArray, AlignmentReader.getBasenames("foobar.txt"));

        assertArrayEquals("Basenames should be unchanged",
                ArrayUtils.addAll(foobarArray, foobarTxtArray),
                AlignmentReader.getBasenames("foobar", "foobar.txt"));

        final String basename = "mybasename";
        final String[] basenameArray = { basename };
        final String[] filenames =
                new String[FileExtensionHelper.COMPACT_ALIGNMENT_FILE_EXTS.length];
        for (int i = 0; i < FileExtensionHelper.COMPACT_ALIGNMENT_FILE_EXTS.length; i++) {
            filenames[i] = basename + FileExtensionHelper.COMPACT_ALIGNMENT_FILE_EXTS[i];

        }
        assertArrayEquals("Basename not stripped properly from " + ArrayUtils.toString(filenames),
                basenameArray, AlignmentReader.getBasenames(filenames));
    }
}
