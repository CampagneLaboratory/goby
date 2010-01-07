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

package edu.cornell.med.icb.goby.aligners;

import edu.cornell.med.icb.goby.config.GobyConfiguration;
import edu.cornell.med.icb.goby.readers.FastXEntry;
import edu.cornell.med.icb.goby.readers.FastXReader;
import edu.cornell.med.icb.goby.reads.ReadsWriter;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.junit.AfterClass;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

/**
 * @author Fabien Campagne
 *         Date: Jul 9, 2009
 *         Time: 3:54:42 PM
 */
public class TestAligner {
    private static final Log LOG = LogFactory.getLog(TestAligner.class);
    private static final String BASE_TEST_DIR = "test-results/aligners";
    private static final String[] reads = {
            "ACTGCGCGCG",
            "AAAAATTTTGGGGGCCCCCCC",
            "AAAAATTTTGGGGGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
    };
    final String[] references = {
            "AGTGCGCGCG",
            "AAAAATTTTGGGGGCCCCCCC",
            "AAAAATTTTGGGGGCCTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
    };

    @Test
    public void testAlignWithLastag() throws IOException, InterruptedException {
        final LastagAligner aligner = new LastagAligner();
        final String databaseDirectory = FilenameUtils.concat(BASE_TEST_DIR, "db-lastag");
        FileUtils.forceMkdir(new File(databaseDirectory));
        aligner.setConfiguration(GobyConfiguration.getConfiguration());

        aligner.setDatabaseDirectory(databaseDirectory);
        aligner.setWorkDirectory(databaseDirectory);

        final File compactReads =
                new File(FilenameUtils.concat(BASE_TEST_DIR, "lastag-input-reads.compact-reads"));
        final File compactReference =
                new File(FilenameUtils.concat(BASE_TEST_DIR, "lastag-input-reference.compact-reads"));

        writeCompact(compactReads, reads);
        writeCompact(compactReference, reads); // SAME AS READS!

        final File fastaFile = aligner.prepareReads(compactReads);
        assertNotNull("Reads fasta file should not be null", fastaFile);
        assertTrue("Reads fasta file should exist", fastaFile.exists());

        int readCount = 0;
        FastXReader fastaReader = null;
        try {
            fastaReader = new FastXReader(fastaFile.getAbsolutePath());
            assertEquals("File should be in fasta format", "fa", fastaReader.getFileType());
            for (final FastXEntry entry : fastaReader) {
                assertEquals("Mismatch for entry " + readCount,
                        reads[readCount], entry.getSequence().toString());
                readCount++;
            }
        } finally {
            if (fastaReader != null) {
                fastaReader.close();
            }
        }
        assertEquals("Number of sequences in fasta file do not match", reads.length, readCount);

        aligner.setDatabaseName(aligner.getDefaultDbNameForReferenceFile(compactReference));
        aligner.align(compactReference, compactReads, FilenameUtils.concat(BASE_TEST_DIR, "lastag-output"));
    }

    @Test
    public void testAlignWithBWA() throws IOException, InterruptedException {
        final BWAAligner aligner = new BWAAligner();
        final String databaseDirectory = FilenameUtils.concat(BASE_TEST_DIR, "db-bwa");
        FileUtils.forceMkdir(new File(databaseDirectory));
        aligner.setConfiguration(GobyConfiguration.getConfiguration());

        aligner.setDatabaseDirectory(databaseDirectory);
        aligner.setWorkDirectory(databaseDirectory);

        final File compactReads =
                new File(FilenameUtils.concat(BASE_TEST_DIR, "bwa-input-reads.compact-reads"));
        final File compactReference =
                new File(FilenameUtils.concat(BASE_TEST_DIR, "bwa-input-reference.compact-reads"));

        writeCompact(compactReads, reads);
        writeCompact(compactReference, references);  // Different from reads

        aligner.setDatabaseName(aligner.getDefaultDbNameForReferenceFile(compactReference));
        aligner.align(compactReference, compactReads, FilenameUtils.concat(BASE_TEST_DIR, "bwa-output"));
    }

    private void writeCompact(final File compactFile, final String[] sequences) throws IOException {
        int index = 1;

        ReadsWriter writer = null;
        OutputStream compactFileStream = null;
        try {
            compactFileStream = FileUtils.openOutputStream(compactFile);
            writer = new ReadsWriter(compactFileStream);
            writer.setNumEntriesPerChunk(9);

            for (final String sequence : sequences) {
                writer.setSequence(sequence);
                writer.setIdentifier(Integer.toString(index++));
                writer.appendEntry();
            }
        } finally {
            if (writer != null) {
                writer.close();
            }
            IOUtils.closeQuietly(compactFileStream);
        }
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
}
