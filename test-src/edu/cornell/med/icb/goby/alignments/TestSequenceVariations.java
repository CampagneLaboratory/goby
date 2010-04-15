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

import edu.cornell.med.icb.goby.modes.AlignMode;
import edu.cornell.med.icb.goby.reads.ReadsWriter;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.junit.AfterClass;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: Apr 30, 2009
 *         Time: 6:15:52 PM
 */
public class TestSequenceVariations {
    private static final Log LOG = LogFactory.getLog(TestSequenceVariations.class);
    private static final String BASE_TEST_DIR = "test-results/sequence-variations";
    private final String referenceFilename =
            FilenameUtils.concat(BASE_TEST_DIR, "reference-sequence-vars.compact-reads");
    private final String readsFilename =
            FilenameUtils.concat(BASE_TEST_DIR, "query-sequence-vars.compact-reads");
    private String lastagAlignmentFilename;
    private String bwaAlignmentFilename;

    @BeforeClass
    public static void initializeTestDirectory() throws IOException {
        if (LOG.isDebugEnabled()) {
            LOG.debug("Creating base test directory: " + BASE_TEST_DIR);
        }

        FileUtils.forceMkdir(new File(BASE_TEST_DIR));
    }

    @AfterClass
    public static void cleanupTestDirectory()  {
        if (LOG.isDebugEnabled()) {
            LOG.debug("Deleting base test directory: " + BASE_TEST_DIR);
        }
        // FileUtils.forceDeleteOnExit(new File(BASE_TEST_DIR));
    }

    @Test
    public void testLastagSequenceVariationParsing() throws IOException {
        final AlignmentReader reader = new AlignmentReader(lastagAlignmentFilename);
        while (reader.hasNext()) {
            final Alignments.AlignmentEntry alignmentEntry = reader.next();

            assertTrue("alignment must have variation", alignmentEntry.getSequenceVariationsCount() > 0);
            int variationIndex = 0;
            for (final Alignments.SequenceVariation var : alignmentEntry.getSequenceVariationsList()) {
                variationIndex++;
                System.out.println(String.format("last entry score=%f referenceIndex=%d  queryIndex=%d variation: %s",
                        alignmentEntry.getScore(),
                        alignmentEntry.getQueryIndex(),
                        alignmentEntry.getTargetIndex(),
                        var.toString()));
                assertLength(alignmentEntry, var);

                switch (alignmentEntry.getQueryIndex()) {
                    case 7:
                        //last finds this alignment. We test the case when the read does not match the reference at the beginning:
                        assertEquals(19, var.getPosition());
                        assertEquals(21, var.getReadIndex());
                        assertEquals("A", var.getFrom());
                        assertEquals("T", var.getTo());

                        break;
                    case 5:
                        //last finds this alignment. We test the case when the read does not match the reference at the beginning:
                        assertEquals(24, var.getPosition());
                        assertEquals(24, var.getReadIndex());
                        assertEquals("T", var.getFrom());
                        assertEquals("-", var.getTo());

                        break;

                    case 1:
                    case 4:
                        //last finds this alignment. We test deletion of TCC from the reference.
                        assertEquals(14, var.getPosition());
                        assertEquals(14, var.getReadIndex());
                        assertEquals("TCC", var.getFrom());
                        assertEquals("---", var.getTo());
                        break;
                    case 3:
                        //last finds this alignment. We test deletion of TCC from the reference.
                        assertEquals(8, var.getPosition());
                        assertEquals(9, var.getReadIndex());
                        assertEquals("-", var.getFrom());
                        assertEquals("C", var.getTo());

                        break;
                    case 0:
                        //last finds this alignment. We test deletion of TCC from the reference.
                        assertEquals(8, var.getPosition());
                        assertEquals(8, var.getReadIndex());
                        assertEquals("--", var.getFrom());
                        assertEquals("CC", var.getTo());
                        break;
                    case 8:
                        switch (variationIndex) {
                            case 1:
                                assertEquals(36, var.getPosition());
                                assertEquals(36, var.getReadIndex());
                                assertEquals("A", var.getFrom());
                                assertEquals("G", var.getTo());
                                break;

                        }
                        break;
                }

            }
        }
    }

    private void assertLength(final Alignments.AlignmentEntry alignmentEntry, final Alignments.SequenceVariation var) {
        final String readRaw = alignments[alignmentEntry.getQueryIndex()].read;
        final String readProcessed = readRaw.replaceAll("-", "");
        assertTrue(String.format("read index %d must be less than read length %d.", var.getReadIndex(), readProcessed.length()),
                var.getReadIndex() < readProcessed.length());
    }


    @Test
    public void testBwaSequenceVariationParsing() throws IOException {
        final AlignmentReader reader = new AlignmentReader(bwaAlignmentFilename);
        while (reader.hasNext()) {
            final Alignments.AlignmentEntry alignmentEntry = reader.next();

            assertTrue("alignment must have variation", alignmentEntry.getSequenceVariationsCount() > 0);

            for (final Alignments.SequenceVariation var : alignmentEntry.getSequenceVariationsList()) {
                System.out.println(String.format("bwa entry score=%f referenceIndex=%d  queryIndex=%d variation: %s",
                        alignmentEntry.getScore(),
                        alignmentEntry.getQueryIndex(),
                        alignmentEntry.getTargetIndex(),
                        var.toString()));
                assertLength(alignmentEntry, var);
                switch (alignmentEntry.getQueryIndex()) {

                    case 5:
                        //last finds this alignment. We test the case when the read does not match the reference at the beginning:
                        assertEquals(24, var.getPosition());
                        assertEquals(24, var.getReadIndex());
                        assertEquals("T", var.getFrom());
                        assertEquals("-", var.getTo());

                        break;

                    case 1:
                    case 4:
                        //last finds this alignment. We test deletion of TCC from the reference.
                        assertEquals(14, var.getPosition());
                        assertEquals(14, var.getReadIndex());
                        assertEquals("TCC", var.getFrom());
                        assertEquals("---", var.getTo());
                        break;
                    case 3:
                        //last finds this alignment. We test deletion of TCC from the reference.
                        assertEquals(8, var.getPosition());
                        assertEquals(9, var.getReadIndex());
                        assertEquals("-", var.getFrom());
                        assertEquals("C", var.getTo());

                        break;
                    case 0:
                        //last finds this alignment. We test deletion of TCC from the reference.
                        assertEquals(8, var.getPosition());
                        assertEquals(8, var.getReadIndex());
                        assertEquals("--", var.getFrom());
                        assertEquals("CC", var.getTo());

                        break;

                }
            }
        }
    }

    @Before
    public void setUp() throws IOException {
        final ReadsWriter referenceWriter = new ReadsWriter(new FileOutputStream(referenceFilename));
        final ReadsWriter queryWriter = new ReadsWriter(new FileOutputStream(readsFilename));
        for (final Alignment entry : alignments) {
            referenceWriter.setDescription("reference:" + entry.description);
            referenceWriter.setIdentifier("reference:" + entry.description);
            referenceWriter.setSequence(entry.reference.replaceAll("-", ""));
            referenceWriter.appendEntry();
            queryWriter.setDescription("read:" + entry.description);
            queryWriter.setSequence(entry.read.replaceAll("-", ""));
            queryWriter.appendEntry();

        }
        referenceWriter.close();
        queryWriter.close();
        // align with Lastag:

        lastagAlignmentFilename = alignWith("lastag");
        bwaAlignmentFilename = alignWith("bwa");
    }

    private String alignWith(final String alignerName) throws IOException {
        final String alignmentFilename = FilenameUtils.concat(BASE_TEST_DIR, alignerName + "-alignment");
        final AlignMode aligner = new AlignMode();
        aligner.setWorkDirectory(new File(BASE_TEST_DIR));
        aligner.setAlignerName(alignerName);
        aligner.setColorSpace(false);
        aligner.setReadsFile(new File(readsFilename));
        aligner.setReferenceFile(new File(referenceFilename));
        aligner.setSearchFlag(true);

        aligner.setOutputBasename(alignmentFilename);
        aligner.setAlignerOptions("matchQuality=BEST_MATCH");
        aligner.setQualityFilterParameters("threshold=1");
        aligner.setKeepTemporaryFiles(true);
        aligner.setDatabaseDirectory(new File(BASE_TEST_DIR));
        aligner.execute();
        return alignmentFilename;
    }

    private class Alignment {
        private String reference = "";
        private String read = "";
        private String description = "";

        public Alignment(final String description, final String reference, final String read) {
            this.description = description;
            this.reference = reference;
            this.read = read;
        }
    }

    private final Alignment[] alignments = {
            new Alignment("0_insertion",
                    //1234567891111111111222
                    //         0123456789012
                    "TTTAAAA--TAAAAAAAAAAAAAAACCCC",
                    "TTTAAAACCTAAAAAAAAAAAAAAACCCC"),
            //  "TTTTAAAACTAAAAAAAAAAAAAAACCCC")
            new Alignment("1_deletion",
                    //1234567891111111111222
                    //         0123456789012
                    "CCAAAAAAAAAAATCCAAAAAAAAAACCCAAAAAAAAAA",
                    "CCAAAAAAAAAAA---AAAAAAAAAACCCAAAAAAAAAA"),
            new Alignment("2_mutations",
                    //1234567891111111111222
                    //         0123456789012
                    "TTTCCCAAATTTCACATCACTACTACTACGGATACAGAACGGGG",
                    "TTTCCCACATTTCCCATCACCACTACTACGGATACAGAACGGGG"),
            //.......M.....M......M.......................
            new Alignment("3_insertion",
                    //1234567891111111111222
                    //         0123456789012
                    "NNNNNNNNNNNNNNNNNNNNNNNNNNNNN", // use Ns when the query already matches another reference.
                    "TTTTAAAACTAAAAAAAAAAAAAAACCCC"),  // has a C inserted (compare to 0_insertion), at position 8 in the reference, and 9 in the read.
            new Alignment("4_deletion",
                    //1234567891111111111222
                    //         0123456789012
                    "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN", // use Ns when the query already matches another reference.
                    "CCAAAAAAAAAAA---AAAAAAAAAACCCAAAAAAAAAA"),
            new Alignment("5_deletion",
                    //1234567891111111111222222
                    //         0123456789012345
                    "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN", // use Ns when the query already matches another reference.
                    "TTTCCCAAATTTCACATCACTAC-ACTACGGATACAGAACGGGG"),     // The reference has a T instead of the gap character at position 24.
            new Alignment("6_mutations-reversed",
                    //1234567891111111111222
                    //         0123456789012
                    "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
                    "CCCCGTTCTGTATCCGTAGTAGTGGTGATGGGAAATGTGGGAAA"),    // reverse complement of sequence 2_mutations, used to check positions for reverse strand matches
            new Alignment("7_mutations_padding",
                    //1234567891111111111222
                    //         0123456789012
                    "--CCTTCCTTCCTTCCTTCCACTATCATTTTAACTACTCATACTATCCCATATA",
                    "AACCTTCCTTCCTTCCTTCCTCTATCATTTTAACTACTCATACTATCCCATATA"),
             new Alignment("8_read_padding",
                    //1234567891111111111222
                    //         0123456789012
                    "NNNN",
                   "TTCCACTATCATTTTAACTACTCATACTATCCCATGTA"),    //   A->G  readIndex=36, position=
    };
}
