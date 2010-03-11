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

import edu.cornell.med.icb.identifier.IndexedIdentifier;
import edu.cornell.med.icb.goby.reads.ReadsWriter;
import edu.cornell.med.icb.goby.modes.AlignMode;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import org.junit.Before;
import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.io.FileOutputStream;
import java.util.regex.Pattern;
import java.util.regex.Matcher;

/**
 * @author Fabien Campagne
 *         Date: Apr 30, 2009
 *         Time: 6:15:52 PM
 */
public class TestSequenceVariations {
    private static final Log LOG = LogFactory.getLog(TestSequenceVariations.class);
    private static final String BASE_TEST_DIR = "test-results/sequence-variations";
    final String referenceFilename = FilenameUtils.concat(BASE_TEST_DIR, "reference-sequence-vars.compact-reads");
    final String readsFilename = FilenameUtils.concat(BASE_TEST_DIR, "query-sequence-vars.compact-reads");
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
    public static void cleanupTestDirectory() throws IOException {
        if (LOG.isDebugEnabled()) {
            LOG.debug("Deleting base test directory: " + BASE_TEST_DIR);
        }
        // FileUtils.forceDeleteOnExit(new File(BASE_TEST_DIR));
    }

    @Test
    public void testLastagSequenceVariationParsing() throws IOException {
        AlignmentReader reader = new AlignmentReader(lastagAlignmentFilename);
        while (reader.hasNext()) {
            Alignments.AlignmentEntry alignmentEntry = reader.next();

            assertTrue("alignment must have variation", alignmentEntry.getSequenceVariationsCount() > 0);
            for (Alignments.SequenceVariation var : alignmentEntry.getSequenceVariationsList()) {
                System.out.println(String.format("last entry score=%f referenceIndex=%d  queryIndex=%d variation: %s",
                        alignmentEntry.getScore(),
                        alignmentEntry.getQueryIndex(),
                        alignmentEntry.getTargetIndex(),
                        var.toString()));
            }
        }
    }


    @Test
    public void testBwaSequenceVariationParsing() throws IOException {
        AlignmentReader reader = new AlignmentReader(bwaAlignmentFilename);
        while (reader.hasNext()) {
            Alignments.AlignmentEntry alignmentEntry = reader.next();
            // TODO enable when BWA variation parsing has been implemented.
            assertTrue("alignment must have variation", alignmentEntry.getSequenceVariationsCount() > 0);
            for (Alignments.SequenceVariation var : alignmentEntry.getSequenceVariationsList()) {
                System.out.println(String.format("bwa entry score=%f referenceIndex=%d  queryIndex=%d variation: %s",
                        alignmentEntry.getScore(),
                        alignmentEntry.getQueryIndex(),
                        alignmentEntry.getTargetIndex(),
                        var.toString()));
            }
        }
    }

    @Before
    public void setUp() throws IOException {
        ReadsWriter referenceWriter = new ReadsWriter(new FileOutputStream(referenceFilename));
        ReadsWriter queryWriter = new ReadsWriter(new FileOutputStream(readsFilename));
        for (alignment entry : alignments) {
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

    private String alignWith(String alignerName) throws IOException {
        String alignmentFilename = FilenameUtils.concat(BASE_TEST_DIR, alignerName + "-alignment");
        AlignMode aligner = new AlignMode();
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

    private class alignment {
        String reference = "";
        String read = "";
        String description = "";

        public alignment(String description, String reference, String read) {
            this.description = description;
            this.reference = reference;
            this.read = read;
        }
    }

    alignment[] alignments = new alignment[]{
            new alignment("0_insertion",
                    //1234567891111111111222
                    //         0123456789012
                    "TTTAAAA--TAAAAAAAAAAAAAAACCCC",
                    "TTTAAAACCTAAAAAAAAAAAAAAACCCC"),
            new alignment("1_deletion",
                    //1234567891111111111222
                    //         0123456789012
                    "CCAAAAAAAAAAATCCAAAAAAAAAACCCAAAAAAAAAA",
                    "CCAAAAAAAAAAA---AAAAAAAAAACCCAAAAAAAAAA"),
            new alignment("2_mutations",
                    //1234567891111111111222
                    //         0123456789012
                    "TTTCCCAAATTTCACATCACTACTACTACGGATACAGAACGGGG",
                    "TTTCCCACATTTCCCATCACCACTACTACGGATACAGAACGGGG"),
            //.......M.....M......M.......................
            new alignment("3_insertion",
                    //1234567891111111111222
                    //         0123456789012
                    "NNNNNNNNNNNNNNNNNNNNNNNNNNNNN", // use Ns when the query already matches another reference.
                    "TTTTAAAACTAAAAAAAAAAAAAAACCCC"),
            new alignment("4_deletion",
                    //1234567891111111111222
                    //         0123456789012
                    "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN", // use Ns when the query already matches another reference.
                    "CCAAAAAAAAAAA---AAAAAAAAAACCCAAAAAAAAAA"),
            new alignment("5_deletion",
                    //1234567891111111111222
                    //         0123456789012
                    "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN", // use Ns when the query already matches another reference.
                    "TTTCCCAAATTTCACATCACTAC-ACTACGGATACAGAACGGGG"),     // The reference has a T instead of the gap character.

    };


}
