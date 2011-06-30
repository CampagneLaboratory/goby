/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
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
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: May 11, 2011
 *         Time: 10:41:31 AM
 */
public class TestSpliceSites {
    @Test
    /**
     * Make a synthetic alignment with one read that spans a splice site:
     */
    public void writeSplicedTestAlignment() throws IOException {
        final String basename = "spliced-1";
        final AlignmentWriter writer =
                new AlignmentWriter(FilenameUtils.concat(BASE_TEST_DIR, basename));
        writer.setNumAlignmentEntriesPerChunk(numEntriesPerChunk);

        final int numTargets = 1;
        final int[] targetLengths = new int[numTargets];

        for (int referenceIndex = 0; referenceIndex < numTargets; referenceIndex++) {
            targetLengths[referenceIndex] = 110;
        }
        writer.setTargetLengths(targetLengths);
        // we write this alignment sorted:
        IndexedIdentifier ids = new IndexedIdentifier();
        ids.registerIdentifier(new MutableString("synth1"));
        writer.setTargetIdentifiers(ids);
        writer.setSorted(true);

        Alignments.SequenceVariation mutation = Alignments.SequenceVariation.newBuilder().setFrom("T").setTo("A").
                setPosition(10).setReadIndex(10).build();
        Alignments.RelatedAlignmentEntry linkForward = Alignments.RelatedAlignmentEntry.newBuilder().
                setFragmentIndex(1).setPosition(70).setTargetIndex(0).build();

        Alignments.RelatedAlignmentEntry linkBackward = Alignments.RelatedAlignmentEntry.newBuilder().
                setFragmentIndex(0).setPosition(10).setTargetIndex(0).build();
        Alignments.AlignmentEntry entry1 = Alignments.AlignmentEntry.newBuilder().setPosition(10).setMatchingReverseStrand(false).
                setQueryLength(50).setQueryIndex(0).setTargetIndex(0).setFragmentIndex(0).setSplicedForwardAlignmentLink(linkForward).
                setQueryAlignedLength(20).addSequenceVariations(mutation).setMappingQuality(255).setScore(20).build();
        Alignments.AlignmentEntry entry2 = Alignments.AlignmentEntry.newBuilder().setPosition(70).setMatchingReverseStrand(false).
                setQueryLength(50).setQueryIndex(0).setTargetIndex(0).setFragmentIndex(1).setSplicedBackwardAlignmentLink(linkBackward).
                setQueryAlignedLength(30).setMappingQuality(255).setScore(30).build();

        writer.appendEntry(entry1);
        writer.appendEntry(entry2);

        // now write on the reverse strand:
        mutation = Alignments.SequenceVariation.newBuilder().setFrom("T").setTo("A").
                setPosition(10).setReadIndex(10).build();
        linkForward = Alignments.RelatedAlignmentEntry.newBuilder().
                setFragmentIndex(0).setPosition(70).setTargetIndex(0).build();

        linkBackward = Alignments.RelatedAlignmentEntry.newBuilder().
                setFragmentIndex(1).setPosition(10).setTargetIndex(0).build();
        entry1 = Alignments.AlignmentEntry.newBuilder().setPosition(10).setMatchingReverseStrand(true).
                setQueryLength(50).setQueryIndex(0).setTargetIndex(0).setFragmentIndex(1).setSplicedForwardAlignmentLink(linkForward).
                setQueryAlignedLength(20).addSequenceVariations(mutation).setMappingQuality(255).setScore(20).build();
        entry2 = Alignments.AlignmentEntry.newBuilder().setPosition(70).setMatchingReverseStrand(true).
                setQueryLength(50).setQueryIndex(0).setTargetIndex(0).setFragmentIndex(0).setSplicedBackwardAlignmentLink(linkBackward).
                setQueryAlignedLength(30).setMappingQuality(255).setScore(30).build();
        writer.appendEntry(entry1);
        writer.appendEntry(entry2);

        writer.close();

    }

    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(TestSkipTo.class);
    private static final String BASE_TEST_DIR = "test-results/alignments-splices";
    private int numEntriesPerChunk = 1000;


    @BeforeClass
    public static void initializeTestDirectory() throws IOException {
        if (LOG.isDebugEnabled()) {
            LOG.debug("Creating base test directory: " + BASE_TEST_DIR);
        }
        FileUtils.forceMkdir(new File(BASE_TEST_DIR));
    }


    /*
    This test is used to create the spliced dataset used for IGV testing. We want to keep the file. Do not delete.
      @AfterClass
    public static void cleanupTestDirectory() throws IOException {
        if (LOG.isDebugEnabled()) {
            LOG.debug("Deleting base test directory: " + BASE_TEST_DIR);
        }
        //    FileUtils.forceDeleteOnExit(new File(BASE_TEST_DIR));
    }  */

}
