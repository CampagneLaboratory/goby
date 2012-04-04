/*
 * Copyright (C) 2009-2012 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.modes;

import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceTestSupport;
import it.unimi.dsi.lang.MutableString;
import junit.framework.Assert;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;

import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

/**
 * @author Fabien Campagne
 *         Date: 3/29/12
 *         Time: 12:11 PM
 */
public class TestSplicedSamHelper {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(TestSplicedSamHelper.class);
    private static final String BASE_TEST_DIR = "test-results/splicedsamhelper";

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


    /*
9068319	0	1	18339	2	28M6371N7M	*	0	0	CCTGCACCTGGCTCCGGCTCTGCTCTACCTGCTGA	aaa^`a``a^^aaaa``aaa_V__`_X]`a`aa_[	MD:Z:35	NH:i:3	NM:i:0	SM:i:2	XQ:i:40	X2:i:40	XS:A:-
11090122	0	1	18345	1	22M6371N13M	*	0	0	CCTGGCTCCGGCTCTGCTCTACCTGCTGAAGATGT	Xa^`U\``]]Y`ZZ\ZYZ\\\Z`ZQ\XJO\VGOQQ	MD:Z:35	NH:i:4	NM:i:0	SM:i:1	XQ:i:40	X2:i:40	XS:A:-
12986491	0	1	18345	2	22M6371N13M	*	0	0	CCTGGCTCCGGCTCTGCTCTACCTGCTGCAGATGT	__a`a_^__^`^`^^^^^^]V[``^\YTFV\XXYS	MD:Z:28A6	NH:i:3	NM:i:1	SM:i:2	XQ:i:40	X2:i:40	XS:A:-

     */

    @Test
    public void testLimits() {
        final SplicedSamHelper samHelper = new SplicedSamHelper();

        SplicedSamHelper.Limits[] limits = samHelper.getLimits(18339, "28M6371N7M", "28A6");
        assertEquals(2, limits.length);
        assertEquals(18339, limits[0].position);
        assertEquals(18339 + 28 + 6371, limits[1].position);

        assertEquals(0, limits[0].cigarStart);
        assertEquals(3, limits[0].cigarEnd);
        assertEquals(8, limits[1].cigarStart);
        assertEquals(10, limits[1].cigarEnd);

        assertEquals(0, limits[0].readStart);
        assertEquals(28, limits[0].readEnd);
        assertEquals(28, limits[1].readStart);
        assertEquals(35, limits[1].readEnd);

        assertEquals("28", limits[0].md);
        assertEquals("A6", limits[1].md);


    }

    @Test
    public void testSpliced() throws IOException {
        final SplicedSamHelper samHelper = new SplicedSamHelper();
        String bases_0_28 = "CCTGCACCTGGCTCCGGCTCTGCTCTAC";
        String quals_0_28 = "aaa^`a``a^^aaaa``aaa_V__`_X]";
        final String bases_28_35 = "CTGCTGA";
        final String sourceRead = bases_0_28 + bases_28_35;
        final String quals_28_35 = "`a`aa_[";
        final String sourceQual = quals_0_28 + quals_28_35;
        samHelper.setSource(0, sourceRead, sourceQual, "28M6371N7M", "35", 18339, false);
        assertEquals(2, samHelper.getNumEntries());
        samHelper.setEntryCursor(0);

        assertEquals(0, samHelper.getNumLeftClipped());
        assertEquals(0, samHelper.getNumRightClipped());
        assertEquals(bases_0_28, samHelper.getQuery().toString());
        assertEquals(bases_0_28, samHelper.getRef().toString());
        assertEquals(quals_0_28, samHelper.getQual().toString());
        assertEquals(18339 - 1, samHelper.getPosition());
        assertEquals(0, samHelper.getQueryPosition());
        assertEquals(28, samHelper.getScore());
        assertEquals(28, samHelper.getAlignedLength());
        assertEquals(28, samHelper.getQueryAlignedLength());
        assertEquals(28, samHelper.getTargetAlignedLength());
        assertEquals(0, samHelper.getNumMisMatches());
        assertEquals(0, samHelper.getNumInsertions());
        assertEquals(0, samHelper.getNumDeletions());
        assertEquals(0, samHelper.getQueryIndex());
        assertEquals(bases_0_28.length(), samHelper.getQueryLength());
        assertFalse(samHelper.isReverseStrand());
        List<SamSequenceVariation> vars = samHelper.getSequenceVariations();
        assertEquals(0, vars.size());

        samHelper.setEntryCursor(1);
        assertEquals(0, samHelper.getNumLeftClipped());
        assertEquals(0, samHelper.getNumRightClipped());
        assertEquals(bases_28_35, samHelper.getQuery().toString());
        assertEquals(bases_28_35, samHelper.getRef().toString());
        assertEquals(quals_28_35, samHelper.getQual().toString());
        assertEquals(18339 - 1 + 6371 + bases_0_28.length(), samHelper.getPosition());
    }

    @Test
    public void testSamToCompactTrickCase1() throws IOException {

        SAMToCompactMode importer = new SAMToCompactMode();
        importer.setInputFile("test-data/splicedsamhelper/tricky-spliced.sam");
        final String outputFilename = FilenameUtils.concat(BASE_TEST_DIR, "spliced-output-alignment-1");
        importer.setOutputFile(outputFilename);
        importer.execute();

        AlignmentReader reader = new AlignmentReaderImpl(outputFilename);
        assertTrue(reader.hasNext());
        Alignments.AlignmentEntry first = reader.next();

        assertEquals(0, first.getQueryIndex());
        assertEquals(0, first.getFragmentIndex());
        assertTrue(first.hasPairAlignmentLink());
        assertTrue(first.hasSplicedForwardAlignmentLink());
        Assert.assertEquals(1, first.getSplicedForwardAlignmentLink().getFragmentIndex());
        Assert.assertEquals(2, first.getPairAlignmentLink().getFragmentIndex());
        assertFalse(first.hasSplicedBackwardAlignmentLink());

        Alignments.AlignmentEntry second = reader.next();
        assertEquals(0, second.getQueryIndex());
        assertEquals(1, second.getFragmentIndex());
        assertTrue(second.hasPairAlignmentLink());

        Assert.assertEquals(2, second.getPairAlignmentLink().getFragmentIndex());
        assertFalse(second.hasSplicedForwardAlignmentLink());
        assertTrue(second.hasSplicedBackwardAlignmentLink());
        Assert.assertEquals(0, second.getSplicedBackwardAlignmentLink().getFragmentIndex());

        Alignments.AlignmentEntry third = reader.next();
        assertEquals(0, third.getQueryIndex());
        assertEquals(2, third.getFragmentIndex());
        assertTrue(third.hasPairAlignmentLink());

        Assert.assertEquals(1, third.getPairAlignmentLink().getFragmentIndex());
        assertFalse(third.hasSplicedBackwardAlignmentLink());
        assertTrue(third.hasSplicedForwardAlignmentLink());
        Assert.assertEquals(3, third.getSplicedForwardAlignmentLink().getFragmentIndex());

        Alignments.AlignmentEntry fourth = reader.next();
        assertEquals(0, fourth.getQueryIndex());
        assertEquals(3, fourth.getFragmentIndex());
        assertTrue(fourth.hasPairAlignmentLink());

        Assert.assertEquals(1, fourth.getPairAlignmentLink().getFragmentIndex());
        assertTrue(fourth.hasSplicedBackwardAlignmentLink());
        Assert.assertEquals(2, fourth.getSplicedBackwardAlignmentLink().getFragmentIndex());

        assertFalse(fourth.hasSplicedForwardAlignmentLink());
    }

    @Test
    public void testSamToCompactTrickCase2() throws IOException {

        SAMToCompactMode importer = new SAMToCompactMode();
        importer.setInputFile("test-data/splicedsamhelper/tricky-spliced-2.sam");
        final String outputFilename = FilenameUtils.concat(BASE_TEST_DIR, "spliced-output-alignment-2");
        importer.setOutputFile(outputFilename);
        importer.execute();

        AlignmentReader reader = new AlignmentReaderImpl(outputFilename);
        assertTrue(reader.hasNext());
        Alignments.AlignmentEntry first = reader.next();

        assertEquals(0, first.getQueryIndex());
        assertEquals(0, first.getFragmentIndex());
        assertTrue(first.hasPairAlignmentLink());
        assertFalse(first.hasSplicedForwardAlignmentLink());
        Assert.assertEquals(1, first.getPairAlignmentLink().getFragmentIndex());
        assertFalse(first.hasSplicedBackwardAlignmentLink());

        Alignments.AlignmentEntry second = reader.next();
        assertEquals(0, second.getQueryIndex());
        assertEquals(1, second.getFragmentIndex());
        assertTrue(second.hasPairAlignmentLink());

        Assert.assertEquals(0, second.getPairAlignmentLink().getFragmentIndex());
        assertFalse(second.hasSplicedForwardAlignmentLink());
        assertFalse(second.hasSplicedBackwardAlignmentLink());


    }

    @Test
    public void testSamToCompactTrickCase3() throws IOException {

        SAMToCompactMode importer = new SAMToCompactMode();
        importer.setInputFile("test-data/splicedsamhelper/tricky-spliced-3.sam");
        final String outputFilename = FilenameUtils.concat(BASE_TEST_DIR, "spliced-output-alignment-3");
        importer.setOutputFile(outputFilename);
        importer.execute();

        AlignmentReader reader = new AlignmentReaderImpl(outputFilename);
        assertTrue(reader.hasNext());
        Alignments.AlignmentEntry first = reader.next();

        assertEquals(0, first.getQueryIndex());
        assertEquals(0, first.getFragmentIndex());
        assertEquals(170769 - 1, first.getPosition());
        assertTrue(first.hasPairAlignmentLink());
        assertFalse(first.hasSplicedForwardAlignmentLink());
        Assert.assertEquals(1, first.getPairAlignmentLink().getFragmentIndex());
        Assert.assertEquals(216048 - 1, first.getPairAlignmentLink().getPosition());
        assertFalse(first.hasSplicedBackwardAlignmentLink());

        Alignments.AlignmentEntry second = reader.next();
        assertEquals(0, second.getQueryIndex());
        assertEquals(1, second.getFragmentIndex());
        assertTrue(second.hasPairAlignmentLink());

        Assert.assertEquals(0, second.getPairAlignmentLink().getFragmentIndex());
        assertTrue(second.hasSplicedForwardAlignmentLink());
        assertFalse(second.hasSplicedBackwardAlignmentLink());


    }

    @Test
    // primary is mapped, but mate is unmapped. Primary must be imported.
    public void testSamToCompactTrickCase4() throws IOException {

        SAMToCompactMode importer = new SAMToCompactMode();
        importer.setInputFile("test-data/splicedsamhelper/tricky-spliced-4.sam");
        final String outputFilename = FilenameUtils.concat(BASE_TEST_DIR, "spliced-output-alignment-4");
        importer.setOutputFile(outputFilename);
        importer.execute();

        AlignmentReader reader = new AlignmentReaderImpl(outputFilename);
        assertTrue(reader.hasNext());
        Alignments.AlignmentEntry first = reader.next();

        assertEquals(0, first.getQueryIndex());
        assertEquals(0, first.getFragmentIndex());
        assertEquals(188966 - 1, first.getPosition());
        assertFalse(first.hasPairAlignmentLink());
        assertFalse(first.hasSplicedForwardAlignmentLink());
        assertFalse(first.hasSplicedBackwardAlignmentLink());


    }

    @Test
    // Test  right soft clip:
    public void testSamToCompactTrickCase5() throws IOException {

        SAMToCompactMode importer = new SAMToCompactMode();
        importer.setInputFile("test-data/splicedsamhelper/tricky-spliced-5.sam");
        final String outputFilename = FilenameUtils.concat(BASE_TEST_DIR, "spliced-output-alignment-5");
        importer.setOutputFile(outputFilename);
        importer.execute();

        AlignmentReader reader = new AlignmentReaderImpl(outputFilename);
        assertTrue(reader.hasNext());
        Alignments.AlignmentEntry first = reader.next();

        assertEquals(0, first.getQueryIndex());
        assertEquals(0, first.getFragmentIndex());
        assertEquals(71428 - 1, first.getPosition());


    }

    @Test
    // Test deletion in the read:
    public void testSamToCompactTrickCase6() throws IOException {

        SAMToCompactMode importer = new SAMToCompactMode();
        importer.setInputFile("test-data/splicedsamhelper/tricky-spliced-6.sam");
        final String outputFilename = FilenameUtils.concat(BASE_TEST_DIR, "spliced-output-alignment-6");
        importer.setOutputFile(outputFilename);
        String[] refs = {"NNNNNNTTAGAAAAACAGAGAGAGAAGGAGAGTAAAGGGAGGAGGCGGAGGAGGAGAAAAGAAGAAAGCAGAGANNNNNN"};

        RandomAccessSequenceTestSupport genomeTestSupport = new RandomAccessSequenceTestSupport(refs);
        importer.setGenome(genomeTestSupport);

        importer.execute();

        AlignmentReader reader = new AlignmentReaderImpl(outputFilename);
        assertTrue(reader.hasNext());
        Alignments.AlignmentEntry first = reader.next();

        assertEquals(0, first.getQueryIndex());
        assertEquals(0, first.getFragmentIndex());
        assertEquals(7 - 1, first.getPosition());


    }

    @Test
    public void testSamToCompactTrickCase7() throws IOException {

        SAMToCompactMode importer = new SAMToCompactMode();
        importer.setInputFile("test-data/splicedsamhelper/tricky-spliced-7.sam");
        final String outputFilename = FilenameUtils.concat(BASE_TEST_DIR, "spliced-output-alignment-7");
        importer.setOutputFile(outputFilename);
        String[] refs = {"NNNNNNTTAGAAAAACAGAGAGAGAAGGAGAGTAAAGGGAGGAGGCGGAGGAGGAGAAAAGAAGAAAGCAGAGANNNNNN"};

        RandomAccessSequenceTestSupport genomeTestSupport = new RandomAccessSequenceTestSupport(refs);
        importer.setGenome(genomeTestSupport);

        importer.execute();

        AlignmentReader reader = new AlignmentReaderImpl(outputFilename);
        assertTrue(reader.hasNext());
        Alignments.AlignmentEntry first = reader.next();

        assertEquals(0, first.getQueryIndex());
        assertEquals(0, first.getFragmentIndex());
        assertEquals(8 - 1, first.getPosition());


    }

    @Test
    public void testSamToCompactTrickCase8() throws IOException {

        SAMToCompactMode importer = new SAMToCompactMode();
        importer.setInputFile("test-data/splicedsamhelper/tricky-spliced-8.sam");
        final String outputFilename = FilenameUtils.concat(BASE_TEST_DIR, "spliced-output-alignment-8");
        importer.setOutputFile(outputFilename);
        MutableString seq = new MutableString();
        for (int i = 0; i < 194407 - 7; i++) {
            seq.append('N');
        }
        seq.append("NNNNNNTTAGAAAAACAGAGAGAGAAGGAGAGTAAAGGGAGGAGGCGGAGGAGGAGAAAAGAAGAAAGCAGAGANNNNNN");
        String[] refs = {seq.toString()};

        RandomAccessSequenceTestSupport genomeTestSupport = new RandomAccessSequenceTestSupport(refs);
        importer.setGenome(genomeTestSupport);

        importer.execute();

        AlignmentReader reader = new AlignmentReaderImpl(outputFilename);
        assertTrue(reader.hasNext());
        Alignments.AlignmentEntry first = reader.next();

        assertEquals(0, first.getQueryIndex());
        assertEquals(0, first.getFragmentIndex());
        assertEquals(194301 - 1, first.getPosition());


    }

    @Test
    public void testSamToCompactTrickCase9() throws IOException {

        SAMToCompactMode importer = new SAMToCompactMode();
        importer.setInputFile("test-data/splicedsamhelper/tricky-spliced-9.sam");
        final String outputFilename = FilenameUtils.concat(BASE_TEST_DIR, "spliced-output-alignment-9");
        importer.setOutputFile(outputFilename);
        MutableString seq = new MutableString();

        seq.append("NNNTTAGAAAAACAGAGAGAGAAGGAGAGTAAAGGGAGGAGGCGGAGGAGGAGAAAAGAAGAAAGCAGAGANNNNNN");
        for (int i = 0; i < 573; i++) {
            seq.insert(25, '-');
        }
        String[] refs = {seq.toString()};

        RandomAccessSequenceTestSupport genomeTestSupport = new RandomAccessSequenceTestSupport(refs);
        importer.setGenome(genomeTestSupport);

        importer.execute();

        AlignmentReader reader = new AlignmentReaderImpl(outputFilename);
        assertTrue(reader.hasNext());
        Alignments.AlignmentEntry first = reader.next();

        assertEquals(0, first.getQueryIndex());
        assertEquals(0, first.getFragmentIndex());
        assertEquals(3 - 1, first.getPosition());


    }

    @Test
    // like 9 no genome
    public void testSamToCompactTrickCase9NoGenome() throws IOException {

        SAMToCompactMode importer = new SAMToCompactMode();
        importer.setInputFile("test-data/splicedsamhelper/tricky-spliced-9.sam");
        final String outputFilename = FilenameUtils.concat(BASE_TEST_DIR, "spliced-output-alignment-9");
        importer.setOutputFile(outputFilename);
        importer.execute();

        AlignmentReader reader = new AlignmentReaderImpl(outputFilename);
        assertTrue(reader.hasNext());
        Alignments.AlignmentEntry first = reader.next();

        assertEquals(0, first.getQueryIndex());
        assertEquals(0, first.getFragmentIndex());
        assertEquals(3 - 1, first.getPosition());


    }
    /*
PATHBIO-SOLEXA2:2:37:931:1658#0	97	chr10	97392943	255	11M10083N29M	=	64636105	0	CTGGATACAATGAGATCTGAAGACGGTTGTACACTTGACC	BBB@BA???BBCBA@>AA=6>A@?B?<<B<>;=@ABA@@<	NM:i:0	XS:A:-	NS:i:0
PATHBIO-SOLEXA2:2:37:931:1658#0	145	chr11	64636105	255	11M447N29M	=	97392943	0	AGGGCCCCTGGGCGCCGCGGCTCTGCTACTGCTGCTGCCC	A?@AB@?@?=A=AAA@<A<BBBBBBBBBA@<@<BBB@BAA	NM:i:0	XS:A:+	NS:i:0

    */
}
