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

package edu.cornell.med.icb.goby.readers.sam;

import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.alignments.PerQueryAlignmentData;
import edu.cornell.med.icb.goby.alignments.TestIteratedSortedAlignment2;
import edu.cornell.med.icb.goby.modes.CompactToSAMMode;
import edu.cornell.med.icb.goby.modes.SAMToCompactMode;
import edu.cornell.med.icb.goby.modes.SamHelper;
import edu.cornell.med.icb.goby.reads.DualRandomAccessSequenceCache;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceTestSupport;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.logging.ProgressLogger;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;
import org.junit.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.Map;

import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.fail;

/**
 * Test SamRecordParser.
 */
public class TestSamRecordParser {

    private static final Logger LOG = Logger.getLogger(TestSamRecordParser.class);

    private SamHelper globalGamHelper = new SamHelper();

    private static final String BASE_TEST_DIR = "test-results/splicedsamhelper";

    //
    //  testSamToCompactTrickCase1-3 fails because this the sam reference builder requires an MD:Z tag.
    //


    @Test
    // like 9 no genome
    public void testSamToCompactTrickCase9NoGenome() throws IOException {
        final String inputFile = "test-data/splicedsamhelper/tricky-spliced-9.sam";
        final SAMFileReader parser = new SAMFileReader(new FileInputStream(inputFile));
        parser.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        final SamRecordParser recordParser = new SamRecordParser();
        for (final SAMRecord samRecord : new SAMRecordIterable(parser.iterator())) {

            final GobySamRecord gobySamRecord = recordParser.processRead(samRecord);

            assertEquals("Incorrect number of segments", 2, gobySamRecord.getNumSegments());
            GobySamSegment first = gobySamRecord.getSegment(0);

            assertEquals(3 - 1, first.getPosition());
        }
    }

    @Test
    // variation after splice
    public void testSamToCompactTrickCase10NoGenome() throws IOException {

        final String inputFile = "test-data/splicedsamhelper/tricky-spliced-10.sam";
        final SAMFileReader parser = new SAMFileReader(new FileInputStream(inputFile));
        parser.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        final SamRecordParser recordParser = new SamRecordParser();
        for (final SAMRecord samRecord : new SAMRecordIterable(parser.iterator())) {

            final GobySamRecord gobySamRecord = recordParser.processRead(samRecord);

            assertEquals("Incorrect number of segments", 2, gobySamRecord.getNumSegments());
            GobySamSegment first = gobySamRecord.getSegment(0);

            if (gobySamRecord.readNum == 0) {
                assertEquals(3 - 1, first.getPosition());
            } else if (gobySamRecord.readNum == 1) {
                assertEquals(3 - 1, first.getPosition());
            }
        }
    }

    @Test
    public void testSamToCompactTrickCase11() throws IOException {
        final String inputFile = "test-data/splicedsamhelper/tricky-spliced-11.sam";
        final SAMFileReader parser = new SAMFileReader(new FileInputStream(inputFile));
        parser.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        final SamRecordParser recordParser = new SamRecordParser();
        for (final SAMRecord samRecord : new SAMRecordIterable(parser.iterator())) {

            final GobySamRecord gobySamRecord = recordParser.processRead(samRecord);

            assertEquals("Incorrect number of segments", 2, gobySamRecord.getNumSegments());
            GobySamSegment first = gobySamRecord.getSegment(0);
            assertEquals(26800015 - 1, first.getPosition());
        }
    }


    @Test
    // like 9 no genome
    public void testSamToCompactTrickCase12NoGenome() throws IOException {

        final String inputFile = "test-data/splicedsamhelper/tricky-spliced-12.sam";
        final SAMFileReader parser = new SAMFileReader(new FileInputStream(inputFile));
        final SamRecordParser recordParser = new SamRecordParser();
        for (final SAMRecord samRecord : new SAMRecordIterable(parser.iterator())) {

            final GobySamRecord gobySamRecord = recordParser.processRead(samRecord);

            assertEquals("Incorrect number of segments", 2, gobySamRecord.getNumSegments());
            GobySamSegment first = gobySamRecord.getSegment(0);
            GobySamSegment second = gobySamRecord.getSegment(1);

            assertEquals(15013, first.getPosition());
            assertEquals(3, first.getQueryPosition());
            assertEquals(25, first.getQueryAlignedLength());

            assertEquals(15795, second.getPosition());
            assertEquals(28, second.getQueryPosition());
            assertEquals(7, second.getQueryAlignedLength());
        }
    }

    @Test
    // To test import of soft clips:
    public void testSamToCompactTrickCase13NoGenomeSoftClips() throws IOException {

        final String inputFile = "test-data/splicedsamhelper/tricky-spliced-13.sam";
        final SAMFileReader parser = new SAMFileReader(new FileInputStream(inputFile));
        final SamRecordParser recordParser = new SamRecordParser();
        for (final SAMRecord samRecord : new SAMRecordIterable(parser.iterator())) {

            final GobySamRecord gobySamRecord = recordParser.processRead(samRecord);

            assertEquals("Incorrect number of segments", 2, gobySamRecord.getNumSegments());
            GobySamSegment first = gobySamRecord.getSegment(0);
            GobySamSegment second = gobySamRecord.getSegment(1);

            assertEquals(15013, first.getPosition());
            assertEquals(3, first.getQueryPosition());
            assertEquals("AAT", first.getSoftClippedBasesLeft());
            assertEquals("", first.getSoftClippedBasesRight());
            assertEquals(25, first.getQueryAlignedLength());

            assertEquals(15795, second.getPosition());
            assertEquals(28, second.getQueryPosition());
            assertEquals(5, second.getQueryAlignedLength());
            assertEquals("", second.getSoftClippedBasesLeft());
            assertEquals("TC", second.getSoftClippedBasesRight());
        }
    }

    @Test
    // To test import of soft clips:
    // NOTE: as this doesn't supply genome, it make alter how the test was intended.
    //       I had to change the softclip. Just adding genome doesn't work because
    //       of how SamRecordParser.allRefBases is constructed to including
    //       splices, etc.
    public void testSamToCompactTrickCase13SoftClipsWithGenome() throws IOException {

        final String inputFile = "test-data/splicedsamhelper/tricky-spliced-14.sam";
        final SAMFileReader parser = new SAMFileReader(new FileInputStream(inputFile));
        parser.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        final SamRecordParser recordParser = new SamRecordParser();

        for (final SAMRecord samRecord : new SAMRecordIterable(parser.iterator())) {
            final GobySamRecord gobySamRecord = recordParser.processRead(samRecord);

            assertEquals("Incorrect number of segments", 2, gobySamRecord.getNumSegments());
            GobySamSegment first = gobySamRecord.getSegment(0);
            GobySamSegment second = gobySamRecord.getSegment(1);

            assertEquals(4, first.getPosition());
            assertEquals(3, first.getQueryPosition());
            assertEquals("AAT", first.getSoftClippedBasesLeft());  // Prev test A=T
            assertEquals("", first.getSoftClippedBasesRight());
            assertEquals(25, first.getQueryAlignedLength());

            assertEquals(40 - 1, second.getPosition());
            assertEquals(28, second.getQueryPosition());
            assertEquals(5, second.getQueryAlignedLength());
            assertEquals("TC", second.getSoftClippedBasesRight());
            assertEquals("", second.getSoftClippedBasesLeft());
        }
    }

    @Test
    public void testSamToCompactTrickCase15NoGenomeThreeSplice() throws IOException {

        final String inputFile = "test-data/splicedsamhelper/tricky-spliced-15.sam";
        final SAMFileReader parser = new SAMFileReader(new FileInputStream(inputFile));
        final SamRecordParser recordParser = new SamRecordParser();
        for (final SAMRecord samRecord : new SAMRecordIterable(parser.iterator())) {
            final GobySamRecord gobySamRecord = recordParser.processRead(samRecord);

            assertEquals("Incorrect number of segments", 3, gobySamRecord.getNumSegments());

            final GobySamSegment first = gobySamRecord.getSegment(0);
            final GobySamSegment second = gobySamRecord.getSegment(1);
            final GobySamSegment third = gobySamRecord.getSegment(2);

            assertEquals(32485524 - 1, first.position);
            assertEquals(2, first.queryPosition);
            assertEquals("CG", first.softClippedBasesLeft.toString());
            assertEquals("", first.softClippedBasesRight.toString());
            assertEquals(6, first.queryAlignedLength);
            assertEquals(6, first.targetAlignedLength);

            assertEquals(32485524 + 6 + 301 - 1, second.position);
            assertEquals(8, second.queryPosition);
            assertEquals("", second.softClippedBasesLeft.toString());
            assertEquals("", second.softClippedBasesRight.toString());
            assertEquals(24, second.queryAlignedLength);
            assertEquals(24, second.targetAlignedLength);

            assertEquals(32485524 + 6 + 301 + 24 + 478 - 1, third.position);
            assertEquals(32, third.queryPosition);
            assertEquals("", third.softClippedBasesLeft.toString());
            assertEquals("", third.softClippedBasesRight.toString());
            assertEquals(3, third.queryAlignedLength);
            assertEquals(3, third.targetAlignedLength);
        }
    }

    @Test
    public void testSamToCompactTrickCase16() throws IOException {

        final String inputFile = "test-data/splicedsamhelper/tricky-spliced-16.sam";
        final SAMFileReader parser = new SAMFileReader(new FileInputStream(inputFile));
        parser.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        final SamRecordParser recordParser = new SamRecordParser();
        for (final SAMRecord samRecord : new SAMRecordIterable(parser.iterator())) {
            final GobySamRecord gobySamRecord = recordParser.processRead(samRecord);

            assertEquals("Incorrect number of segments", 1, gobySamRecord.getNumSegments());
            final GobySamSegment segment = gobySamRecord.getSegment(0);
            if (gobySamRecord.readNum == 0) {
                final GobySamSegment first = segment;
                assertEquals(190077 - 1, first.getPosition());
                assertEquals(13, first.getQueryPosition());
                assertEquals("AGTGGCAGCACGA", first.getSoftClippedBasesLeft());
                assertEquals("TGCT", first.getSoftClippedBasesRight());
                assertEquals(51, first.getQueryAlignedLength());
                assertEquals(50, first.getTargetAlignedLength());

                assertEquals(5, first.getSequenceVariationsCount());

                GobyQuickSeqvar seqvar = first.getSequenceVariations(0);
                assertEquals("T", seqvar.getFrom());
                assertEquals("C", seqvar.getTo());
                assertEquals(22, seqvar.getReadIndex());
                assertEquals(9, seqvar.getPosition());
                assertArrayEquals(byteArray(3), seqvar.getToQualitiesAsBytes());  // 8 returned in the old test, 3 is right

                seqvar = first.getSequenceVariations(1);
                assertEquals("T", seqvar.getFrom());
                assertEquals("G", seqvar.getTo());    // Original test said A
                assertEquals(26, seqvar.getReadIndex());
                assertEquals(13, seqvar.getPosition());
                assertArrayEquals(byteArray(34), seqvar.getToQualitiesAsBytes());  // 24 in the old test, this is right

                seqvar = first.getSequenceVariations(2);
                assertEquals("C", seqvar.getFrom());
                assertEquals("A", seqvar.getTo());
                assertEquals(33, seqvar.getReadIndex());
                assertEquals(20, seqvar.getPosition());
                assertArrayEquals(byteArray(28), seqvar.getToQualitiesAsBytes());  // 14 in the old test

                seqvar = first.getSequenceVariations(3);
                assertEquals("-", seqvar.getFrom());
                assertEquals("C", seqvar.getTo());   // Original test said CG
                assertEquals(35, seqvar.getReadIndex());
                assertEquals(21, seqvar.getPosition());
                assertArrayEquals(byteArray(27), seqvar.getToQualitiesAsBytes());  // 3, 15 in the old test

                seqvar = first.getSequenceVariations(4);
                assertEquals("T", seqvar.getFrom());
                assertEquals("A", seqvar.getTo());   // Original test said CG
                assertEquals(36, seqvar.getReadIndex());
                assertEquals(22, seqvar.getPosition());
                assertArrayEquals(byteArray(35), seqvar.getToQualitiesAsBytes());  // 3, 15 in the old test

            } else if (gobySamRecord.readNum == 1) {
                final GobySamSegment second = segment;
                //second's CIGAR is 20S48M
                assertEquals(190246 - 1, second.getPosition());
                assertEquals(20, second.getQueryPosition());
                assertEquals("CAGTGTCGTGGCTGCACGCC", second.getSoftClippedBasesLeft());
                assertEquals("", second.getSoftClippedBasesRight());
                assertEquals(48, second.getQueryAlignedLength());
                assertEquals(48, second.getTargetAlignedLength());

                assertEquals(1, second.getSequenceVariationsCount());

                GobyQuickSeqvar seqvar = second.getSequenceVariations(0);
                assertEquals("A", seqvar.getFrom());
                assertEquals("C", seqvar.getTo());
                assertEquals(18, seqvar.getReadIndex());
                assertEquals(31, seqvar.getPosition());
                assertArrayEquals(byteArray(35), seqvar.getToQualitiesAsBytes());
            }
        }
    }

    /**
     * This comes from
     * Quality score difference.
     * See XAAOBVT  [Open Session]
     *   chr1:45,881,903-45,881,942
     *   T variation with qual score 2 should be 36
     * Sample = XAAOBVT
     * Read group = 1
     * ----------------------
     * Read name = 509.6.68.19057.157284
     * @throws IOException
     */
    @Test
    public void testSamToCompactTrickCase17() throws IOException {
        final String inputFile = "test-data/splicedsamhelper/tricky-spliced-17.sam";
        final SAMFileReader parser = new SAMFileReader(new FileInputStream(inputFile));
        parser.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        final SamRecordParser recordParser = new SamRecordParser();
        for (final SAMRecord samRecord : new SAMRecordIterable(parser.iterator())) {
            final GobySamRecord gobySamRecord = recordParser.processRead(samRecord);

            assertEquals("Incorrect number of segments", 1, gobySamRecord.getNumSegments());

            final GobySamSegment first = gobySamRecord.getSegment(0);
            
            assertEquals(45881869 - 1, first.getPosition());

            //509.6.68.19057.157284	83	chr1	45881869	29	6S23M1I6M1D16M1I47M	=	45881519	-443	TTACCCGCTTTCCTTGCCCAAATTTTAAGTTTCNGGAAAAGGGGAGGGAAATGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGTGACAGAGTGTCAC	#######################################################ECGGGGGGGGGGGGGGGGGGGGGGGGHHHHHHHHHHHHHHHHHHH	MD:Z:3G4A7C4C1A2T2^T2C4T3A0G5C44	RG:Z:1	XG:i:3	AM:i:29	NM:i:14	SM:i:29	XM:i:10	XO:i:3	XT:A:M

            assertEquals(6, first.getQueryPosition());
            assertEquals("TTACCC", first.getSoftClippedBasesLeft());
            assertEquals("", first.getSoftClippedBasesRight());

            assertEquals(13, first.getSequenceVariationsCount());

            GobyQuickSeqvar last = first.getSequenceVariations(12);
            assertEquals("T", last.getTo());
            assertEquals("C", last.getFrom());
            assertArrayEquals(byteArray(36), last.getToQualitiesAsBytes());
        }
    }

    /**
     * Same test as above, but write the sam to compact and then read the sequence variation form the compact.
     * @throws IOException
     */
    @Test
    public void testSamToCompactTrickCase17ViaWriter() throws IOException {
        SAMToCompactMode importer = new SAMToCompactMode();
        importer.setInputFile("test-data/splicedsamhelper/tricky-spliced-17.sam");
        final String outputFilename = FilenameUtils.concat(BASE_TEST_DIR, "tricky-spliced-17");
        importer.setPreserveSoftClips(true);
        importer.setOutputFile(outputFilename);
        importer.execute();

        final AlignmentReader reader = new AlignmentReaderImpl(outputFilename);
        assertTrue(reader.hasNext());
        final Alignments.AlignmentEntry first = reader.next();


        assertEquals(45881869 - 1, first.getPosition());

        //509.6.68.19057.157284	83	chr1	45881869	29	6S23M1I6M1D16M1I47M	=	45881519	-443	TTACCCGCTTTCCTTGCCCAAATTTTAAGTTTCNGGAAAAGGGGAGGGAAATGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGTGACAGAGTGTCAC	#######################################################ECGGGGGGGGGGGGGGGGGGGGGGGGHHHHHHHHHHHHHHHHHHH	MD:Z:3G4A7C4C1A2T2^T2C4T3A0G5C44	RG:Z:1	XG:i:3	AM:i:29	NM:i:14	SM:i:29	XM:i:10	XO:i:3	XT:A:M

        assertEquals(6, first.getQueryPosition());
        assertEquals("TTACCC", first.getSoftClippedBasesLeft());
        assertEquals("", first.getSoftClippedBasesRight());

        assertEquals(13, first.getSequenceVariationsCount());

        Alignments.SequenceVariation last = first.getSequenceVariations(12);
        assertEquals("T", last.getTo());
        assertEquals("C", last.getFrom());
        assertArrayEquals(byteArray(36), last.getToQuality().toByteArray());
    }

    /**
     * Find the local celegans random genome from an array of options.
     * @return the local celegans random genome.
     */
    private String findCelegansGenome() {
        final String[] dirs = {
                "/tmp/celegans",
                "/scratchLocal/gobyweb/input-data/reference-db/goby-benchmark-paper/cElegans",
                "/home/ccontrol/goby-data/celegans" };
        for (final String dir : dirs) {
            final String testRootFilename = dir + "/" + "random-access-genome";
            final String testFilename = testRootFilename + ".names";
            System.out.println("Looking for :" + testFilename);
            final File testFile = new File(testFilename);
            if (testFile.exists()) {
                return testRootFilename;
            }
        }
        return null;
    }


    /**
     * Find the local 1000g random genome from an array of options.
     * @return the local celegans random genome.
     */
    private String findThousandGenome() {
        final String[] dirs = {
                "/tmp/1000g",
                "/scratchLocal/gobyweb/input-data/reference-db/1000GENOMES.37/homo_sapiens/reference",
                "/home/ccontrol/goby-data/1000g-random-access" };
        for (final String dir : dirs) {
            final String testRootFilename = dir + "/" + "random-access-genome";
            final String testFilename = testRootFilename + ".names";
            System.out.println("Looking for :" + testFilename);
            final File testFile = new File(testFilename);
            if (testFile.exists()) {
                return testRootFilename;
            }
        }
        return null;
    }

    /**
     * Test that DOES fail, for local testing, not for server testing.
     * The 1M doesn't exist on the testing server.
     * @throws IOException error
     */
    // @Test
    public void testRoundTripFail() throws IOException {
        final RoundTripConfig rtc = new RoundTripConfig();
        rtc.inputGenomeFilename = findThousandGenome();
        rtc.sourceBamFilename = "test-data/splicedsamhelper/HZFWPTI-first-500.sam";
        rtc.destGobyBasename = FilenameUtils.concat(BASE_TEST_DIR, "1M");
        rtc.destBamFilename = FilenameUtils.concat(BASE_TEST_DIR, "1M.bam");
        rtc.convertBamToGoby = false;
        rtc.convertGobyToBam = false;
        testRoundTripAny(rtc);
    }

    /**
     * Test round trip of tricky-spliced-18.sam that input and output
     * sam compare.
     * @throws IOException error
     */
    @Test
    public void testRoundTripTrickySpliced18() throws IOException {
        final RoundTripConfig rtc = new RoundTripConfig();
        rtc.inputGenomeFilename = findThousandGenome();
        rtc.sourceBamFilename = "test-data/splicedsamhelper/tricky-spliced-18.sam";
        rtc.destGobyBasename = FilenameUtils.concat(BASE_TEST_DIR, "tricky-spliced-18");
        rtc.destBamFilename = FilenameUtils.concat(BASE_TEST_DIR, "tricky-spliced-18.sam");
        testRoundTripAny(rtc);
        rtc.keepQualityScores = false;
        testRoundTripAny(rtc);
        rtc.keepSoftClips = false;
        testRoundTripAny(rtc);
    }

    /**
     * The first 500 alignments of HZ. Round trip test.
     * @throws IOException error
     */
    @Test
    public void testRoundTripHZFirst500() throws IOException {
        final RoundTripConfig rtc = new RoundTripConfig();
        rtc.inputGenomeFilename = findThousandGenome();
        rtc.sourceBamFilename = "test-data/splicedsamhelper/HZFWPTI-first-500.sam";
        rtc.destGobyBasename = FilenameUtils.concat(BASE_TEST_DIR, "HZFWPTI-first-500");
        rtc.destBamFilename = FilenameUtils.concat(BASE_TEST_DIR, "HZFWPTI-first-500.sam");
        testRoundTripAny(rtc);
        rtc.keepQualityScores = false;
        testRoundTripAny(rtc);
        rtc.keepSoftClips = false;
        testRoundTripAny(rtc);
    }

    /**
     * The first 1M alignments of UAN. This large dataset does not exist on the testing server.
     * @throws IOException error
     */
    // @Test
    public void testRoundTrip1M() throws IOException {
        final RoundTripConfig rtc = new RoundTripConfig();
        rtc.inputGenomeFilename = findThousandGenome();
        rtc.sourceBamFilename = "test-data/splicedsamhelper/1M.bam";
        rtc.destGobyBasename = FilenameUtils.concat(BASE_TEST_DIR, "1M");
        rtc.destBamFilename = FilenameUtils.concat(BASE_TEST_DIR, "1M.bam");
        rtc.canonicalMdzForComparison = false;
        testRoundTripAny(rtc);
    }

    /**
     * The whole of HENGLIT. This large dataset does not exist on the testing server.
     * @throws IOException error
     */
    // @Test
    public void testRoundTripHenglit() throws IOException {
        final RoundTripConfig rtc = new RoundTripConfig();
        rtc.inputGenomeFilename = findCelegansGenome();
        rtc.sourceBamFilename = "/tmp/HENGLIT.bam";
        rtc.destGobyBasename = "/tmp/HENGLIT-to-goby";
        rtc.destBamFilename = "/tmp/HENGLIT-from-goby.bam";
        // rtc.convertBamToGoby = false;
        // rtc.convertGobyToBam = false;
        rtc.canonicalMdzForComparison = false;
        testRoundTripAny(rtc);
    }

    /**
     * Configuration for round trip comparison.
     */
    private static class RoundTripConfig {
        String inputGenomeFilename;
        String sourceBamFilename;
        String destGobyBasename;
        String destBamFilename;
        boolean convertBamToGoby = true;
        boolean convertGobyToBam = true;
        boolean keepQualityScores = true;
        boolean keepSoftClips = true;
        boolean stopAtOneFailure = false;
        boolean canonicalMdzForComparison = true;
    }

    public void testRoundTripAny(final RoundTripConfig rtc) throws IOException {

        // IMPORTANT!!
        // ** These two should always be set to true unless you are doing MANUAL testing and want to not do
        // ** one or both of the conversions.

        if (!new File(rtc.destGobyBasename + ".entries").exists()) {
            rtc.convertBamToGoby = true;
        }
        if (!new File(rtc.destBamFilename).exists()) {
            rtc.convertGobyToBam = true;
        }

        assertNotNull("Could not locate random-access-genome in specified locations", rtc.inputGenomeFilename);
        RandomAccessSequenceInterface genome = null;
        if (rtc.convertBamToGoby || rtc.convertGobyToBam) {
            genome = new DualRandomAccessSequenceCache();
            try {
                ((DualRandomAccessSequenceCache)genome).load(rtc.inputGenomeFilename);
            } catch (ClassNotFoundException e) {
                throw new IOException("Could not load genome", e);
            }
        }

        if (rtc.convertBamToGoby) {
            LOG.info("Converting bam to compact alignment");
            final SAMToCompactMode importer = new SAMToCompactMode();
            importer.setInputFile(rtc.sourceBamFilename);
            importer.setPreserveSoftClips(rtc.keepSoftClips);
            importer.setPreserveAllTags(true);
            importer.setOutputFile(rtc.destGobyBasename);
            importer.setGenome(genome);
            importer.setPreserveReadQualityScores(rtc.keepQualityScores);
            importer.execute();
        }

        if (rtc.convertGobyToBam) {
            LOG.info("Converting compact alignment to bam");
            final CompactToSAMMode exporter = new CompactToSAMMode();
            exporter.setGenome(genome);
            exporter.setInputBasename(rtc.destGobyBasename);
            exporter.setOutput(rtc.destBamFilename);
            exporter.execute();
        }

        LOG.info("Comparing source bam and destination bam");
        final SAMFileReader sourceParser = new SAMFileReader(new FileInputStream(rtc.sourceBamFilename));
        final SAMFileReader destParser = new SAMFileReader(new FileInputStream(rtc.destBamFilename));
        // We need to set the validation to silent because an incomplete file (if the source isn't the entire file)
        // we can see errors that wouldn't exist in a real conversion.
        sourceParser.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        destParser.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        final SAMRecordIterator sourceIterator = sourceParser.iterator();
        final SAMRecordIterator destIterator = destParser.iterator();
        AlignmentReaderImpl gobyReader = null;
        if (rtc.destGobyBasename != null) {
            gobyReader = new AlignmentReaderImpl(rtc.destGobyBasename);
        }
        final ProgressLogger progress = new ProgressLogger(LOG);
        progress.displayFreeMemory = true;
        final SamComparison samComparison = new SamComparison();
        progress.start();
        samComparison.setMappedQualitiesPreserved(rtc.keepQualityScores);
        samComparison.setSoftClipsPreserved(rtc.keepSoftClips);
        samComparison.setCheckMate(false);
        samComparison.setCanonicalMdzForComparison(rtc.canonicalMdzForComparison);
        samComparison.reset();
        while (sourceIterator.hasNext()) {
            samComparison.setExpectedSamRecord(sourceIterator.next());
            if (samComparison.getExpectedSamRecord().getReadUnmappedFlag()) {
                // We don't store unmapped reads, so skip this source record
                continue;
            }
            assertTrue("Not enough records in destination SAM/BAM file", destIterator.hasNext());
            samComparison.setActualSamRecord(destIterator.next());
            if (gobyReader != null) {
                assertTrue("Not enough records in goby compact-alignment file", gobyReader.hasNext());
                samComparison.setGobyAlignment(gobyReader.next());
            }
            if (rtc.stopAtOneFailure) {
                assertTrue("sam comparison failed", samComparison.compareSamRecords());
            } else {
                samComparison.compareSamRecords();
            }
            progress.lightUpdate();
        }
        progress.stop();
        if (!rtc.stopAtOneFailure && samComparison.getComparisonFailureCount() > 0) {
            fail("Number of comparison failures: " + samComparison.getComparisonFailureCount());
        }
    }

    @Test
    public void testDelNonSplice1() throws IOException {
        final String inputFile = "test-data/splicedsamhelper/del-nonsplice-1.sam";
        final SAMFileReader parser = new SAMFileReader(new FileInputStream(inputFile));
        parser.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        final SamRecordParser recordParser = new SamRecordParser();
        for (final SAMRecord samRecord : new SAMRecordIterable(parser.iterator())) {
            final GobySamRecord gobySamRecord = recordParser.processRead(samRecord);

            assertEquals("Incorrect number of segments", 1, gobySamRecord.getNumSegments());
            final GobySamSegment segment = gobySamRecord.getSegment(0);
            if (gobySamRecord.readNum == 0) {
                final GobySamSegment first = segment;
                assertEquals(31 - 1, first.getPosition());
                assertEquals(0, first.getQueryPosition());
                assertEquals("", first.getSoftClippedBasesLeft());
                assertEquals("", first.getSoftClippedBasesRight());
                assertEquals(47, first.getQueryAlignedLength());
                assertEquals(50, first.getTargetAlignedLength());

                assertEquals(1, first.getSequenceVariationsCount());

                GobyQuickSeqvar seqvar = first.getSequenceVariations(0);
                assertEquals("TCC", seqvar.getFrom());
                assertEquals("---", seqvar.getTo());
                assertEquals(26, seqvar.getReadIndex());
                assertEquals(22, seqvar.getPosition());
                assertEquals(0, seqvar.getToQualitiesAsBytes().length);
            }
        }
    }

    @Test
    public void testLeftPadding1() throws IOException {
        final String inputFile = "test-data/splicedsamhelper/leftpad-nosplice-1.sam";
        final SAMFileReader parser = new SAMFileReader(new FileInputStream(inputFile));
        parser.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        final SamRecordParser recordParser = new SamRecordParser();
        for (final SAMRecord samRecord : new SAMRecordIterable(parser.iterator())) {
            final GobySamRecord gobySamRecord = recordParser.processRead(samRecord);

            assertEquals("Incorrect number of segments", 1, gobySamRecord.getNumSegments());
            final GobySamSegment segment = gobySamRecord.getSegment(0);
            if (gobySamRecord.readNum == 0) {
                final GobySamSegment first = segment;
                assertEquals(6 - 1, first.getPosition());
                assertEquals(5, first.getQueryPosition());
                assertEquals("AAAAA", first.getSoftClippedBasesLeft());
                assertEquals("", first.getSoftClippedBasesRight());
                assertEquals(45, first.getQueryAlignedLength());
                assertEquals(45, first.getTargetAlignedLength());

                assertEquals(1, first.getSequenceVariationsCount());

                GobyQuickSeqvar seqvar = first.getSequenceVariations(0);
                assertEquals("A", seqvar.getFrom());
                assertEquals("G", seqvar.getTo());
                assertEquals(20, seqvar.getReadIndex());
                assertEquals(15, seqvar.getPosition());
                assertEquals(1, seqvar.getToQualitiesAsBytes().length);
                assertArrayEquals(byteArray(20), seqvar.getToQualitiesAsBytes());
            }
        }
    }

    /**
     * Compare the gsnap->sam created sequence variations we find when we parse the sam file with
     * SamRecordParser exactly match those generated with gsnap->compactAlignment->displaySequenceVariations(per-base).
     * @throws IOException
     */
    @Test
    public void testSeqVarReads() throws IOException {
        final String inputFile = "test-data/seq-var-test/seq-var-reads-gsnap.sam";
        final SAMFileReader parser = new SAMFileReader(new FileInputStream(inputFile));
        parser.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        final SamRecordParser recordParser = new SamRecordParser();
        final Int2ObjectMap<PerQueryAlignmentData> seqvarDataMap =  TestIteratedSortedAlignment2.readSeqVarFile(
                "test-data/seq-var-test/seq-var-reads-gsnap.seqvar");
        final int[] seqvarQueryIndexes = seqvarDataMap.keySet().toIntArray();
        Arrays.sort(seqvarQueryIndexes);

        final Int2ObjectMap<PerQueryAlignmentData> samSeqvarDataMap = new Int2ObjectOpenHashMap<PerQueryAlignmentData>();

        for (final SAMRecord samRecord : new SAMRecordIterable(parser.iterator())) {
            final GobySamRecord gobySamRecord = recordParser.processRead(samRecord);
            if (gobySamRecord == null) {
                continue;
            }
            if (gobySamRecord.getSequenceVariationsCount() == 0) {
                continue;
            }
            final GobySamSegment segment = gobySamRecord.getSegment(0);

            final int queryIndex = gobySamRecord.getReadNum();
            for (final GobyQuickSeqvar seqvar : gobySamRecord.getSequenceVariations()) {

                // convert variation position to position on the reference:
                final int positionOnReference = segment.getPosition() + seqvar.getPosition();
                final int readIndex = seqvar.getReadIndex();
                final String from = seqvar.getFrom();
                final String to = seqvar.getTo();
                final int fromLength = from.length();

                int fromOffset = 0;
                int toOffset = 0;
                final int readIndexIncrementValue = (gobySamRecord.isReverseStrand() ? -1 : 1);
                for (int i = 0; i < fromLength; i++) {
                    final char fromChar = from.charAt(i);
                    final char toChar = to.charAt(i);

                    PerQueryAlignmentData variationsForIndex = samSeqvarDataMap.get(queryIndex);
                    if (variationsForIndex == null) {
                        variationsForIndex = new PerQueryAlignmentData();
                        variationsForIndex.queryPosition = segment.getSoftClippedBasesLeft().length();
                        variationsForIndex.reverseStrand = gobySamRecord.isReverseStrand();
                        samSeqvarDataMap.put(queryIndex, variationsForIndex);
                    }
                    variationsForIndex.observe(positionOnReference + fromOffset, readIndex + toOffset, fromChar, toChar);
                    if (fromChar != '-') {
                        fromOffset += 1;
                    }
                    if (toChar != '-') {
                        toOffset += readIndexIncrementValue;
                    }
                }
            }
        }

        final int[] samSeqvarQueryIndexes = samSeqvarDataMap.keySet().toIntArray();
        Arrays.sort(samSeqvarQueryIndexes);

        /*
        System.out.println("Sequence variations from .tsv file");
        dumpSeqVars(seqvarQueryIndexes, seqvarDataMap);
        System.out.println();
        System.out.println("Sequence variations from .sam file");
        dumpSeqVars(seqvarQueryIndexes, samSeqvarDataMap);
        */
        verifySequenceVariationsMatch(seqvarQueryIndexes, samSeqvarQueryIndexes, seqvarDataMap, samSeqvarDataMap);
    }

    public static void dumpSeqVars(int[] seqvarQueryIndexes, final Int2ObjectMap<PerQueryAlignmentData> seqvarDataMap) {
        for (int queryIndex : seqvarQueryIndexes) {
            final PerQueryAlignmentData align = seqvarDataMap.get(queryIndex);
            System.out.println("queryIndex=" + queryIndex);
            System.out.println(align.toString());
        }
    }

    /**
     * Compare sequence varations.
     * @param seqvarQueryIndexes the sorted query indexes for the seqvarDataMap
     * @param alignmentQueryIndexes the sorted query indexes for the alignmentDataMap
     * @param seqvarDataMap the query index to PerQueryAlignmentData for the "expected" values
     * @param alignmentDataMap the query index to PerQueryAlignmentData for the "actual" values
     * @throws IOException error reading input files
     */
    public static void verifySequenceVariationsMatch(
            final int[] seqvarQueryIndexes,
            final int[] alignmentQueryIndexes,
            final Int2ObjectMap<PerQueryAlignmentData> seqvarDataMap,
            final Int2ObjectMap<PerQueryAlignmentData> alignmentDataMap) throws IOException {
        assertArrayEquals("queryIndexes are not the same", seqvarQueryIndexes, alignmentQueryIndexes);
        for (final int queryIndex : seqvarQueryIndexes) {
            final PerQueryAlignmentData align = alignmentDataMap.get(queryIndex);
            final PerQueryAlignmentData seqvar = seqvarDataMap.get(queryIndex);

            final Map<String, String> alignSeqVarsMap = align.refPositionReadIndexToBaseMap;
            final Map<String, String> varSeqVarsMap = seqvar.refPositionReadIndexToBaseMap;

            assertEquals(String.format("queryIndex=%d alignSeqVarsMap.size()(%d) should equal varSeqVarsMap.size()(%d)",
                    queryIndex,
                    alignSeqVarsMap.size(), varSeqVarsMap.size()),
                    alignSeqVarsMap.size(), varSeqVarsMap.size());
            for (final Map.Entry<String, String> varEntry : varSeqVarsMap.entrySet()) {
                // Make sure the sequence variations match
                final String varEntryBases = varEntry.getValue();
                final String alignEntryBases = alignSeqVarsMap.get(varEntry.getKey());
                assertNotNull(String.format("queryIndex=%d Could not find alignSeqVarsMap entry for %s",
                        queryIndex, varEntry.getKey()),
                        alignEntryBases);
                assertEquals(String.format("queryIndex=%d alignEntryBases(%s) should equal varEntryBases(%s)",
                        queryIndex, alignEntryBases, varEntryBases),
                        alignEntryBases, varEntryBases);
            }
        }
    }

    private byte[] byteArray(final int... bytes) {
        final byte[] result = new byte[bytes.length];
        for (int i = 0; i < bytes.length; i++) {
            result[i] = (byte) bytes[i];
        }
        return result;
    }


    /**
     * Test SamToCompact convertBases() can handle requesting
     * to convert bases before and after what is provided by
     * the reference genome. And thet = is returned at the appropriate
     * times, otherwise the READ base is returned.
     */
    @Test
    public void samToCompactConvertBases() {
        final MutableString seq = new MutableString();
        //          01234567
        seq.append("CAGTGTAC");
        String[] refs = {seq.toString()};

        RandomAccessSequenceTestSupport genome = new RandomAccessSequenceTestSupport(refs);


        SAMToCompactMode samToCompact = new SAMToCompactMode();
        samToCompact.setGenome(genome);

        // AT the left side of the reference
        //           01234567890
        String readBases = "CTGTGTAC";
        String convertBases = samToCompact.convertBases(0, 0, readBases.getBytes(), 0, 3);
        assertEquals("=T=", convertBases);

        // AT the right side of the reference
        //           01234567890
        readBases = "GTGC";
        convertBases = samToCompact.convertBases(0, 5, readBases.getBytes(), 1, 1 + 3);
        assertEquals("=G=", convertBases);
        
        // Before the left side of the reference
        readBases = "AAACAG";
        convertBases = samToCompact.convertBases(0, -3, readBases.getBytes(), 0, 0 + 6);
        assertEquals("AAA===", convertBases);

        // Beyond the right side of the reference
        //           01234567890
        readBases = "ACGGTACGCAD";
        convertBases = samToCompact.convertBases(0, 4, readBases.getBytes(), 3, 3 + 7);
        assertEquals("====GCA", convertBases);


        // Beyond left AND right
        readBases = "GGCAGTGTACTT";
        convertBases = samToCompact.convertBases(0, -2, readBases.getBytes(), 0, 12);
        assertEquals("GG========TT", convertBases);
    }



}
