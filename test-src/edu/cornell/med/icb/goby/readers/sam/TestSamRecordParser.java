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

import edu.cornell.med.icb.goby.alignments.*;
import edu.cornell.med.icb.goby.modes.SAMToCompactMode;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceTestSupport;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.lang.MutableString;
import junit.framework.Assert;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.Map;

import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.*;

/**
 * Test SamRecordParser.
 */
public class TestSamRecordParser {

    private static final Logger LOG = Logger.getLogger(TestSamRecordParser.class);

    private static final String BASE_TEST_INPUT_DIR = "test-data/splicedsamhelper";
    private static final String BASE_TEST_OUTPUT_DIR = "test-results/splicedsamhelper";

    //
    //  testSamToCompactTrickCase1-3 fails because this the sam reference builder requires an MD:Z tag.
    //

    @BeforeClass
    public static void beforeClass() {
        new File(BASE_TEST_OUTPUT_DIR).mkdirs();
    }

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
            final GobySamSegment first = gobySamRecord.getSegment(0);

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
            final GobySamSegment first = gobySamRecord.getSegment(0);

            if (gobySamRecord.readNum == 0 || gobySamRecord.readNum == 1) {
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
            final GobySamSegment first = gobySamRecord.getSegment(0);
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
            final GobySamSegment first = gobySamRecord.getSegment(0);
            final GobySamSegment second = gobySamRecord.getSegment(1);

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
            final GobySamSegment first = gobySamRecord.getSegment(0);
            final GobySamSegment second = gobySamRecord.getSegment(1);

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
            final GobySamSegment first = gobySamRecord.getSegment(0);
            final GobySamSegment second = gobySamRecord.getSegment(1);

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
                assertEquals(190077 - 1, segment.getPosition());
                assertEquals(13, segment.getQueryPosition());
                assertEquals("AGTGGCAGCACGA", segment.getSoftClippedBasesLeft());
                assertEquals("TGCT", segment.getSoftClippedBasesRight());
                assertEquals(51, segment.getQueryAlignedLength());
                assertEquals(50, segment.getTargetAlignedLength());

                assertEquals(5, segment.getSequenceVariationsCount());

                GobyQuickSeqvar seqvar = segment.getSequenceVariations(0);
                assertEquals("T", seqvar.getFrom());
                assertEquals("C", seqvar.getTo());
                assertEquals(22, seqvar.getReadIndex());
                assertEquals(9, seqvar.getPosition());
                assertArrayEquals(RoundTripAlignment.byteArray(3), seqvar.getToQualitiesAsBytes());  // 8 returned in the old test, 3 is right

                seqvar = segment.getSequenceVariations(1);
                assertEquals("T", seqvar.getFrom());
                assertEquals("G", seqvar.getTo());    // Original test said A
                assertEquals(26, seqvar.getReadIndex());
                assertEquals(13, seqvar.getPosition());
                assertArrayEquals(RoundTripAlignment.byteArray(34), seqvar.getToQualitiesAsBytes());  // 24 in the old test, this is right

                seqvar = segment.getSequenceVariations(2);
                assertEquals("C", seqvar.getFrom());
                assertEquals("A", seqvar.getTo());
                assertEquals(33, seqvar.getReadIndex());
                assertEquals(20, seqvar.getPosition());
                assertArrayEquals(RoundTripAlignment.byteArray(28), seqvar.getToQualitiesAsBytes());  // 14 in the old test

                seqvar = segment.getSequenceVariations(3);
                assertEquals("-", seqvar.getFrom());
                assertEquals("C", seqvar.getTo());   // Original test said CG
                assertEquals(35, seqvar.getReadIndex());
                assertEquals(21, seqvar.getPosition());
                assertArrayEquals(RoundTripAlignment.byteArray(27), seqvar.getToQualitiesAsBytes());  // 3, 15 in the old test

                seqvar = segment.getSequenceVariations(4);
                assertEquals("T", seqvar.getFrom());
                assertEquals("A", seqvar.getTo());   // Original test said CG
                assertEquals(36, seqvar.getReadIndex());
                assertEquals(22, seqvar.getPosition());
                assertArrayEquals(RoundTripAlignment.byteArray(35), seqvar.getToQualitiesAsBytes());  // 3, 15 in the old test

            } else if (gobySamRecord.readNum == 1) {
                //second's CIGAR is 20S48M
                assertEquals(190246 - 1, segment.getPosition());
                assertEquals(20, segment.getQueryPosition());
                assertEquals("CAGTGTCGTGGCTGCACGCC", segment.getSoftClippedBasesLeft());
                assertEquals("", segment.getSoftClippedBasesRight());
                assertEquals(48, segment.getQueryAlignedLength());
                assertEquals(48, segment.getTargetAlignedLength());

                assertEquals(1, segment.getSequenceVariationsCount());

                final GobyQuickSeqvar seqvar = segment.getSequenceVariations(0);
                assertEquals("A", seqvar.getFrom());
                assertEquals("C", seqvar.getTo());
                assertEquals(18, seqvar.getReadIndex());
                assertEquals(31, seqvar.getPosition());
                assertArrayEquals(RoundTripAlignment.byteArray(35), seqvar.getToQualitiesAsBytes());
            }
        }
    }

    /**
     * This comes from
     * Quality score difference.
     * See XAAOBVT  [Open Session]
     * chr1:45,881,903-45,881,942
     * T variation with qual score 2 should be 36
     * Sample = XAAOBVT
     * Read group = 1
     * ----------------------
     * Read name = 509.6.68.19057.157284
     *
     * @throws IOException error reading
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

            final GobyQuickSeqvar last = first.getSequenceVariations(12);
            assertEquals("T", last.getTo());
            assertEquals("C", last.getFrom());
            assertArrayEquals(RoundTripAlignment.byteArray(36), last.getToQualitiesAsBytes());
        }
    }

    /**
     * Same test as above, but write the sam to compact and then read the sequence variation form the compact.
     *
     * @throws IOException error reading
     */
    @Test
    public void testSamToCompactTrickCase17ViaWriter() throws IOException {
        final SAMToCompactMode importer = new SAMToCompactMode();
        importer.setInputFile("test-data/splicedsamhelper/tricky-spliced-17.sam");
        final String outputFilename = FilenameUtils.concat(BASE_TEST_OUTPUT_DIR, "tricky-spliced-17");
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

        final Alignments.SequenceVariation last = first.getSequenceVariations(12);
        assertEquals("T", last.getTo());
        assertEquals("C", last.getFrom());
        assertArrayEquals(RoundTripAlignment.byteArray(36), last.getToQuality().toByteArray());
    }

    /**
     * Test round trip of tricky-spliced-18.sam that input and output
     * sam compare.
     *
     * @throws IOException error
     */
    @Test
    public void testRoundTripTrickySpliced18() throws IOException {
        final RoundTripAlignment rtc = new RoundTripAlignment();
        rtc.inputGenomeFilename = TestGobyPaperTop5000s.findThousandGenome();
        rtc.sourceBamFilename = "test-data/splicedsamhelper/tricky-spliced-18.sam";
        rtc.destGobyBasename = FilenameUtils.concat(BASE_TEST_OUTPUT_DIR, "tricky-spliced-18");
        rtc.destBamFilename = FilenameUtils.concat(BASE_TEST_OUTPUT_DIR, "tricky-spliced-18.sam");
        rtc.testRoundTripAny();
        rtc.keepQualityScores = false;
        rtc.testRoundTripAny();
        rtc.keepSoftClips = false;
        rtc.testRoundTripAny();
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
                assertEquals(31 - 1, segment.getPosition());
                assertEquals(0, segment.getQueryPosition());
                assertEquals("", segment.getSoftClippedBasesLeft());
                assertEquals("", segment.getSoftClippedBasesRight());
                assertEquals(47, segment.getQueryAlignedLength());
                assertEquals(50, segment.getTargetAlignedLength());

                assertEquals(1, segment.getSequenceVariationsCount());

                final GobyQuickSeqvar seqvar = segment.getSequenceVariations(0);
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
                assertEquals(6 - 1, segment.getPosition());
                assertEquals(5, segment.getQueryPosition());
                assertEquals("AAAAA", segment.getSoftClippedBasesLeft());
                assertEquals("", segment.getSoftClippedBasesRight());
                assertEquals(45, segment.getQueryAlignedLength());
                assertEquals(45, segment.getTargetAlignedLength());

                assertEquals(1, segment.getSequenceVariationsCount());

                final GobyQuickSeqvar seqvar = segment.getSequenceVariations(0);
                assertEquals("A", seqvar.getFrom());
                assertEquals("G", seqvar.getTo());
                assertEquals(20, seqvar.getReadIndex());
                assertEquals(15, seqvar.getPosition());
                assertEquals(1, seqvar.getToQualitiesAsBytes().length);
                assertArrayEquals(RoundTripAlignment.byteArray(20), seqvar.getToQualitiesAsBytes());
            }
        }
    }

    /**
     * Compare the gsnap->sam created sequence variations we find when we parse the sam file with
     * SamRecordParser exactly match those generated with gsnap->compactAlignment->displaySequenceVariations(per-base).
     *
     * @throws IOException error
     */
    @Test
    public void testSeqVarReads() throws IOException {
        final String inputFile = "test-data/seq-var-test/seq-var-reads-gsnap.sam";
        final SAMFileReader parser = new SAMFileReader(new FileInputStream(inputFile));
        parser.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        final SamRecordParser recordParser = new SamRecordParser();
        final Int2ObjectMap<PerQueryAlignmentData> seqvarDataMap = TestIteratedSortedAlignment2.readSeqVarFile(
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
                final int readIndexIncrementValue = gobySamRecord.isReverseStrand() ? -1 : 1;
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

        verifySequenceVariationsMatch(seqvarQueryIndexes, samSeqvarQueryIndexes, seqvarDataMap, samSeqvarDataMap);
    }

    /**
     * Compare sequence varations.
     *
     * @param seqvarQueryIndexes    the sorted query indexes for the seqvarDataMap
     * @param alignmentQueryIndexes the sorted query indexes for the alignmentDataMap
     * @param seqvarDataMap         the query index to PerQueryAlignmentData for the "expected" values
     * @param alignmentDataMap      the query index to PerQueryAlignmentData for the "actual" values
     */
    public static void verifySequenceVariationsMatch(
            final int[] seqvarQueryIndexes,
            final int[] alignmentQueryIndexes,
            final Int2ObjectMap<PerQueryAlignmentData> seqvarDataMap,
            final Int2ObjectMap<PerQueryAlignmentData> alignmentDataMap) {
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
        final String[] refs = {seq.toString()};

        final RandomAccessSequenceTestSupport genome = new RandomAccessSequenceTestSupport(refs);


        final SAMToCompactMode samToCompact = new SAMToCompactMode();
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

    /**
     * We had a problem with the target indexes/lengths being written twice in SamToCompact. The second
     * time was wrong. This checks that it is right.
     *
     * @throws java.io.IOException error
     */
    @Test
    public void testTestTargetIndexCreation() throws IOException {
        final RoundTripAlignment rtc = new RoundTripAlignment();
        rtc.inputGenomeFilename = TestGobyPaperTop5000s.findMM9();
        rtc.sourceBamFilename = FilenameUtils.concat(BASE_TEST_INPUT_DIR, "JRODTYG-5000.sam.gz");
        rtc.destGobyBasename = FilenameUtils.concat(BASE_TEST_OUTPUT_DIR, "JRODTYG-5000");
        rtc.createCompactFromSam();

        final AlignmentReaderImpl gobyReader = new AlignmentReaderImpl(rtc.destGobyBasename);
        gobyReader.readHeader();
        final int[] targetsLengths = gobyReader.getTargetLength();
        final IndexedIdentifier identifiers = gobyReader.getTargetIdentifiers();
        final int chr1Index = identifiers.getInt(new MutableString("chr1"));
        Assert.assertEquals("Incorrect sequence length", 197195432, targetsLengths[chr1Index]);
    }

}
