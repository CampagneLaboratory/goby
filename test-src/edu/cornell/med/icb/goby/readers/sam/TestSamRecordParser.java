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

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import org.junit.Test;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.List;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

/**
 * Test SamRecordParser.
 */
public class TestSamRecordParser {


    //
    //  testSamToCompactTrickCase1-3 fails because this the sam reference builder requires an MD:Z tag.
    //

    /*
    @Test
    public void testSamToCompactTrickCase3() throws IOException {
        final String inputFile = "test-data/splicedsamhelper/tricky-spliced-3.sam";
        final SAMFileReader parser = new SAMFileReader(new FileInputStream(inputFile));
        parser.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        final SamRecordParser recordParser = new SamRecordParser();
        for (final SAMRecord samRecord : new SAMRecordIterable(parser.iterator())) {

            final List<GobySamRecord> segments = recordParser.processRead(samRecord);
            if (segments.get(0).readNum == 0) {
                assertEquals("Incorrect number of segments", 1, segments.size());
                GobySamRecord first = segments.get(0);
                assertEquals(170769 - 1, first.getPosition());
                assertTrue(first.hasMate);
                assertEquals(216048 - 1, first.mateStartPosition);
            } else {
                assertEquals("Incorrect number of segments", 2, segments.size());
                GobySamRecord second = segments.get(0);
                assertEquals(216048 - 1, second.getPosition());
                assertTrue(second.hasMate);
                assertEquals(170769 - 1, second.mateStartPosition);
            }
        }
    }
    */

    @Test
    // like 9 no genome
    public void testSamToCompactTrickCase9NoGenome() throws IOException {
        final String inputFile = "test-data/splicedsamhelper/tricky-spliced-9.sam";
        final SAMFileReader parser = new SAMFileReader(new FileInputStream(inputFile));
        parser.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        final SamRecordParser recordParser = new SamRecordParser();
        for (final SAMRecord samRecord : new SAMRecordIterable(parser.iterator())) {

            final List<GobySamRecord> segments = recordParser.processRead(samRecord);

            assertEquals("Incorrect number of segments", 2, segments.size());
            GobySamRecord first = segments.get(0);

            assertEquals(0, first.getFragmentIndex());
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

            final List<GobySamRecord> segments = recordParser.processRead(samRecord);

            assertEquals("Incorrect number of segments", 2, segments.size());
            GobySamRecord first = segments.get(0);

            if (first.readNum == 0) {
                assertEquals(0, first.getFragmentIndex());
                assertEquals(3 - 1, first.getPosition());
            } else if (first.readNum == 1) {
                assertEquals(0, first.getFragmentIndex());
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

            final List<GobySamRecord> segments = recordParser.processRead(samRecord);

            assertEquals("Incorrect number of segments", 2, segments.size());
            GobySamRecord first = segments.get(0);
            assertEquals(0, first.getFragmentIndex());
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

            final List<GobySamRecord> segments = recordParser.processRead(samRecord);

            assertEquals("Incorrect number of segments", 2, segments.size());
            GobySamRecord first = segments.get(0);
            GobySamRecord second = segments.get(1);

            assertEquals(15013, first.getPosition());
            assertEquals(3, first.getQueryPosition());
            assertEquals(0, first.getFragmentIndex());
            assertEquals(25, first.getQueryAlignedLength());

            assertEquals(15795, second.getPosition());
            assertEquals(28, second.getQueryPosition());
            assertEquals(1, second.getFragmentIndex());
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

            final List<GobySamRecord> segments = recordParser.processRead(samRecord);

            assertEquals("Incorrect number of segments", 2, segments.size());
            GobySamRecord first = segments.get(0);
            GobySamRecord second = segments.get(1);

            assertEquals(15013, first.getPosition());
            assertEquals(3, first.getQueryPosition());
            assertEquals("AAT", first.getSoftClippedBasesLeft());
            assertEquals("", first.getSoftClippedBasesRight());
            assertEquals(25, first.getQueryAlignedLength());

            assertEquals(15795, second.getPosition());
            assertEquals(28, second.getQueryPosition());
            assertEquals(1, second.getFragmentIndex());
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
            final List<GobySamRecord> segments = recordParser.processRead(samRecord);

            assertEquals("Incorrect number of segments", 2, segments.size());
            GobySamRecord first = segments.get(0);
            GobySamRecord second = segments.get(1);

            assertEquals(4, first.getPosition());
            assertEquals(3, first.getQueryPosition());
            assertEquals("AAT", first.getSoftClippedBasesLeft());  // Prev test A=T
            assertEquals("", first.getSoftClippedBasesRight());
            assertEquals(25, first.getQueryAlignedLength());

            assertEquals(40 - 1, second.getPosition());
            assertEquals(28, second.getQueryPosition());
            assertEquals(1, second.getFragmentIndex());
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
            final List<GobySamRecord> segments = recordParser.processRead(samRecord);

            assertEquals("Incorrect number of segments", 3, segments.size());

            final GobySamRecord first = segments.get(0);
            final GobySamRecord second = segments.get(1);
            final GobySamRecord third = segments.get(2);

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
            final List<GobySamRecord> segments = recordParser.processRead(samRecord);

            assertEquals("Incorrect number of segments", 1, segments.size());
            final GobySamRecord segment = segments.get(0);
            if (segment.readNum == 0) {
                final GobySamRecord first = segment;
                assertEquals(190077 - 1, first.getPosition());
                assertEquals(13, first.getQueryPosition());
                assertEquals("AGTGGCAGCACGA", first.getSoftClippedBasesLeft());
                assertEquals("TGCT", first.getSoftClippedBasesRight());
                assertEquals(51, first.getQueryAlignedLength());
                assertEquals(50, first.getTargetAlignedLength());

                assertEquals(4, first.getSequenceVariationsCount());

                GobyQuickSeqvar seqvar = first.getSequenceVariations(0);
                assertEquals("T", seqvar.getFrom());
                assertEquals("C", seqvar.getTo());
                assertEquals(22, seqvar.getReadIndex());
                assertEquals(9, seqvar.getPosition());
                assertArrayEquals(byteArray(3), seqvar.getToQuality().toByteArray());  // 8 returned in the old test, 3 is right

                seqvar = first.getSequenceVariations(1);
                assertEquals("T", seqvar.getFrom());
                assertEquals("G", seqvar.getTo());    // Original test said A
                assertEquals(26, seqvar.getReadIndex());
                assertEquals(13, seqvar.getPosition());
                assertArrayEquals(byteArray(34), seqvar.getToQuality().toByteArray());  // 24 in the old test, this is right

                seqvar = first.getSequenceVariations(2);
                assertEquals("C", seqvar.getFrom());
                assertEquals("A", seqvar.getTo());
                assertEquals(33, seqvar.getReadIndex());
                assertEquals(20, seqvar.getPosition());
                assertArrayEquals(byteArray(28), seqvar.getToQuality().toByteArray());  // 14 in the old test

                seqvar = first.getSequenceVariations(3);
                assertEquals("-T", seqvar.getFrom());
                assertEquals("CA", seqvar.getTo());   // Original test said CG
                assertEquals(35, seqvar.getReadIndex());
                assertEquals(21, seqvar.getPosition());
                assertArrayEquals(byteArray(27, 35), seqvar.getToQuality().toByteArray());  // 3, 15 in the old test

            } else if (segment.readNum == 1) {
                final GobySamRecord second = segment;
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
                assertArrayEquals(byteArray(35), seqvar.getToQuality().toByteArray());
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
}
