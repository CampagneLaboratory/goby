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

package edu.cornell.med.icb.goby.reads;

import static junit.framework.Assert.assertEquals;
import org.apache.commons.lang.math.IntRange;
import org.apache.commons.lang.math.Range;
import org.junit.Test;

/**
 * Validates the functionality of {@link edu.cornell.med.icb.goby.reads.QualityEncoding}.
 * Sanger format can encode a Phred quality score from 0 to 93 using ASCII 33 to 126
 * (although in raw read data the Phred quality score rarely exceeds 60, higher scores are
 * possible in assemblies or read maps).  Illumina 1.3+ format can encode a Phred quality
 * score from 0 to 62 using ASCII 64 to 126 (although in raw read data Phred scores from
 * 0 to 40 only are expected). Solexa/Illumina 1.0 format can encode a Solexa/Illumina
 * quality score from -5 to 62 using ASCII 59 to 126 (although in raw read data Solexa
 * scores from -5 to 40 only are expected)
 * <p/>
 * <pre>
 * SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
 * ^ ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
 * ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 * !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
 * |                         |    |        |                              |                     |
 * 33                        59   64       73                            104                   126
 * <p/>
 * S - Sanger       Phred+33,  41 values  (0, 40)
 * I - Illumina 1.3 Phred+64,  41 values  (0, 40)
 * X - Solexa       Solexa+64, 68 values (-5, 62)
 * </pre>
 */
public class TestQualityEncoding {
    /**
     * Expected character set for quality encoded values.
     */
    private final String encodedQualityCharacters =
            "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                    + "[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

    /**
     * Validate score conversion for {@link edu.cornell.med.icb.goby.reads.QualityEncoding#SANGER}.
     */
    @Test
    public void sanger() {
        final Range qualityRange = new IntRange(0, 40);
        final String qualityChars = encodedQualityCharacters.substring(0, 40);
        final QualityEncoding encoding = QualityEncoding.SANGER;
        for (int i = 0; i < 40; i++) {
            final byte score = (byte) (qualityRange.getMinimumInteger() + i);
            assertEquals("wrong character for score " + score,
                    qualityChars.charAt(i), encoding.phredQualityScoreToAsciiEncoding(score));
            assertEquals("wrong score for character " + qualityChars.charAt(i),
                    score, encoding.asciiEncodingToPhredQualityScore(qualityChars.charAt(i)));
        }

        final QualityEncoding solexaEncoding = QualityEncoding.ILLUMINA;
        for (byte qPhred = 0; qPhred < 93; qPhred++) {
            final char sanger = solexaEncoding.phredQualityScoreToAsciiEncoding(qPhred);
            final byte roundTrip = solexaEncoding.asciiEncodingToPhredQualityScore(sanger);
            assertEquals(qPhred, roundTrip);

        }
    }

    /**
     * Validate score conversion for
     * {@link edu.cornell.med.icb.goby.reads.QualityEncoding#ILLUMINA}.
     */
    @Test
    public void illumina() {
        final Range qualityRange = new IntRange(0, 40);
        final String qualityChars = encodedQualityCharacters.substring(31, 71);
        final QualityEncoding encoding = QualityEncoding.ILLUMINA;
        for (int i = 0; i < 40; i++) {
            final byte score = (byte) (qualityRange.getMinimumInteger() + i);
            assertEquals("wrong character for score " + score,
                    qualityChars.charAt(i), encoding.phredQualityScoreToAsciiEncoding(score));
            assertEquals("wrong score for character " + qualityChars.charAt(i),
                    score, encoding.asciiEncodingToPhredQualityScore(qualityChars.charAt(i)));
        }

        final QualityEncoding solexaEncoding = QualityEncoding.ILLUMINA;
        for (byte qPhred = 0; qPhred < 62; qPhred++) {
            final char illumina13 = solexaEncoding.phredQualityScoreToAsciiEncoding(qPhred);
            final byte roundTrip = solexaEncoding.asciiEncodingToPhredQualityScore(illumina13);
            assertEquals(qPhred, roundTrip);

        }
    }

    /**
     * Validate score conversion for {@link edu.cornell.med.icb.goby.reads.QualityEncoding#SOLEXA}.
     */
    @Test
    public void solexa() {
        final QualityEncoding solexaEncoding = QualityEncoding.SOLEXA;


        for (byte qPhred = 1; qPhred < 62; qPhred++) {
            System.out.println("qPhred="+qPhred);
            final char solexa = solexaEncoding.phredQualityScoreToAsciiEncoding(qPhred);

            final byte roundTripSolexa = solexaEncoding.asciiEncodingToPhredQualityScore(solexa);

            assertEquals(qPhred, roundTripSolexa);

        }
    }
}
