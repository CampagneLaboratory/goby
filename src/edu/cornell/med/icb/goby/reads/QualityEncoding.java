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

/**
 * Quality score encodings supported by Goby. See
 * <a href="http://en.wikipedia.org/wiki/FASTQ_format">http://en.wikipedia.org/wiki/FASTQ_format</a>
 * for a description of these various encodings.
 *
 * Sanger format can encode a Phred quality score from 0 to 93 using ASCII 33 to 126
 * (although in raw read data the Phred quality score rarely exceeds 60, higher scores are
 * possible in assemblies or read maps).  Illumina 1.3+ format can encode a Phred quality
 * score from 0 to 62 using ASCII 64 to 126 (although in raw read data Phred scores from
 * 0 to 40 only are expected). Solexa/Illumina 1.0 format can encode a Solexa/Illumina
 * quality score from -5 to 62 using ASCII 59 to 126 (although in raw read data Solexa
 * scores from -5 to 40 only are expected)
 *
 * <pre>
 * SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
 ^ ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
 * ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 * !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
 * |                         |    |        |                              |                     |
 * 33                        59   64       73                            104                   126
 *
 * S - Sanger       Phred+33,  41 values  (0, 40)
 * I - Illumina 1.3 Phred+64,  41 values  (0, 40)
 * X - Solexa       Solexa+64, 68 values (-5, 62)
 * </pre>
 */
public enum QualityEncoding {
    /**
     * Sanger encoding.
     */
    SANGER(33),
    /**
     * Illumina encoding.
     */
    ILLUMINA(64),
    /**
     * Solexa encoding.
     */
    SOLEXA(64);

    /**
     * The offset used to convert quality scores to/from ASCII characters.
     */
    private final int phredOffset;

    /**
     * Create a new QualityEncoding with the specified offset.
     * @param phredOffset The offset used to convert quality scores to/from ASCII characters.
     */
    private QualityEncoding(final int phredOffset) {
        this.phredOffset = phredOffset;
    }

    /**
     * Convert a quality score value to the ascii encoded equivalent.
     * @param qualityScore The quality score to convert.
     * @return The ASCII character representation of the score
     */
    public char qualityScoreToAsciiEncoding(final byte qualityScore) {
        return (char) (qualityScore + phredOffset);
    }

    /**
     * Convert a quality score value from the ascii encoded equivalent.
     * @param asciiCharacter the character to convert
     * @return The score value converted from the ascii encoding
     */
    public byte asciiEncodingToQualityScore(final char asciiCharacter) {
        return (byte) (asciiCharacter - phredOffset);
    }
}
