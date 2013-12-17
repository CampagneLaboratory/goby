/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 *  This file is part of the Goby IO API.
 *
 *     The Goby IO API is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     The Goby IO API is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with the Goby IO API.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.cornell.med.icb.goby.reads;

/**
 * Quality score encodings supported by Goby. See
 * <a href="http://en.wikipedia.org/wiki/FASTQ_format">http://en.wikipedia.org/wiki/FASTQ_format</a>
 * for a description of these various encodings. A much better description appeared in
 * Cock PJ et a, Nucleic Acids Research 2010, Vol 38, No 6.
 * <p/>
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
 *
 * S - Sanger       Phred+33,  41 values  (0, 40)
 * I - Illumina 1.3 Phred+64,  41 values  (0, 40)
 * X - Solexa       Solexa+64, 68 values (-5, 62)
 * </pre>
 */
public enum QualityEncoding {
    /**
     * BAM encoding.
     */
    BAM(0, false, 0, 93),
    /**
     * Sanger encoding.
     */
    SANGER(33, false, 0, 93),
    /**
     * Illumina encoding.
     */
    ILLUMINA(64, false, 0, 62),
    /**
     * Solexa encoding.
     */
    SOLEXA(59, true, -5, 62),

    /**
     * NO encoding.
     */
    PHRED(0, false, 0, 127);

    /**
     * The offset used to convert quality scores to/from ASCII characters.
     */
    private final int asciiOffset;
    /**
     * A flag to indicate that solexa scale conversion must occur.
     * See Cock PJ et a, Nucleic Acids Research 2010, Vol 38, No 6.
     */
    private final boolean solexaEncoding;
    private int minPhredScore;
    private int maxPhredScore;
    /**
     * When the force flag is set, force encoding values within range.
     */
    private boolean force = false;

    /**
     * Create a new QualityEncoding with the specified offset.
     *
     * @param asciiOffset    The offset used to convert quality scores to/from ASCII characters.
     * @param solexaEncoding Indicates whether solexa scale conversion must occur.
     * @param minPhredScore  the minimum allowed phred score
     * @param maxPhredScore  the maximum allowed phred score
     */
    private QualityEncoding(final int asciiOffset, final boolean solexaEncoding, final int minPhredScore, final int maxPhredScore) {
        this.asciiOffset = asciiOffset;
        this.solexaEncoding = solexaEncoding;
        this.minPhredScore = minPhredScore;
        this.maxPhredScore = maxPhredScore;
    }

    /**
     * Convert a quality score value to the ascii encoded equivalent.
     * Conversion implemented as described in Cock PJ et a, Nucleic Acids Research 2010, Vol 38, No 6.
     *
     * @param qPhred The quality score to convert in Phred scale.
     * @return The ASCII character representation of the score
     */
    public char phredQualityScoreToAsciiEncoding(final byte qPhred) {
        if (solexaEncoding) {
            final int qSolexa = (int) Math.round((10 * Math.log10(
                    Math.pow(10d,
                            ((double) qPhred) / 10d)
                            - 1
            )
            ));
          //  System.out.printf("qSolexa=%d%n", qSolexa);
            return (char) (qSolexa + asciiOffset);
        } else {
            return (char) (qPhred + asciiOffset);
        }


    }

    /**
     * Convert a quality score value from the ascii encoded equivalent to Phred scale.
     * Conversion implemented as described in Cock PJ et a, Nucleic Acids Research 2010, Vol 38, No 6.
     *
     * @param asciiCharacter the character to convert
     * @return The score value converted to Phred quality score scale.
     */
    public byte asciiEncodingToPhredQualityScore(final char asciiCharacter) {

        if (solexaEncoding) {
            final int qSolexa = (asciiCharacter - asciiOffset);
            // System.out.printf("qSolexa=%d%n", qSolexa);
            // convert Qsolexa to Qphred:
            byte qPhredScore = (byte) Math.round((10 * Math.log10(
                    Math.pow(10d, (((double) qSolexa) / 10d))
                            + 1)
            ));
          //  System.out.printf("qPhredScore=%d%n", qPhredScore);
            qPhredScore = forceInRange(qPhredScore);
            return qPhredScore;

        } else {
            return forceInRange((byte) (asciiCharacter - asciiOffset));
        }

    }

    private byte forceInRange(byte value) {
        if (!force) {
            return value;
        } else {

            value = (byte) Math.min(maxPhredScore, value);
            value = (byte) Math.max(minPhredScore, value);
            return value;
        }
    }

    /**
     * Returns true if the phredScore is within valid range for this encoding.
     *
     * @param phredScore quality score on Phred scale.
     * @return True  if the phredScore is within valid range for this encoding, false otherwise.
     */
    public boolean isWithinValidRange(final byte phredScore) {
        if (force) return true;
        return (phredScore <= maxPhredScore && phredScore >= minPhredScore);
    }

    public void setForce(boolean force) {
        this.force = force;
    }
}
