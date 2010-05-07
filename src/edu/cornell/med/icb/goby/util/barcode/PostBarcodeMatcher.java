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

package edu.cornell.med.icb.goby.util.barcode;

import it.unimi.dsi.lang.MutableString;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

/**
 * This is used to determine the barcode on files where the barcode should be found at the END of the sequence.
 *
 * @author Kevin Dorff
 */
public class PostBarcodeMatcher extends BarcodeMatcher{

    /** Logging. */
    private static final Log LOG = LogFactory.getLog(BarcodeMatcherResult.class);

    /**
     * Create a PostBarcodeMatcher.
     * @param barcodesStrArray MutableString version of the barcodes strings (including trailing adapter).
     * @param barcodeLength the length of the barcode (not including the adapter)
     * @param allowedMismatches the number of allowed mismatches when matching to be acceptable
     */
    public PostBarcodeMatcher(final String[] barcodesStrArray, final int barcodeLength, int allowedMismatches) {
        init(barcodesStrArray, barcodeLength, allowedMismatches);
    }

    /**
     * Searches for toFind in sequence. It is assumed toFind will be at the end of sequence
     * and that the beginning of toFind is the most important part.
     * @param sequence the sequence to search
     * @param toFind the string to find
     * @param barcodeIndex the index of the barcode within all of the barcodes
     * @param minMatchLength the minimum match length to consider
     * @return returns the number of differences between sequence and toFind
     */
    BarcodeMatcherResult bestMatch(final MutableString sequence, final MutableString toFind, final int barcodeIndex, final int minMatchLength) {
        final OverlapResult overlapResult = overlapPortion(sequence, toFind);
        int leastNumDiffs = Integer.MAX_VALUE;
        int leastNumDiffsBarcodeStartPos = 0;
        int pos = 0;
        for (int matchLength = overlapResult.length; matchLength >= minMatchLength; matchLength--, pos++) {
            final int numDiffs = numDifferences(sequence, toFind, overlapResult.start + pos, 0, matchLength);
            if (numDiffs <= leastNumDiffs) {
                leastNumDiffs = numDiffs;
                leastNumDiffsBarcodeStartPos = sequence.length() - matchLength;
                if (leastNumDiffs == 0) {
                    break;
                }
            }
        }
        final int preMatchSequenceLength = sequence.length();
        final int barcodeMatchedLength = preMatchSequenceLength - leastNumDiffsBarcodeStartPos;
        final int actualSequenceLength = preMatchSequenceLength - barcodeMatchedLength;

        return new BarcodeMatcherResult(
                barcodeIndex,
                leastNumDiffs,
                0, actualSequenceLength,
                leastNumDiffsBarcodeStartPos, barcodeMatchedLength);
    }

    /**
     * Determine the overlap portion of the two strings given their lengths
     * @param sequence the string we are searching
     * @param toFind the string we are looking for (at the end of search)
     * @return OverlapResult which specifies start and length
     */
    OverlapResult overlapPortion(final MutableString sequence, final MutableString toFind) {
        OverlapResult overlapResult = new OverlapResult();
        int sequenceLength = sequence.length();
        int toFindLength = toFind.length();
        if (sequenceLength >= toFindLength) {
            overlapResult.start = (sequenceLength - toFindLength);
            overlapResult.length = toFindLength;
        } else if (sequenceLength < toFindLength) {
            overlapResult.start = 0;
            overlapResult.length = sequenceLength;
        }
        return overlapResult;
    }

}
