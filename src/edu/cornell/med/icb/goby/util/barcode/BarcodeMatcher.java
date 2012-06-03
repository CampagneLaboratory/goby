/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
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

import it.unimi.dsi.fastutil.ints.Int2LongLinkedOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2LongMap;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

/**
 * Interface for classes to match Barcoded sequences.
 *
 * @author Kevin Dorff
 */
public abstract class BarcodeMatcher {

    /** Logging. */
    private static final Log LOG = LogFactory.getLog(BarcodeMatcher.class);

    /** Map of barcode index to the number of entries found for that barcode. */
    protected final Int2LongMap barcodeIndexToHitsMap = new Int2LongLinkedOpenHashMap();

    /** The barcode strings (including adapters). */
    protected MutableString[] barcodes;

    /** The length of the actual barcode (not including the adapter). */
    protected int barcodeLength;

    /** The number of allowed mismatchines when matching to a barcode. */
    protected int allowedMismatches;

    /**
     * Used by constructors to make a BarcodeMatcher.
     * @param barcodesStrArray MutableString version of the barcodes strings (including trailing adapter).
     * @param barcodeLength the length of the barcode (not including the adapter)
     * @param allowedMismatches the number of allowed mismatches when matching to be acceptable
     */
    public void init(final String[] barcodesStrArray, final int barcodeLength, final int allowedMismatches) {
        if (barcodesStrArray != null) {
            barcodes = new MutableString[barcodesStrArray.length];
            for (int i = 0; i < barcodesStrArray.length; i++) {
                barcodes[i] = new MutableString(barcodesStrArray[i]);
                barcodeIndexToHitsMap.put(i, 0);
            }
        } else {
            barcodes = null;
        }
        this.barcodeLength = barcodeLength;
        this.allowedMismatches = allowedMismatches;
    }

    BarcodeMatcherResult bestMatch(final MutableString sequence, final MutableString toFind, final int barcodeNum) {
        return bestMatch(sequence, toFind, barcodeNum, barcodeLength);
    }

    abstract BarcodeMatcherResult bestMatch(final MutableString sequence, final MutableString toFind, final int barcodeIndex, final int minMatchLength);

    /**
     * Determines which barcode this sequence matches or returns null if a barcode isn't found for
     * the sequence (withing the number of allowedMismatches).
     * @param sequence the sequence to get the barcode for
     * @return the barcode for the sequence or null if not found
     */
    public BarcodeMatcherResult matchSequence(final MutableString sequence) {
        BarcodeMatcherResult bestMatch = null;
        int numAtBestMatch = 0;
        for (int barcodeNum = 0; barcodeNum < barcodes.length; barcodeNum++) {
            final MutableString barcode = barcodes[barcodeNum];
            final BarcodeMatcherResult result = bestMatch(sequence, barcode, barcodeNum);
            if (bestMatch == null) {
                numAtBestMatch = 1;
                bestMatch = result;
            } else if (result.getNumberOfDiffs() < bestMatch.getNumberOfDiffs()) {
                numAtBestMatch = 1;
                bestMatch = result;
            } else if (result.getNumberOfDiffs() == bestMatch.getNumberOfDiffs()) {
                numAtBestMatch++;
            }
            if (bestMatch.getNumberOfDiffs() == 0) {
                // We found a perfect match. Just finish..
                break;
            }
        }
        if (bestMatch == null || bestMatch.getNumberOfDiffs() > allowedMismatches) {
            return null;
        }
        if (numAtBestMatch > 1) {
            bestMatch.setAmbiguous(true);
        }
        final int bestBarcodeNum = bestMatch.getBarcodeIndex();
        barcodeIndexToHitsMap.put(bestBarcodeNum, barcodeIndexToHitsMap.get(bestBarcodeNum) + 1);
        return bestMatch;
    }

    /**
     * Get the map of barcode index to the number of entries found for that barcode.
     * @return the map of barcode index to the number of entries found for that barcode
     */
    public Int2LongMap getBarcodeIndexToHitsMap() {
        return barcodeIndexToHitsMap;
    }

    /**
     * Get the barcode with adapter for the given barcode index.
     * @param index barcode index
     * @return the barcode with adapter
     */
    public MutableString getBarcodeWithAdapterAtIndex(final int index) {
        return barcodes[index];
    }

    /**
     * Get the barcode WITHOUT the adapter for the given barcode index.
     * @param index barcode index
     * @return the barcode without the adapter
     */
    public MutableString getBarcodeOnlyAtIndex(final int index) {
        return barcodes[index].substring(0, barcodeLength);
    }

    /**
     * The length of the actual barcode (without the adapter).
     * @return length of barcode
     */
    public int getBarcodeLength() {
        return barcodeLength;
    }

    /**
     * The number of allowed mismatches when matching a sequence to the barcodes.
     * @return the number of allowed mismatches
     */
    public int getAllowedMismatches() {
        return allowedMismatches;
    }

    /**
     * Determine the number of differences between two strings a and b. Start comparing
     * a from index aStart and b from index bStart. Compare for matchLength characters.
     * @param a the first string to compare
     * @param b the second string to compare
     * @param aStart the index of a to start the comparison
     * @param bStart the index of b to start the comparison
     * @param matchLength the number of characters to compare
     * @return the number of differences between a and b
     */
    int numDifferences(
            final MutableString a, final MutableString b, final int aStart, final int bStart, final int matchLength) {
        int numDifferences = 0;
        for (int i = 0; i < matchLength; i++) {
            if (a.charAt(aStart + i) != b.charAt(bStart + i)) {
                numDifferences++;
                if (numDifferences > allowedMismatches) {
                    break;
                }
            }
        }
        if (LOG.isDebugEnabled()) {
            LOG.debug(String.format("Comparing %s and %s (diff=%d)",
                    a.substring(aStart, aStart + matchLength),
                    b.substring(bStart, bStart + matchLength), numDifferences));
        }
        return numDifferences;
    }
}
