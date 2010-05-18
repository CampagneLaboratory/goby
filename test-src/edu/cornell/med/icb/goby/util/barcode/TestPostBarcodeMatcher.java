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

import edu.cornell.med.icb.goby.readers.FastXEntry;
import edu.cornell.med.icb.goby.readers.FastXReader;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.IOUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import static org.junit.Assert.assertEquals;
import org.junit.Test;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Random;

/**
 * Describe class here.
 *
 * @author Kevin Dorff
 */
public class TestPostBarcodeMatcher {
    /** Logging. */
    private static final Log LOG = LogFactory.getLog(TestPostBarcodeMatcher.class);

    private static final Random rnd = new Random(System.currentTimeMillis());

    public static final char[] ACGT = new char[] {'A', 'C', 'G', 'T'};

    public static String SPACES = null;

    public static String[] TEST_BARCODES = new String[] {
            "TCACTTCGTATGCCGTCTTCTGCTTG",
            "TCATCTCGTATGCCGTCTTCTGCTTG",
            "TCCACTCGTATGCCGTCTTCTGCTTG",
            "TCCGTTCGTATGCCGTCTTCTGCTTG",
            "TCCTATCGTATGCCGTCTTCTGCTTG",
            "TCGATTCGTATGCCGTCTTCTGCTTG",
            "TCGCGTCGTATGCCGTCTTCTGCTTG",
            "TCTAGTCGTATGCCGTCTTCTGCTTG",
            "TCTCCTCGTATGCCGTCTTCTGCTTG",
            "TCTGATCGTATGCCGTCTTCTGCTTG",
            "TTAAGTCGTATGCCGTCTTCTGCTTG",
            "TAACGTCGTATGCCGTCTTCTGCTTG"
    };

    @Test
    public void testNumDifferences() {
        final PostBarcodeMatcher matcher = new PostBarcodeMatcher(null, 0, 100);
        assertEquals(3, matcher.numDifferences(new MutableString("ABC"), new MutableString("123"), 0, 0, 3));
        assertEquals(0, matcher.numDifferences(new MutableString("ABC"), new MutableString("ABC"), 0, 0, 3));
        assertEquals(1, matcher.numDifferences(new MutableString("ABC"), new MutableString("BBC"), 0, 0, 3));
        assertEquals(1, matcher.numDifferences(new MutableString("ABC"), new MutableString("ACC"), 0, 0, 3));
        assertEquals(1, matcher.numDifferences(new MutableString("ABC"), new MutableString("ABB"), 0, 0, 3));
        assertEquals(2, matcher.numDifferences(new MutableString("ABC"), new MutableString("ACB"), 0, 0, 3));


        assertEquals(0, matcher.numDifferences(new MutableString("ABCDEF"), new MutableString("DEFABC"), 3, 0, 3));
        assertEquals(0, matcher.numDifferences(new MutableString("DEFABC"), new MutableString("ACBDEF"), 0, 3, 3));

        assertEquals(1, matcher.numDifferences(new MutableString("ABBDEF"), new MutableString("DEFABC"), 0, 3, 3));
        assertEquals(1, matcher.numDifferences(new MutableString("DEFACC"), new MutableString("ACBDEF"), 3, 0, 3));

        assertEquals(0, matcher.numDifferences(new MutableString("ABCDEF"), new MutableString("XBCDEX"), 1, 1, 3));

    }

    @Test
    public void testNumDifferencesStartLength() {
        final PostBarcodeMatcher matcher = new PostBarcodeMatcher(null, 0, 100);
        assertEquals(0, matcher.numDifferences(new MutableString("DEFABC"), new MutableString("ABC"), 3, 0, 3));
        assertEquals(1, matcher.numDifferences(new MutableString("ABC"), new MutableString("ACCDEF"), 0, 0, 3));
    }

    @Test
    public void testNumDifferences2() {
        final PostBarcodeMatcher matcher = new PostBarcodeMatcher(null, 0, 100);
        // Perfect match when considering 5
        assertEquals(0, matcher.numDifferences(m("ABCDEFGHIJKL"), m("HIJKLMNOPQRS"), 7, 0, 5));
        // Almost perfect match, just the X is different
        assertEquals(1, matcher.numDifferences(m("ABCDEFGHIJKL"), m("HIJKXMNOPQRS"), 7, 0, 5));
        // Almost perfect match, just the X is different
        assertEquals(1, matcher.numDifferences(m("ABCDEFGHIJKL"), m("HIJXLMNOPQRS"), 7, 0, 5));
        // X is different but it is outside of the 5-match we want
        assertEquals(0, matcher.numDifferences(m("ABCDEFGHIJKL"), m("HIJKLXNOPQRS"), 7, 0, 5));
    }

    @Test
    public void testBestMatch() {
        final PostBarcodeMatcher matcher = new PostBarcodeMatcher(null, 0, 100);
        assertEquals(new BarcodeMatcherResult(new HashMap<String, Object>() {{
                    put("barcodeIndex", 1);
                    put("numberOfDiffs", 3);
                    put("sequenceStartPosition", 0);
                    put("sequenceLength", 0);
                    put("barcodeStartPosition", 0);
                    put("barcodeMatchLength", 3);
                }}), matcher.bestMatch(m("ABC"), m("123"), 1, 3));

        assertEquals(new BarcodeMatcherResult(new HashMap<String, Object>() {{
                    put("barcodeIndex", 2);
                    put("numberOfDiffs", 0);
                    put("sequenceStartPosition", 0);
                    put("sequenceLength", 2);
                    put("barcodeStartPosition", 2);
                    put("barcodeMatchLength", 4);
                }}), matcher.bestMatch(m("ABCDEF"), m("CDEFGHIJKL"), 2, 3));

        assertEquals(new BarcodeMatcherResult(new HashMap<String, Object>() {{
                    put("barcodeIndex", 3);
                    put("numberOfDiffs", 1);
                    put("sequenceStartPosition", 0);
                    put("sequenceLength", 3);
                    put("barcodeStartPosition", 3);
                    put("barcodeMatchLength", 3);
                }}), matcher.bestMatch(m("ABCDEF"), m("DGFGHIJKL"), 3, 3));
    }

    @Test
    public void testOverlapPortion() {
        final PostBarcodeMatcher matcher = new PostBarcodeMatcher(null, 0, 100);
        assertEquals(new OverlapResult(0, 3), matcher.overlapPortion(m("ABC"), m("123")));
        assertEquals(new OverlapResult(0, 6), matcher.overlapPortion(m("ABCDEF"), m("CDEFGHIJKL")));
        assertEquals(new OverlapResult(4, 6), matcher.overlapPortion(m("ABCDEFEFGH"), m("DGFGHI")));
    }

    @Test
    public void testMatchBarcodes() throws IOException {
        final long time = System.currentTimeMillis();
        final PostBarcodeMatcher matcher = new PostBarcodeMatcher(TEST_BARCODES, 5, 2);

        int i = 0;
        int numMatches = 0;
        int numAmbiguous = 0;
        int numNoMatches = 0;
        for (final FastXEntry entry :
                new FastXReader("test-data/barcode/post-test.fa")) {
            final BarcodeMatcherResult bestMatch = matcher.matchSequence(entry.getSequence());
            if (bestMatch != null) {
                final int[] result = parseTestHeader(entry.getEntryHeader());
                final int sequenceLength = result[0];
                final int barcodeLen = result[1];
                final int whichBarcodeIndex = result[2];
                assertEquals(sequenceLength, bestMatch.getSequenceLength());
                assertEquals(barcodeLen, bestMatch.getBarcodeMatchLength());
                assertEquals(whichBarcodeIndex, bestMatch.getBarcodeIndex());
                if (LOG.isDebugEnabled()) {
                    LOG.info(String.format("Entry %d", i));
                    LOG.info(String.format("    %s %s",
                            bestMatch.sequenceOf(entry.getSequence()),
                            bestMatch.sequenceBarcodeOf(entry.getSequence())));
                    LOG.info(String.format("    %s %s (%d)",
                            spaces(bestMatch.getBarcodeStartPosition()),
                            bestMatch.barcodeWithMatchingAdapterOf(matcher),
                            bestMatch.getNumberOfDiffs()));
                    LOG.info(String.format("    %s %s",
                            spaces(bestMatch.getBarcodeStartPosition()),
                            bestMatch.barcodeOnlyOf(matcher)));
                }
                numMatches++;
                if (bestMatch.isAmbiguous()) {
                    numAmbiguous++;
                }
            } else {
                if (LOG.isDebugEnabled()) {
                    LOG.debug(String.format("No match for entry %d", i));
                }
                numNoMatches++;
            }
            i++;
            if (i % 100000 == 0) {
                LOG.info(i);
            }
        }
        LOG.info(String.format("Num matches = %d, Num Ambiguous = %d, Num no matches = %d",
                numMatches, numAmbiguous, numNoMatches));
        final long duration = (System.currentTimeMillis() - time) / 1000;
        LOG.info(String.format("Time to parse 8 million reads %d seconds", duration));
    }

    public static int[] parseTestHeader(final MutableString testHeader) {
        final String[] parts = testHeader.toString().split(":");
        final int[] result = new int[parts.length];
        for (int i = 0; i < parts.length; i++) {
            result[i] = Integer.parseInt(parts[i]);
        }
        return result;
    }

    public void createPostBarcodedFile() throws IOException {
        PrintWriter out = null;
        try {
            out = new PrintWriter("test-data/barcode/post-test.fa");
            out.println("# >sequenceLength:barcodeLen:whichBarcodeIndex");
            final int minBarcodeLen = 5;
            final int maxBarcodeLen = TEST_BARCODES[0].length();
            final int randomSpan = maxBarcodeLen - minBarcodeLen + 1;
            final int totalSequenceLength = maxBarcodeLen + 5;
            for (int i = 0; i < 200; i++) {
                final int barcodeLen = rnd.nextInt(randomSpan) + minBarcodeLen;
                final int sequenceLength = totalSequenceLength - barcodeLen;
                final int whichBarcodeIndex = rnd.nextInt(TEST_BARCODES.length);
                out.printf(">%d:%d:%d%n", sequenceLength, barcodeLen, whichBarcodeIndex);
                for (int j = 0; j < sequenceLength; j++) {
                    out.print(ACGT[rnd.nextInt(ACGT.length)]);
                }
                out.print(TEST_BARCODES[whichBarcodeIndex].substring(0, barcodeLen));
                out.println();
            }
        } finally {
            IOUtils.closeQuietly(out);
        }
    }

    private static MutableString perturbBarcode(final int numDiffs, final MutableString barcode) {
        if (numDiffs == 0) {
            return barcode;
        }
        final IntList perturbLocations = new IntArrayList(numDiffs);
        for (int i = 0; i < numDiffs; i++) {
            while (true) {
                final int position = rnd.nextInt(barcode.length());
                if (!perturbLocations.contains(position)) {
                    perturbLocations.add(position);
                    break;
                }
            }
        }
        for (final int perturbLocation : perturbLocations) {
            while (true) {
                final char newChar = TestPostBarcodeMatcher.ACGT[rnd.nextInt(TestPostBarcodeMatcher.ACGT.length)];
                if (barcode.charAt(perturbLocation) != newChar) {
                    barcode.setCharAt(perturbLocation, newChar);
                    break;
                }
            }
        }
        return barcode;
    }

    public static synchronized String spaces(final int num) {
        if (SPACES == null) {
            final StringBuffer sb = new StringBuffer();
            for (int i = 0; i < 2048; i++) {
                sb.append(' ');
            }
            SPACES = sb.toString();
        }
        return SPACES.substring(0, num);
    }

    public MutableString m(final String str) {
        return new MutableString(str);
    }
}
