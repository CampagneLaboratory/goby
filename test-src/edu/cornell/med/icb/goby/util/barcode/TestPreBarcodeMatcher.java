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

import org.junit.Test;
import static org.junit.Assert.assertEquals;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.commons.io.IOUtils;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntArrayList;

import java.io.*;
import java.util.HashMap;
import java.util.Random;

import edu.cornell.med.icb.goby.readers.FastXReader;
import edu.cornell.med.icb.goby.readers.FastXEntry;

/**
 * Describe class here.
 *
 * @author Kevin Dorff
 */
public class TestPreBarcodeMatcher {
    /** Logging. */
    private static final Log LOG = LogFactory.getLog(TestPreBarcodeMatcher.class);

    private final static Random rnd = new Random(System.currentTimeMillis());

    @Test
    public void testMatchBarcodes() throws IOException {
        final long time = System.currentTimeMillis();
        PreBarcodeMatcher matcher = new PreBarcodeMatcher(TestPostBarcodeMatcher.TEST_BARCODES, 5, 2);

        int i = 0;
        int numMatches = 0;
        int numAmbiguous = 0;
        int numNoMatches = 0;
        for (final FastXEntry entry :
                new FastXReader("test-data/barcode/pre-test.fa")) {
            BarcodeMatcherResult bestMatch = matcher.matchSequence(entry.getSequence());
            if (bestMatch != null) {
                int[] result = TestPostBarcodeMatcher.parseTestHeader(entry.getEntryHeader());
                final int sequenceLength = result[0];
                final int barcodeLen = result[1];
                final int whichBarcodeIndex = result[2];
                LOG.info("Testing " + i);
                assertEquals(sequenceLength, bestMatch.getSequenceLength());
                assertEquals(barcodeLen, bestMatch.getBarcodeMatchLength());
                assertEquals(whichBarcodeIndex, bestMatch.getBarcodeIndex());
                if (LOG.isDebugEnabled()) {
                    LOG.info(String.format("Entry %d (Ambuguous? %s)", i, bestMatch.isAmbiguous()));
                    LOG.info(String.format("    Header=%s", entry.getEntryHeader()));
                    LOG.info(String.format("    %s %s",
                            bestMatch.sequenceBarcodeOf(entry.getSequence()),
                            bestMatch.sequenceOf(entry.getSequence())));
                    LOG.info(String.format("    %s (%d)",
                            bestMatch.barcodeWithMatchingAdapterOf(matcher),
                            bestMatch.getNumberOfDiffs()));
                    LOG.info(String.format("    %s",
                            bestMatch.barcodeOnlyOf(matcher)));
                }
                numMatches++;
                if (bestMatch.isAmbiguous()) {
                    numAmbiguous++;
                }
            } else {
                if (LOG.isInfoEnabled()) {
                    LOG.info(String.format("No match for entry %d", i));
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
        long duration = (System.currentTimeMillis() - time) / 1000;
        LOG.info(String.format("Time to parse 8 million reads %d seconds", duration));
    }

    public void createPreBarcodedFile() throws IOException {
        PrintWriter out = null;
        try {
            out = new PrintWriter("test-data/barcode/pre-test.fa");
            out.println("# >sequenceLength:barcodeLen:whichBarcodeIndex");
            final int minBarcodeLen = 5;
            final int maxBarcodeLen = TestPostBarcodeMatcher.TEST_BARCODES[0].length();
            final int randomSpan = maxBarcodeLen - minBarcodeLen + 1;
            final int totalSequenceLength = maxBarcodeLen + 5;
            for (int i = 0; i < 200; i++) {
                int barcodeLen = rnd.nextInt(randomSpan) + minBarcodeLen;
                int sequenceLength = totalSequenceLength - barcodeLen;
                int whichBarcodeIndex = rnd.nextInt(TestPostBarcodeMatcher.TEST_BARCODES.length);
                out.printf(">%d:%d:%d%n", sequenceLength, barcodeLen, whichBarcodeIndex);
                out.print(TestPostBarcodeMatcher.TEST_BARCODES[whichBarcodeIndex].substring(0, barcodeLen));
                for (int j = 0; j < sequenceLength; j++) {
                    out.print(TestPostBarcodeMatcher.ACGT[rnd.nextInt(TestPostBarcodeMatcher.ACGT.length)]);
                }
                out.println();
            }
        } finally {
            IOUtils.closeQuietly(out);
        }
    }

    public MutableString m(final String str) {
        return new MutableString(str);
    }
}