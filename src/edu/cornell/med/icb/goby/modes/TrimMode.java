/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.modes;

import com.google.protobuf.ByteString;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.reads.Reads;
import edu.cornell.med.icb.goby.reads.ReadsReader;
import it.unimi.dsi.bits.BitVector;
import it.unimi.dsi.bits.LongArrayBitVector;
import it.unimi.dsi.fastutil.booleans.BooleanListIterator;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.logging.ProgressLogger;
import it.unimi.dsi.util.IntHyperLogLogCounterArray;
import org.apache.commons.io.IOUtils;
import org.apache.log4j.Logger;

import java.io.IOException;
import java.io.Writer;

/**
 * Analyse a compact-reads file to determine what to trim.
 *
 * @author Fabien Campagne
 *         Date: June 6 2011
 *         Time: 9:40 PM
 */
public class TrimMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "trim";
    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Analyze compact reads to determine what to trim";

    private String inputFilename;
    private String outputFilename;
    private static final Logger LOG = Logger.getLogger(TrimMode.class);
    private int kgramSize = 10;
    private int numScanned = 0;


    /**
     * {@inheritDoc}
     */
    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }


    /**
     * Configure.
     *
     * @param args command line arguments
     * @return this object for chaining
     * @throws java.io.IOException error parsing
     * @throws com.martiansoftware.jsap.JSAPException
     *                             error parsing
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args)
            throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        inputFilename = jsapResult.getString("input");
        outputFilename = jsapResult.getString("output");


        return this;
    }

    BitVector vector = LongArrayBitVector.getInstance();
    double maxCount;

    @Override
    public void execute() throws IOException {

        ReadsReader reader = null;
        Writer writer = null;
        ProgressLogger progress = new ProgressLogger(LOG);
        try {

            reader = new ReadsReader(inputFilename);
            final long numDistinctElements = 10000; // any k-gram is expected at most 1 million times.
            final double estimateStandardDeviationSought = 0.1; // we don't need a lot of precision on the estimate of counts
            final int numberOfCounters = (int) Math.pow(2, (double) kgramSize * 2);
            System.out.println("numberOfCounters " + numberOfCounters);
            final IntHyperLogLogCounterArray counters = new IntHyperLogLogCounterArray(numberOfCounters,
                    IntHyperLogLogCounterArray.registerSize(numDistinctElements),
                    IntHyperLogLogCounterArray.log2NumberOfRegisters(estimateStandardDeviationSought));

            for (Reads.ReadEntry entry : reader) {
                observe(counters, entry.getSequence(), entry.getReadIndex());
            }
        } finally {
            IOUtils.closeQuietly(writer);

        }

        progress.stop();
    }

    private void observe(IntHyperLogLogCounterArray counters, ByteString sequence, int readIndex) {
        ++numScanned;
        int length = sequence.size();
        boolean maxCountChanged = false;
        int maxKgramCode=-1;
        for (int index = 0; index < length - kgramSize; ++index) {
            vector.clear();
            boolean containsN = false;
            for (int j = 0; j < kgramSize; ++j) {
                char c = (char) sequence.byteAt(j);
                switch (c) {
                    case 'A':
                        vector.add(0);
                        vector.add(0);
                        break;
                    case 'C':
                        vector.add(1);
                        vector.add(0);
                        break;
                    case 'T':
                        vector.add(0);
                        vector.add(1);
                        break;
                    case 'G':
                        vector.add(1);
                        vector.add(1);
                        break;
                    default:
                        containsN = true;
                        break;
                }
            }
            if (!containsN) {
                int kgramCode = (int) vector.bits()[0];
                //     System.out.printf("adding %d %d", kgramCode, readIndex);
                counters.add(kgramCode, readIndex);
                final double count = counters.count(kgramCode);
                if (count > maxCount) {
                    maxCount=count;
                    maxCountChanged = true;
                    maxKgramCode = kgramCode;

                }
            }
        }

        if (maxCountChanged && maxCount>10) {

            System.out.printf("count: %g %s read-index %d %n", maxCount, decode(maxKgramCode),readIndex);
        }
    }

    private MutableString decode
            (
                    int kgramCode) {
        LongArrayBitVector vector = LongArrayBitVector.wrap(new long[]{kgramCode});
        final BooleanListIterator it = vector.iterator();
        MutableString sequence = new MutableString();
        int pos = 0;
        while (it.hasNext()) {
            boolean a = it.nextBoolean();
            boolean b = it.nextBoolean();
            if (a) {
                if (b) {
                    sequence.append('G');
                } else {
                    sequence.append('C');
                }

            } else {
                if (b) {
                    sequence.append('T');
                } else {
                    sequence.append('A');
                }
            }
            ++pos;
            if (pos == kgramSize) {
                return sequence;
            }
        }
        return null;
    }

    public static void main
            (
                    final String[] args) throws IOException, JSAPException {
        new TrimMode().configure(args).execute();
    }
}
