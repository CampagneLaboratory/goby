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

package edu.cornell.med.icb.goby.modes;

import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.reads.Reads;
import edu.cornell.med.icb.goby.reads.ReadsReader;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.io.FastBufferedOutputStream;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.io.IOUtils;

import java.io.*;

/**
 * Evaluate heptamer weights following the method of Hansen et al, NAR April 2010.
 *
 * @author Fabien Campagne
 *         Date: May 4 2009
 *         Time: 12:28 PM
 */
public class HeptamerFrequenciesMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "heptamer-frequencies";
    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Estimates heptamer frequencies.";


    private String inputFilenames[];
    private String outputFilename;
    private PrintWriter writer;
    private int heptamerLength = 7;
    private Int2ObjectOpenHashMap<int[]> heptamerCounts = new Int2ObjectOpenHashMap<int[]>();
    private String weightFilename;
    private int[] heptamerTotalCounts;


    public String getModeName() {
        return MODE_NAME;
    }

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

        inputFilenames = jsapResult.getStringArray("input");
        outputFilename = jsapResult.getString("output");
        weightFilename = jsapResult.getString("weights");


        return this;
    }


    @Override
    public void execute() throws IOException {
        // output file extension is based on the output format type
        IndexedIdentifier heptamerToIndices = new IndexedIdentifier();

        final ProgressLogger progress = new ProgressLogger();
        progress.start();
        progress.displayFreeMemory = true;
        int readIndices[] = {1, 2, -5, -4, -3, -2, -1, 0};
        heptamerTotalCounts = new int[readIndices.length];
        for (String inputFilename : inputFilenames) {
            ReadsReader reader = new ReadsReader(new FileInputStream(inputFilename));
            try {
                int count = 0;
                final MutableString sequence = new MutableString();

                for (final Reads.ReadEntry readEntry : reader) {
                    ReadsReader.decodeSequence(readEntry, sequence);
                    //   if (count++ > 10000) break;
                    int item = 0;
                    for (int readIndex : readIndices) {
                        int recodedReadIndex = recodeReadIndex(sequence, readIndex);

                        final int end = recodedReadIndex - 1 + heptamerLength;
                        final int start = recodedReadIndex - 1;
                        //     System.out.printf("%d %d %d %d%n", readIndex,sequence.length(),start, end );
                        final MutableString heptamer = sequence.substring(start, end);
                        if (heptamer.indexOf('N') == -1) {
                            int heptamerIndex = heptamerToIndices.registerIdentifier(heptamer);
                            accumulateFrequency(heptamerIndex, readIndex, item++, readIndices.length);
                        }
                    }
                    progress.lightUpdate();
                }

            } finally

            {
                if (reader != null) {
                    try {
                        reader.close();
                    } catch (IOException e) { // NOPMD
                        // silently ignore
                    }
                }
            }
        }
        DoubleIndexedIdentifier indicesToHeptamer = new DoubleIndexedIdentifier(heptamerToIndices);
        PrintWriter weightWriter = null;
        try {
            writer = new PrintWriter(new FileWriter(outputFilename));
            weightWriter = new PrintWriter(new FileWriter(weightFilename));
            for (int heptamerIndex = 0; heptamerIndex < heptamerCounts.size(); heptamerIndex++) {
                final MutableString heptamer = indicesToHeptamer.getId(heptamerIndex);
                int itemIndex = 0;
                for (int readIndex : readIndices) {

                    writer.printf("%s\t%d\t%d\t%d%n",
                            heptamer,
                            heptamerIndex,
                            readIndex,
                            heptamerCounts.get(heptamerIndex)[itemIndex++]);
                }
                weightWriter.printf("%s\t%g%n", heptamer, calculateWeight(heptamerIndex, readIndices));
            }
            writer.close();
        }


        finally {
            IOUtils.closeQuietly(writer);
            IOUtils.closeQuietly(weightWriter);

        }

        progress.stop();
    }

    private int recodeReadIndex(MutableString sequence, int readIndex) {
        int recodedReadIndex = readIndex;
        if (readIndex <= 0) {
            recodedReadIndex = sequence.length() - heptamerLength + readIndex;
        }
        return recodedReadIndex;
    }

    private double calculateWeight(int heptamerIndex, int readIndices[]) {
        double numerator = 0;
        double denominator = 0;
        int itemIndex = 0;
        int numEndPositions = 0;
        int numStartPositions = 0;
        for (int readIndex : readIndices) {

            double proportion = divide(heptamerCounts.get(heptamerIndex)[itemIndex], heptamerTotalCounts[itemIndex]);
            if (readIndex <= 0) {
                numerator += proportion;
                numEndPositions++;
            } else {
                denominator += proportion;
                numStartPositions++;
            }
            itemIndex++;
        }
        numerator /= numEndPositions;
        denominator /= numStartPositions;
        return numerator / denominator;

    }

    private double divide(int numerator, int denominator) {
        return (double) numerator / (double) denominator;
    }

    private void accumulateFrequency(int heptamerIndex, int readIndex, int i, int maxItems) {
        int[] countByReadIndex = heptamerCounts.get(heptamerIndex);
        if (countByReadIndex == null) {
            countByReadIndex = new int[maxItems];
            heptamerCounts.put(heptamerIndex, countByReadIndex);
        }

        countByReadIndex[i] += 1;
        heptamerTotalCounts[i] += 1;
    }


    public static void main(final String[] args) throws IOException, JSAPException {
        new HeptamerFrequenciesMode().configure(args).execute();
    }
}