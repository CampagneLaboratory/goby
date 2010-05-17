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
import it.unimi.dsi.fastutil.io.BinIO;
import it.unimi.dsi.fastutil.shorts.ShortList;
import it.unimi.dsi.fastutil.shorts.ShortArrayList;
import it.unimi.dsi.fastutil.shorts.Short2FloatMap;
import it.unimi.dsi.fastutil.shorts.Short2FloatOpenHashMap;
import it.unimi.dsi.fastutil.floats.FloatArrayList;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.io.IOUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

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

    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(HeptamerFrequenciesMode.class);


    private String inputFilenames[];
    private String outputFilename;
    private PrintWriter writer;
    private int heptamerLength = 7;
    private Int2ObjectOpenHashMap<int[]> heptamerCounts = new Int2ObjectOpenHashMap<int[]>();
    private String weightFilename;
    private int[] heptamerTotalCounts;
    private ShortList readIndexToHeptamerIndex = new ShortArrayList();
    private String mapFilename;


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
        mapFilename = jsapResult.getString("map");


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
                int numberOfReads = 0;
                for (final Reads.ReadEntry readEntry : reader) {
                    ReadsReader.decodeSequence(readEntry, sequence);
                    // if (count++ > 1000000) break;
                    int item = 0;

                    for (int positionInRead : readIndices) {
                        int recodedReadIndex = recodeReadIndex(sequence, positionInRead);

                        final int end = recodedReadIndex - 1 + heptamerLength;
                        final int start = recodedReadIndex - 1;
                        //     System.out.printf("%d %d %d %d%n", positionInRead,sequence.length(),start, end );
                        final MutableString heptamer = sequence.substring(start, end);

                        if (heptamer.indexOf('N') == -1) {
                            // heptamers that include any number of Ns are ignored.
                            short heptamerIndex = (short) heptamerToIndices.registerIdentifier(heptamer);
                            accumulateFrequency(heptamerIndex, positionInRead, item, readIndices.length);
                            if (positionInRead == 1) {
                                // this is the heptamer that starts at position 1 of the read,
                                // associate this read index to the heptamer index:
                                int compactReadIndex = readEntry.getReadIndex();

                                // readIndexToHeptamerIndex.add(heptamerIndex);
                                if (readIndexToHeptamerIndex.size() - 1 < compactReadIndex) {
                                    readIndexToHeptamerIndex.size((readIndexToHeptamerIndex.size() + 10) * 2);
                                }

                                readIndexToHeptamerIndex.set(compactReadIndex, heptamerIndex);
                            }

                        }
                        item++;


                    }

                    progress.lightUpdate();
                    numberOfReads++;
                }
                readIndexToHeptamerIndex.size(numberOfReads);
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
        Short2FloatMap heptamerIndexToWeight = new Short2FloatOpenHashMap();
        StringBuffer basenameBuffer = new StringBuffer();
        int last = 0;
        for (String inputFile : inputFilenames) {
            basenameBuffer.append(FilenameUtils.getBaseName(inputFile));
            if ((++last) != inputFilenames.length) basenameBuffer.append(",");
        }

        LOG.info("Writing frequencies and weights.");

        String basename = basenameBuffer.toString();
        try {
            writer = new PrintWriter(new FileWriter(outputFilename));
            weightWriter = new PrintWriter(new FileWriter(weightFilename));
            for (short heptamerIndex = 0; heptamerIndex < heptamerCounts.size(); heptamerIndex++) {
                final MutableString heptamer = indicesToHeptamer.getId(heptamerIndex);
                int itemIndex = 0;
                for (int readIndex : readIndices) {

                    writer.printf("%s\t%d\t%d\t%d%n",
                            heptamer,
                            heptamerIndex,
                            readIndex,
                            heptamerCounts.get(heptamerIndex)[itemIndex++]);
                }
                final float weight = calculateWeight(heptamerIndex, readIndices);
                weightWriter.printf("%s\t%s\t%g%n", basename,
                        heptamer, weight);
                heptamerIndexToWeight.put(heptamerIndex, weight);
            }
            LOG.info("writting binary weight file.");
            FloatArrayList weights = new FloatArrayList();
            weights.size(readIndexToHeptamerIndex.size());
            for (int readIndex = 0; readIndex < readIndexToHeptamerIndex.size(); readIndex++) {
                final short heptamerIndex = readIndexToHeptamerIndex.get(readIndex);
                weights.set(readIndex, heptamerIndex == -1 ? 1 : heptamerIndexToWeight.get(heptamerIndex));
            }

            BinIO.storeObject(weights, mapFilename);
        }


        finally {
            IOUtils.closeQuietly(writer);
            IOUtils.closeQuietly(weightWriter);

        }

        progress.stop();
    }

    private String totabs(MutableString heptamer) {
        MutableString result = new MutableString();
        for (int i = 0; i < heptamer.length(); i++) {
            result.append(heptamer.charAt(i));
            result.append('\t');
        }
        return result.toString();
    }


    private int recodeReadIndex(MutableString sequence, int readIndex) {
        int recodedReadIndex = readIndex;
        if (readIndex <= 0) {
            recodedReadIndex = sequence.length() - heptamerLength + readIndex;
        }
        return recodedReadIndex;
    }

    private float calculateWeight(int heptamerIndex, int readIndices[]) {
        float numerator = 0;
        float denominator = 0;
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