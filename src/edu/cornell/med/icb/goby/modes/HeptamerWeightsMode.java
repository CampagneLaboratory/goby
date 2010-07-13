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
import edu.cornell.med.icb.goby.algorithmic.data.HeptamerInfo;
import edu.cornell.med.icb.goby.algorithmic.data.WeightsInfo;
import edu.cornell.med.icb.goby.reads.Reads;
import edu.cornell.med.icb.goby.reads.ReadsReader;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.shorts.ShortArrayList;
import it.unimi.dsi.fastutil.shorts.ShortList;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Evaluate heptamer weights following the method of Hansen et al, NAR April 2010.
 *
 * @author Fabien Campagne
 *         Date: May 4 2009
 *         Time: 12:28 PM
 */
public class HeptamerWeightsMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "heptamer-weights";
    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Estimates heptamer frequencies and weigths. Independent re-implementation of the weight estimates described in Hansen KD et al NAR 2010. See http://campagnelab.org/software/goby/tutorials/estimate-heptamer-weights/";

    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(HeptamerWeightsMode.class);

    private HeptamerInfo heptamers;
    private String[] inputFilenames;
    private String heptamerCountFilename;
    private PrintWriter writer;
    private Int2ObjectOpenHashMap<int[]> heptamerCounts = new Int2ObjectOpenHashMap<int[]>();
    private String tabWeightFilename;
    private int[] heptamerTotalCounts;
    private ShortList readIndexToHeptamerIndex = new ShortArrayList();
    private String mapFilename;
    private String heptamerInfoFilename;
    private boolean colorSpace;


    @Override
    public String getModeName() {
        return MODE_NAME;
    }

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

        inputFilenames = jsapResult.getStringArray("input");
        heptamerCountFilename = jsapResult.getString("heptamer-counts");
        tabWeightFilename = jsapResult.getString("weights");
        mapFilename = jsapResult.getString("map");
        heptamerInfoFilename = jsapResult.getString("heptamer-info");
        colorSpace = jsapResult.getBoolean("color-space");

        return this;
    }


    @Override
    public void execute() throws IOException {
        heptamers = new HeptamerInfo();
        // output file extension is based on the output format type


        final ProgressLogger progress = new ProgressLogger();
        progress.start();
        progress.displayFreeMemory = true;
        final int[] readIndices = {1, 2, -5, -4, -3, -2, -1, 0};
        heptamerTotalCounts = new int[readIndices.length];
        for (final String inputFilename : inputFilenames) {
            LOG.info("Now scanning " + inputFilename);
            final ReadsReader reader = new ReadsReader(new FileInputStream(inputFilename));
            try {
                final int count = 0;
                final MutableString sequence = new MutableString();
                int numberOfReads = 0;
                heptamers.colorSpace = colorSpace;
                for (final Reads.ReadEntry readEntry : reader) {
                    ReadsReader.decodeSequence(readEntry, sequence);
                    if (colorSpace) {
                        sequence.delete(0, 1);

                    }
                    // if (count++ > 100000) break;
                    int item = 0;

                    for (final int positionInRead : readIndices) {
                        final int recodedReadIndex = recodeReadIndex(sequence, positionInRead);

                        final int end = recodedReadIndex - 1 + heptamers.heptamerLength;
                        final int start = recodedReadIndex - 1;
                        //     System.out.printf("%d %d %d %d%n", positionInRead,sequence.length(),start, end );
                        final MutableString heptamer = sequence.substring(start, end);

                        if (heptamer.indexOf('N') == -1) {
                            // heptamers that include any number of Ns are ignored.
                            final short heptamerIndex = (short) heptamers.heptamerToIndices.registerIdentifier(heptamer);
                            accumulateFrequency(heptamerIndex, positionInRead, item, readIndices.length);
                            if (positionInRead == 1) {
                                // this is the heptamer that starts at position 1 of the read,
                                // associate this read index to the heptamer index:
                                final int compactReadIndex = readEntry.getReadIndex();

                                // readIndexToHeptamerIndex.add(heptamerIndex);
                                if (readIndexToHeptamerIndex.size() - 1 < compactReadIndex) {
                                    readIndexToHeptamerIndex.size((compactReadIndex + 10) * 2);
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
            } finally {
                if (reader != null) {
                    try {
                        reader.close();
                    } catch (IOException e) { // NOPMD
                        // silently ignore
                    }
                }
            }
        }
        final DoubleIndexedIdentifier indicesToHeptamer = new DoubleIndexedIdentifier(heptamers.heptamerToIndices);
        PrintWriter tabWeightWriter = null;

        final StringBuffer basenameBuffer = new StringBuffer();
        int last = 0;
        for (final String inputFile : inputFilenames) {
            basenameBuffer.append(FilenameUtils.getBaseName(inputFile));
            if ((++last) != inputFilenames.length) {
                basenameBuffer.append(",");
            }
        }

        LOG.info("Writing frequencies and weights.");

        final String basename = basenameBuffer.toString();
        try {

            writer = heptamerCountFilename == null ? null : new PrintWriter(new FileWriter(heptamerCountFilename));
            if (writer != null) {
                writer.printf("heptamer\t%s\theptamerIndex\treadPositionIndex\theptamerCount%n",
                        totabs(new MutableString("ABCDEFGHIJKLMN")));
            }
            tabWeightWriter = tabWeightFilename == null ? null : new PrintWriter(new FileWriter(tabWeightFilename));
            if (tabWeightWriter != null) {
                tabWeightWriter.printf("basename\theptamer\tweight%n");
            }

            for (short heptamerIndex = 0; heptamerIndex < heptamerCounts.size(); heptamerIndex++) {
                final MutableString heptamer = indicesToHeptamer.getId(heptamerIndex);
                int itemIndex = 0;
                for (final int readIndex : readIndices) {

                    if (writer != null) {
                        writer.printf("%s\t%s\t%d\t%d\t%d%n",
                                heptamer,
                                totabs(heptamer),
                                heptamerIndex,
                                readIndex,
                                heptamerCounts.get(heptamerIndex)[itemIndex++]);
                    }
                }
                final float weight = calculateWeight(heptamerIndex, readIndices);
                if (tabWeightWriter != null) {
                    tabWeightWriter.printf("%s\t%s\t%g%n", basename,
                            heptamer, weight);
                }
                heptamers.heptamerIndexToWeight.put(heptamerIndex, weight);
            }
            LOG.info("writing heptamer info.");
            heptamers.save(heptamerInfoFilename);

            LOG.info("writing binary weight file.");
            final WeightsInfo weights = new WeightsInfo();
            weights.size(readIndexToHeptamerIndex.size());
            for (int readIndex = 0; readIndex < readIndexToHeptamerIndex.size(); readIndex++) {
                final short heptamerIndex = readIndexToHeptamerIndex.get(readIndex);
                weights.setWeight(readIndex, heptamerIndex == -1 ? 1 : heptamers.heptamerIndexToWeight.get(heptamerIndex));
            }

            weights.save(mapFilename);
        } finally {
            IOUtils.closeQuietly(writer);
            IOUtils.closeQuietly(tabWeightWriter);

        }

        progress.stop();
    }

    private String totabs(final MutableString heptamer) {
        final MutableString result = new MutableString();
        for (int i = 0; i < heptamer.length(); i++) {
            result.append(heptamer.charAt(i));
            final boolean notEnd = i < heptamer.length() - 1;
            if (notEnd) {
                result.append('\t');
            }
        }
        return result.toString();
    }


    private int recodeReadIndex(final MutableString sequence, final int readIndex) {
        int recodedReadIndex = readIndex;
        if (readIndex <= 0) {
            recodedReadIndex = sequence.length() - heptamers.heptamerLength + readIndex;
        }
        return recodedReadIndex;
    }

    private float calculateWeight(final int heptamerIndex, final int[] readIndices) {
        float numerator = 0;
        float denominator = 0;
        int itemIndex = 0;
        int numEndPositions = 0;
        int numStartPositions = 0;
        for (final int readIndex : readIndices) {

            final double proportion = divide(heptamerCounts.get(heptamerIndex)[itemIndex], heptamerTotalCounts[itemIndex]);
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

    private double divide(final int numerator, final int denominator) {
        return (double) numerator / (double) denominator;
    }

    private void accumulateFrequency(final int heptamerIndex, final int readIndex, final int i, final int maxItems) {
        int[] countByReadIndex = heptamerCounts.get(heptamerIndex);
        if (countByReadIndex == null) {
            countByReadIndex = new int[maxItems];
            heptamerCounts.put(heptamerIndex, countByReadIndex);
        }

        countByReadIndex[i] += 1;
        heptamerTotalCounts[i] += 1;
    }


    public static void main(final String[] args) throws IOException, JSAPException {
        new HeptamerWeightsMode().configure(args).execute();
    }
}
