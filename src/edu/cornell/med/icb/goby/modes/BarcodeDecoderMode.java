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

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;

import java.io.IOException;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.util.prefs.Preferences;

import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.io.TSVReader;
import edu.cornell.med.icb.goby.reads.ReadsReader;
import edu.cornell.med.icb.goby.reads.Reads;
import edu.cornell.med.icb.goby.reads.ReadsWriter;
import edu.cornell.med.icb.goby.util.barcode.BarcodeMatcher;
import edu.cornell.med.icb.goby.util.barcode.PostBarcodeMatcher;
import edu.cornell.med.icb.goby.util.barcode.BarcodeMatcherResult;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectLinkedOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.lang.MutableString;

/**
 * @author Fabien Campagne
 *         Date: May 18, 2010
 *         Time: 3:10:59 PM
 */
public class BarcodeDecoderMode extends AbstractGobyMode {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(BarcodeDecoderMode.class);

    /**
     * The mode name.
     */
    private static final String MODE_NAME = "barcode-decoder";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Process a compact reads file and decode barcodes. " +
            "Will either produce a large compact-reads file where matching reads are copied (in this " +
            "case the barcodeIndex attribute is set appropriately on each read), or will produce a " +
            "set of compact reads files, one for each sample indicates in barcode-info. The second " +
            "option is used if no output file is provided on the command line.";
    /**
     * The compact reads file to process.
     */
    private String inputFilename;
    /**
     * The filename to the tab delimited file with barcodes and sample ids.
     */
    private String barcodeInfoFilename;
    private Int2ObjectMap<String> barcodeIndexToSampleId;
    private String[] barcodes;
    private int maxMismatches = 0;

    private String outputFilename;
    private int minimalMatchLength;

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
    public AbstractCommandLineMode configure(final String[] args) throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        inputFilename = jsapResult.getString("input");
        outputFilename = jsapResult.getString("output");
        barcodeInfoFilename = jsapResult.getString("barcode-info");
        minimalMatchLength = jsapResult.getInt("minimal-match-length");
        maxMismatches = jsapResult.getInt("max-mismatches");

        return this;
    }

    public void execute() throws IOException {
        loadBarcodeInfo(barcodeInfoFilename);


        //   assert inputFilenames.length == 1 : "only one read file supported for now.";
        ReadsWriter singleWriter = null;
        ReadsWriter writers[] = new ReadsWriter[barcodeIndexToSampleId.size()];
        if (outputFilename == null) {
            for (int i = 0; i < writers.length; i++) {
                writers[i] = new ReadsWriter(new FileOutputStream(barcodeIndexToSampleId.get(i).trim() + ".compact-reads"));
            }
        } else {
            singleWriter = new ReadsWriter(new FileOutputStream(outputFilename));
        }

        String inputReadsFilename = inputFilename;
        BarcodeMatcher matcher = new PostBarcodeMatcher(barcodes, minimalMatchLength, maxMismatches);
        MutableString sequence = new MutableString();
        MutableString sequenceNoBarcode = new MutableString();
        try {
            int countMatched = 0;
            int countNoMatch = 0;
            int countAmbiguous = 0;
            for (Reads.ReadEntry readEntry : new ReadsReader(inputReadsFilename)) {
                ReadsReader.decodeSequence(readEntry, sequence);
                BarcodeMatcherResult match = matcher.matchSequence(sequence);
                if (match != null) {
                    // remove the barcode from the sequence:
                    int readIndex = readEntry.getReadIndex();
                    sequenceNoBarcode.setLength(0);
                    sequenceNoBarcode.append(sequence.subSequence(0, match.getBarcodeStartPosition() - 1));
                    final int barcodeIndex = match.getBarcodeIndex();
                    if (match.isAmbiguous()) {
                        ++countAmbiguous;
                    }
                    /* System.out.printf("read %d barcode match: %n%s%n%s%n",
                           readIndex,
                           sequence.subSequence(match.getBarcodeStartPosition(),
                                   sequence.length()),
                           barcodes[barcodeIndex]);
                    */
                    ReadsWriter writer = outputFilename == null ? writers[barcodeIndex] : singleWriter;
                    writer.setSequence(sequenceNoBarcode);
                    writer.setBarcodeIndex(barcodeIndex);

                    if (readEntry.hasDescription()) {
                        writer.setDescription(readEntry.getDescription());
                    }
                    if (readEntry.hasReadIdentifier()) {
                        writer.setIdentifier(readEntry.getReadIdentifier());
                    }
                    if (readEntry.hasQualityScores()) {

                        writer.setQualityScores(readEntry.getQualityScores().toByteArray());
                    }
                    writer.appendEntry(readEntry.getReadIndex());
                    ++countMatched;
                } else {
                    ++countNoMatch;

                }
            }
            System.out.format("barcode found in %g %% of the reads %n", percent(countMatched, countMatched + countNoMatch));
            System.out.format("Found %g %% ambiguous matches %n", percent(countAmbiguous, countMatched));
        }
        finally {
            for (int i = 0; i < writers.length; i++) {
                if (writers[i] != null) {
                    writers[i].close();
                }
            }
            if (outputFilename != null) {
                singleWriter.close();
            }

        }

    }

    private double percent(int countMatched, int total) {
        return (double) countMatched / (double) total * 100d;
    }

    private void loadBarcodeInfo
            (String
                    barcodeInfoFilename) {
        try {
            ObjectArrayList<String> barcodes = new ObjectArrayList<String>();
            barcodeIndexToSampleId = new Int2ObjectOpenHashMap<String>();
            TSVReader reader = new TSVReader(new FileReader(barcodeInfoFilename), '\t');
            while (reader.hasNext()) {
                reader.next();
                String sampleId = reader.getString();
                int barcodeIndex = reader.getInt();
                String barcode = reader.getString();
                barcodes.add(barcode);
                barcodeIndexToSampleId.put(barcodeIndex, sampleId);
            }
            this.barcodes = barcodes.toArray(new String[barcodes.size()]);
        } catch (Exception e) {
            System.err.println("Cannot load barcode information from file " + barcodeInfoFilename);
            System.exit(1);
        }
    }


    /**
     * Main method.
     *
     * @param args command line args.
     * @throws com.martiansoftware.jsap.JSAPException
     *                             error parsing
     * @throws java.io.IOException error parsing or executing.
     */
    public static void main
            (
                    final String[] args) throws JSAPException, IOException {
        new BarcodeDecoderMode().configure(args).execute();
    }
}
