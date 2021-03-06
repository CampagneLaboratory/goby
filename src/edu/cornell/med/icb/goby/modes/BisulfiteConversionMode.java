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

package edu.cornell.med.icb.goby.modes;

import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.reads.ReadsWriterImpl;
import it.unimi.dsi.fastutil.io.FastBufferedOutputStream;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.logging.ProgressLogger;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

import edu.cornell.med.icb.goby.reads.ReadsWriter;
import edu.cornell.med.icb.goby.reads.ReadsReader;
import edu.cornell.med.icb.goby.reads.Reads;

/**
 * Convert sequences to mimick bisulfite conversion. Every cytosine is assumed to be unmethylated and is therefore
 * converted to thymine.
 *
 * @author Fabien Campagne
 *         Date: Oct 20 0210
 *         Time: 12:28 PM
 */
public class BisulfiteConversionMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "bisulfite-conversion";
    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Treat DNA sequences with bisulfite, in silico. The output compact " +
            "reads file contains sequences where cytosine bases have been substituted with thimine (C->T conversion)." +
            "Please note that this mode does not support color space at this time. ";


    private String inputFilename;
    private String outputFilename;

    /**
     * The number of sequences that will be written in each compressed chunk. Th default is
     * suitable for very many short sequences but should be reduced to a few sequences per
     * chunk if each sequence is very large.
     */
    private int sequencePerChunk = 10000;

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
        sequencePerChunk = jsapResult.getInt("sequence-per-chunk");

        return this;
    }

    @Override
    public void execute() throws IOException {

        System.out.printf("Converting %s to %s, sequences per chuck=%d%n",
                inputFilename, outputFilename, sequencePerChunk);

        final ProgressLogger progress = new ProgressLogger();
        progress.start();
        progress.displayFreeMemory = true;
        ReadsWriter writer = null;
        ReadsReader reader = null;
        try {

            writer = new ReadsWriterImpl(new FastBufferedOutputStream(new FileOutputStream(outputFilename)));
            writer.setNumEntriesPerChunk(sequencePerChunk);

            final MutableString sequence = new MutableString();
            final MutableString sequencePair = new MutableString();
            byte[] byteBuffer = new byte[0];

            reader = new ReadsReader(new FileInputStream(inputFilename));
            for (final Reads.ReadEntry readEntry : reader) {
                Reads.ReadEntry.Builder builder = readEntry.toBuilder();
                ReadsReader.decodeSequence(readEntry, sequence, false);
                if (sequence.length() > byteBuffer.length) {
                    byteBuffer = new byte[sequence.length()];
                }
                convert(sequence);
                builder.setSequence(ReadsWriterImpl.encodeSequence(sequence, byteBuffer));

                if (readEntry.hasSequencePair()) {
                    ReadsReader.decodeSequence(readEntry, sequencePair, true);
                    if (sequencePair.length() > byteBuffer.length) {
                        byteBuffer = new byte[sequencePair.length()];
                    }
                    convert(sequencePair);
                    builder.setSequencePair(ReadsWriterImpl.encodeSequence(sequencePair, byteBuffer));

                }
                writer.appendEntry(builder);

                progress.lightUpdate();
            }
        } finally {

            if (writer != null) {
                try {
                    writer.close();
                } catch (IOException e) { // NOPMD
                    // silently ignore
                }
            }
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e) { // NOPMD
                    // silently ignore
                }
            }
        }

        progress.stop();
    }

    // Convert each C to T
    private void convert(final MutableString sequence) {
        final int length = sequence.length();
        for (int i = 0; i < length; i++) {
            char base = sequence.charAt(i);
            if (base == 'C') base = 'T';

            sequence.setCharAt(i, base);
        }

    }

    public static void main(final String[] args) throws IOException, JSAPException {
        new BisulfiteConversionMode().configure(args).execute();
    }
}