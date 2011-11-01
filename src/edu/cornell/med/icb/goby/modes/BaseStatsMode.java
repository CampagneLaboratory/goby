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

import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.reads.Reads;
import edu.cornell.med.icb.goby.reads.ReadsReader;
import edu.cornell.med.icb.goby.reads.ReadsWriter;
import edu.cornell.med.icb.goby.util.motifs.MotifMatcher;
import edu.cornell.med.icb.goby.util.motifs.SubSequenceMotifMatcher;
import it.unimi.dsi.fastutil.io.FastBufferedOutputStream;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.log4j.Logger;

import java.io.*;

/**
 * Provide base frequency statistics for a set of reads.
 *
 * @author Fabien Campagne
 *         Date: Sept 30 2011
 *         Time: 12:28 PM
 */
public class BaseStatsMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "base-stats";
    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Calculates frequencies for bases and motifs for a set of reads.";


    private String outputFilename;
    private boolean doCpX;
    /**
     * Set of input filenames.
     */
    private String[] inputFilenames;


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

    private static final Logger LOG = Logger.getLogger(BaseStatsMode.class);

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
        doCpX = jsapResult.getBoolean("cpX");
        return this;
    }

    @Override
    public void execute() throws IOException {


        final PrintStream writer = outputFilename == null || "-".equals(outputFilename) ?
                System.out : new PrintStream(new FileOutputStream(outputFilename));
        MutableString sequence = new MutableString();
        for (String filename: inputFilenames) {

            ReadsReader reader = new ReadsReader(filename);
// an array of longs that will hold the count in elements 'A' 'C' 'G' T'
            long[] tallies = new long['T' + 1];
            ProgressLogger progress = new ProgressLogger(LOG);
            progress.start("Starting to analyze base stats.");
            progress.itemsName = "reads";
            // stores frequencies of the CpX motifs
            long[] cpFreq = new long[4];

            final MotifMatcher[] matchers = {
                    new SubSequenceMotifMatcher("CG"),
                    new SubSequenceMotifMatcher("CA"),
                    new SubSequenceMotifMatcher("CT"),
                    new SubSequenceMotifMatcher("CC")
            };

            for (Reads.ReadEntry entry : reader) {
                for (final MotifMatcher matcher : matchers) {
                    matcher.newSequence();
                }
                ReadsReader.decodeSequence(entry, sequence);
                for (int i = 0; i < sequence.length(); i++) {
                    char c = sequence.charAt(i);
                    tallies[c]++;
                    int index = 0;
                    if (doCpX) {
                        for (final MotifMatcher matcher : matchers) {
                            matcher.accept(c);
                            if (matcher.matched()) {
                                ++cpFreq[index];
                            }
                            index++;
                        }
                    }
                }

                progress.lightUpdate();

            }
            progress.done();

            long sum = tallies['A'] + tallies['C'] + tallies['T'] + tallies['G'];
            writer.printf(
                    "BASES\t%s\tA:\t%d\t%g%%%n" +
                    "BASES\t%s\tC:\t%d\t%g%% %n" +
                    "BASES\t%s\tT:\t%d\t%g%% %n" +
                    "BASES\t%s\tG:\t%d\t%g%% %n",
                    filename,tallies['A'], percent(tallies['A'], sum),
                    filename,tallies['C'], percent(tallies['C'], sum),
                    filename,tallies['T'], percent(tallies['T'], sum),
                    filename,tallies['G'], percent(tallies['G'], sum)
            );
            if (doCpX) {
                long sumcpX = cpFreq[0] + cpFreq[1] + cpFreq[2] + cpFreq[3];
                writer.printf(
                        "CPX\t%s\tCpG:\t%d\t%g%%\t%n" +
                        "CPX\t%s\tCpA:\t%d\t%g%%\t%n" +
                        "CPX\t%s\tCpT:\t%d\t%g%%\t%n" +
                        "CPX\t%s\tCpC:\t%d\t%g%%\t%n" +
                        "CPX\t%s\tCpT+CpC:\t%d\t%g%%\t%n",
                        filename,cpFreq[0], percent(cpFreq[0], sumcpX),
                        filename,cpFreq[1], percent(cpFreq[1], sumcpX),
                        filename,cpFreq[2], percent(cpFreq[2], sumcpX),
                        filename,cpFreq[3], percent(cpFreq[3], sumcpX),
                        filename,cpFreq[2] + cpFreq[3], percent(cpFreq[2] + cpFreq[3], sumcpX)
                );
            }

            writer.flush();
            progress.stop();

        }
    }

    private static double percent
            (
                    double tally,
                    double sum) {
        return tally / sum * 100.0;
    }

    public static void main
            (
                    final String[] args) throws IOException, JSAPException {
        new BaseStatsMode().configure(args).execute();
    }
}