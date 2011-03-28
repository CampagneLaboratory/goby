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
import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.IterateSortedAlignmentsToPileup;
import org.apache.log4j.Logger;

import java.io.IOException;
import java.io.PrintWriter;

/**
 * This mode writes a region of an alignment as a sequence alignemnt in FASTA or other format.
 *
 * @author Fabien Campagne
 *         Date: January 28 2011
 */
public class AlignmentToPileupMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "alignment-to-pileup";
    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION =
            "This mode writes a region of an alignment as a sequence alignemnt in FASTA or other outputFormat. (Since Goby 1.9.2).";

    private static final Logger LOG = Logger.getLogger(AlignmentToPileupMode.class);
    private String[] inputFilenames;
    private String outputFile;
    private int[] readerIndexToGroupIndex;
    private int numberOfGroups;
    private CharSequence currentReferenceId;
    private int thresholdDistinctReadIndices = 3;
    private int minimumVariationSupport = 10;
    private PrintWriter outWriter;

    private String[] groups;
    /**
     * The maximum value of read index, indexed by readerIndex;
     */
    private int numberOfReadIndices[];
    private int startFlapSize;

    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }

    OutputFormat outputFormat;

    public void setInputFilename(String basename1) {
        inputFilenames = new String[]{basename1};
    }

    public void setOutputWriter(PrintWriter outWriter) {
        this.outWriter = outWriter;
    }

    public void setFormat(OutputFormat format) {
        outputFormat = format;
    }

    

    public enum OutputFormat {
        FASTA,
        ONE_PER_LINE
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
        inputFilenames = jsapResult.getStringArray("input");

        outputFile = jsapResult.getString("output");
        outWriter = "-".equals(outputFile) ? new PrintWriter(System.out) : new PrintWriter(outputFile);
        outputFormat = OutputFormat.valueOf(jsapResult.getString("format").toUpperCase());
        readerIndexToGroupIndex = new int[inputFilenames.length];

        startFlapSize = jsapResult.getInt("start-flap-size");
        sortedPositionIterator = new IterateSortedAlignmentsToPileup();
        sortedPositionIterator.setStartFlapLength(startFlapSize);
        sortedPositionIterator.parseIncludeReferenceArgument(jsapResult);

        return this;
    }

    public void initializeIterator(String startPositionDefinition,String endPositionDefinition) {
        sortedPositionIterator = new IterateSortedAlignmentsToPileup();
        sortedPositionIterator.setStartFlapLength(startFlapSize);
        sortedPositionIterator.setStartPositionArgument(startPositionDefinition);
        sortedPositionIterator.setEndPositionArgument(endPositionDefinition);

    }

    IterateSortedAlignmentsToPileup sortedPositionIterator;


    /**
     * Perform the job.
     *
     * @throws java.io.IOException
     */
    @Override
    public void execute() throws IOException {
        final String outputFilename = outputFile;

        final String[] basenames = AlignmentReader.getBasenames(inputFilenames);
        final boolean allSorted = ConcatenateAlignmentMode.isAllSorted(basenames);
        if (!allSorted) {
            System.out.println("Each input alignment must be sorted. Aborting.");
            System.exit(10);
        }


        sortedPositionIterator.initialize(outWriter, inputFilenames, startFlapSize);

        sortedPositionIterator.iterate(basenames);
        sortedPositionIterator.finish(null, outputFormat);
        outWriter.close();
    }


    /**
     * Main method.
     *
     * @param args command line args.
     * @throws com.martiansoftware.jsap.JSAPException
     *                             error parsing
     * @throws java.io.IOException error parsing or executing.
     */
    public static void main(final String[] args) throws JSAPException, IOException {

        new AlignmentToPileupMode().configure(args).execute();
        System.exit(0);
    }

}