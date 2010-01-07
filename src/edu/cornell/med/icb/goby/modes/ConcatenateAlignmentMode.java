/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                         Weill Medical College of Cornell University
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
import edu.cornell.med.icb.goby.aligners.AbstractAligner;
import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.AlignmentWriter;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.alignments.ConcatAlignmentReader;
import edu.cornell.med.icb.goby.alignments.Merge;
import it.unimi.dsi.logging.ProgressLogger;

import java.io.File;
import java.io.IOException;

/**
 * Concatenate compact alignment files.
 * <p/>
 * Reference sequences must match exactly across the input alignments.
 * Query are assumed to be entirely distinct and will be treated as independent observations (e.g.,
 * reads from multiple independent samples). To this effect, alignment entries read from
 * different input basenames, which would otherwise share an identical query index,
 * are renumbered with distinct query indices.
 *
 * @author Fabien Campagne
 *         Date: Apr 28, 2009
 *         Time: 6:03:56 PM
 */
public class ConcatenateAlignmentMode extends AbstractGobyMode {

    /**
     * The mode name.
     */
    public static final String MODE_NAME = "concatenate-alignments";
    public static final String MODE_DESCRIPTION = "Concatenate compact alignment files.  <p/> Reference sequences must match exactly across the input alignments.  Query are assumed to be entirely distinct and will be treated as independent observations (e.g., reads from multiple independent samples). To this effect, alignment entries read from different input basenames, which would otherwise share an identical query index, are renumbered with distinct query indices.";

    private String[] inputFilenames;
    private String outputFile;

    private int sequencePerOutput = Integer.MAX_VALUE;

    // default is true.
    private boolean adjustQueryIndices=true;

    public String getModeName() {
        return MODE_NAME;
    }

    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }

    public void setAdjustQueryIndices(boolean adjustQueryIndices) {
        this.adjustQueryIndices = adjustQueryIndices;
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
        adjustQueryIndices = jsapResult.getBoolean("adjust-query-indices",true);
        return this;
    }

    /**
     * Perform the concatenation.
     *
     * @throws java.io.IOException
     */
    @Override
    public void execute() throws IOException {


        final String outputFilename                 = outputFile;
        final AlignmentWriter writer                = new AlignmentWriter(outputFilename);
        final String[] basenames                    = AlignmentReader.getBasenames(inputFilenames);
        final ConcatAlignmentReader alignmentReader = new ConcatAlignmentReader(basenames);
        final ProgressLogger progress               = new ProgressLogger();

        alignmentReader.setAdjustQueryIndices(adjustQueryIndices);

        int entriesInOutputFile = 0;
        long numLogicalEntries  = 0;
        long numEntries         = 0;
        int numQueries          = alignmentReader.getNumberOfQueries();

        progress.start("concatenating entries");

        for (final Alignments.AlignmentEntry entry : alignmentReader) {

            writer.appendEntry(entry);
            numLogicalEntries += entry.getMultiplicity();
            numEntries++;
            entriesInOutputFile++;
            progress.lightUpdate();

        }
        alignmentReader.getStatistics();
        progress.stop();
        // too many hits is prepared as for Merge:
        Merge.prepareMergedTooManyHits(outputFile, alignmentReader.getNumberOfQueries(), basenames);
        writer.setNumQueries(alignmentReader.getNumberOfQueries());
        writer.setNumTargets(alignmentReader.getNumberOfTargets());
        if (alignmentReader.getTargetIdentifiers() != null) {
            writer.setTargetIdentifiers(alignmentReader.getTargetIdentifiers());
        }
        if (alignmentReader.getQueryIdentifiers() != null) {
            writer.setQueryIdentifiers(alignmentReader.getQueryIdentifiers());
        }
        writer.setStatistics(alignmentReader.getStatistics());
        writer.putStatistic("overall.matched.percent", String.format("%3.3g",divide(numLogicalEntries,numQueries ) * 100d));
        writer.close();

        writer.printStats(System.out);
        System.out.printf("Wrote a total of %d alignment entries.%n", entriesInOutputFile);
    }

    private double divide(long a, long b) {
        return ((double) a) / ((double) b);
    }


    public static void main(final String[] args) throws IOException, JSAPException {
        new ConcatenateAlignmentMode().configure(args).execute();
    }


    public void setInputFileNames(String[] inputFilenames) {
      this.  inputFilenames=inputFilenames;
    }

    public void setOutputFilename(String output) {
        this.outputFile=output;
    }

    public File[] getOutputFiles() {
        return AbstractAligner.buildResults(outputFile);
    }
}
