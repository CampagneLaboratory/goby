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
import edu.cornell.med.icb.goby.util.SimulateBisulfiteReads;
import edu.mssm.crover.cli.CLI;

import java.io.File;
import java.io.IOException;

/**
 * Generate simulated reads.
 *
 * @author Fabien Campagne
 */
public class SimulateReadsMode extends AbstractGobyMode {

    /**
     * The mode name.
     */
    private static final String MODE_NAME = "simulate-reads";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Generate simulated reads to test other applications.";

     @Override
    public String getModeName() {
        return MODE_NAME;
    }

    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }
    
    private SimulateBisulfiteReads processor;
    private String fastaReference;
    private String refChoice;
    private int from;
    private int to;
    private String methylationRateFilename;




    /**
     * Configure the mode arguements.
     *
     * @param args the arguments
     * @return this object for chaining
     * @throws java.io.IOException error parsing arguments
     * @throws com.martiansoftware.jsap.JSAPException
     *                             error parsing arguments
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args) throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        fastaReference = jsapResult.getString("input");
        String outputFilename = jsapResult.getString("output");
        String regionTrueRates = jsapResult.getString("methylation-rates");
        // reference sequence to use
        refChoice = jsapResult.getString( "chromosome");
        from = jsapResult.getInt( "start");
        to = jsapResult.getInt( "end");
        if (to < from) {
            System.err.println("argument to -e must be larger than argument to -s");
            System.exit(1);
        }
        int readLength = jsapResult.getInt("read-length");
        methylationRateFilename = jsapResult.getString("methylation-rates");
        final boolean bisulfite = jsapResult.getBoolean("bisulfite");
        String strandChoice = jsapResult.getString("strand");

        processor = new SimulateBisulfiteReads();

        processor.configure(bisulfite, strandChoice);

        processor.setBisulfiteTreatment(bisulfite);
        processor.setReadLength(readLength);
        processor.setOutputFilename(outputFilename);
        processor.setRegionTrueRates(regionTrueRates);
        processor.setNumRepeats(jsapResult.getInt("num-repeats"));

        return this;
    }

    /**
     * Split a fasta / fastq file by (a) readlength and (b) the maximum number of
     * entries per file. This will output the files that are written to stdout
     *
     * @throws java.io.IOException error reading / writing files.
     */
    @Override
    public void execute() throws IOException {

        processor.process(refChoice, fastaReference, from, to, methylationRateFilename);
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
        new SimulateReadsMode().configure(args).execute();
    }
}