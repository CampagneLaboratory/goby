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
import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.alignments.IterateAlignments;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Display the sequence variations found in alignments.
 *
 * @author Fabien Campagne
 */
public class DisplaySequenceVariationsMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "display-sequence-variations";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Display the sequence variations found in an alignment";

    /**
     * The input filenames.
     */
    private String[] inputFilenames;

    /**
     * The output file.
     */
    private String outputFilename;
    /**
     * The input basenames.
     */
    private String[] basenames;
    private MyIterateAlignments alignmentIterator;

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
    public AbstractCommandLineMode configure(final String[] args) throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        inputFilenames = jsapResult.getStringArray("input");
        basenames = AlignmentReader.getBasenames(inputFilenames);
        outputFilename = jsapResult.getString("output");

        alignmentIterator = new MyIterateAlignments();
        alignmentIterator.parseIncludeReferenceArgument(jsapResult);
        return this;
    }

    /**
     * Display sequence variations.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        final PrintWriter writer = outputFilename == null ? new PrintWriter(System.out) :
                new PrintWriter(new FileWriter(outputFilename));
        writer.println("# query-index target-id [position-on-reference:var-from/var-to,]+");
        try {
            alignmentIterator.setOutputWriter(writer);

            // Iterate through each alignment and write sequence variations to output file:
            alignmentIterator.iterate(basenames);
        }


        finally {

            writer.close();
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
    public static void main(final String[] args) throws JSAPException, IOException {
        new DisplaySequenceVariationsMode().configure(args).execute();
    }

    private static class MyIterateAlignments extends IterateAlignments {
        PrintWriter outputWriter;

        public void setOutputWriter(PrintWriter outputWriter) {
            this.outputWriter = outputWriter;
        }

        public void processAlignmentEntry(Alignments.AlignmentEntry alignmentEntry) {
            outputWriter.print(String.format("%d %s ",

                    alignmentEntry.getQueryIndex(),
                    getReferenceId(alignmentEntry.getTargetIndex())));
            boolean variations = false;
            for (Alignments.SequenceVariation var : alignmentEntry.getSequenceVariationsList()) {
                variations = true;

                outputWriter.print(String.format("%d:%s/%s,",

                        var.getPosition(),
                        var.getFrom(),
                        var.getTo()));
            }
            if (variations) {
                outputWriter.println();

            }
        }
    }
}