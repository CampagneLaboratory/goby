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

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.apache.commons.io.FilenameUtils;

/**
 * Display the sequence variations found in alignments.
 *
 * @author Fabien Campagne
 */
// TODO report the position of the variant relative to the start of the read. This is useful to
// software that need to assess the signigicance of a variation, because sequencing errors occur
// with different probabilities as sequencing progresses through a read.
// We will need to store this position in the SequenceVariations data structure or calculate it
// here, using read length if this is possible. Need to assess how the calculation should be done,
// it can be tricky when considering combinations of insertions, deletions and mutations.
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

    enum OutputFormat {
        CONCISE,
        TSV,
        TAB_DELIMITED,
    }

    private OutputFormat outputFormat;

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
        outputFormat = OutputFormat.valueOf(jsapResult.getString("format").toUpperCase());

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
        switch (outputFormat) {
            case CONCISE:

                break;
            case TAB_DELIMITED:
            case TSV:
                writer.println("basename target-id query-index position-on-reference var-from var-to");
                break;
        }

        try {
            alignmentIterator.setOutputWriter(writer, outputFormat);

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
        private OutputFormat outputFormat;

        public void setOutputWriter(PrintWriter outputWriter, OutputFormat outputFormat) {
            this.outputWriter = outputWriter;
            this.outputFormat = outputFormat;
        }

        public void processAlignmentEntry(AlignmentReader alignmentReader, Alignments.AlignmentEntry alignmentEntry) {
            String basename = alignmentReader.basename();
            // remove the path:
            basename = FilenameUtils.getBaseName(basename);
            switch (outputFormat) {

                case CONCISE: {
                    if (alignmentEntry.getSequenceVariationsCount() > 0) {
                        outputWriter.print(String.format("%d %s ",

                                alignmentEntry.getQueryIndex(),
                                getReferenceId(alignmentEntry.getTargetIndex())));
                        boolean variations = false;
                        for (Alignments.SequenceVariation var : alignmentEntry.getSequenceVariationsList()) {
                            variations = true;
                            // convert variation position to position on the reference:
                            final int positionOnReference = alignmentEntry.getPosition() + var.getPosition();

                            outputWriter.print(String.format("%d:%s/%s,",

                                    positionOnReference,
                                    var.getFrom(),
                                    var.getTo()));
                        }
                        if (variations) {
                            outputWriter.println();

                        }
                    }
                }
                break;
                case TSV:
                case TAB_DELIMITED: {
                    boolean variations = false;

                    for (Alignments.SequenceVariation var : alignmentEntry.getSequenceVariationsList()) {
                        variations = true;
                        // convert variation position to position on the reference:
                        final int positionOnReference = alignmentEntry.getPosition() + var.getPosition();
                        outputWriter.println(String.format("%s\t%s\t%d\t%d\t%s\t%s",
                                basename,
                                getReferenceId(alignmentEntry.getTargetIndex()),
                                alignmentEntry.getQueryIndex(),
                                positionOnReference,
                                var.getFrom(),
                                var.getTo()));
                    }

                    if (variations) {
                        outputWriter.println();

                    }
                }
                break;
            }


        }
    }
}