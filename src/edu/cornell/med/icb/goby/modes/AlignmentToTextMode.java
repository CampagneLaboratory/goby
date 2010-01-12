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
import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Converts a compact alignment to plain text.
 *
 * @author Fabien Campagne
 */
public class AlignmentToTextMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "alignment-to-text";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Converts a compact alignment to plain text.";

    /**
     * The input file.
     */
    private String inputFile;

    /**
     * The output file.
     */
    private String outputFile;
    private String basename;

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

        inputFile = jsapResult.getString("input");
        basename = AlignmentReader.getBasename(inputFile);
        outputFile = jsapResult.getString("output");
        return this;
    }

    /**
     * Run the map2text mode.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        final AlignmentReader reader = new AlignmentReader(basename);
        reader.readHeader();
        final int numberOfReferences = reader.getNumberOfTargets();

        final DoubleIndexedIdentifier referenceIds = new DoubleIndexedIdentifier(reader.getTargetIdentifiers());
        reader.close();
        System.out.println(String.format("Alignment contains %d reference sequences", numberOfReferences));

        //  CountsWriter writers[] = new CountsWriter[numberOfReferences];
        final IntSet referencesToProcess = new IntOpenHashSet();

        // create count writers, one for each reference sequence in the alignment:
        for (int referenceIndex = 0; referenceIndex < numberOfReferences; referenceIndex++) {
            referencesToProcess.add(referenceIndex);
        }

        final AlignmentReader referenceReader = new AlignmentReader(inputFile);
        referenceReader.readHeader();

        final PrintWriter writer = outputFile == null ? new PrintWriter(System.out) :
                new PrintWriter(new FileWriter(outputFile));

        // read the alignment:
        System.out.println("Converting the alignment..");
        for (final Alignments.AlignmentEntry alignmentEntry : referenceReader) {
            final int referenceIndex = alignmentEntry.getTargetIndex();
            final String referenceName = referenceIds.getId(referenceIndex).toString();
            //           System.out.println(referenceName);
            if (referencesToProcess.contains(referenceIndex)) {
                final int startPosition = alignmentEntry.getPosition();
                final int alignmentLength = alignmentEntry.getQueryAlignedLength();
                for (int i = 0; i < alignmentEntry.getMultiplicity(); ++i) {
//                    System.out.println(startPosition+"   "+ startPosition + alignmentLength+
//                            alignmentEntry.getMatchingReverseStrand());
                    final int queryIndex = alignmentEntry.getQueryIndex();
                    writer.write(String.format("%d\t%s\t%d\t%d\t%g\t%d\t%d\t%b%n",
                            queryIndex,
                            referenceName,
                            alignmentEntry.getNumberOfIndels(),
                            alignmentEntry.getNumberOfMismatches(),
                            alignmentEntry.getScore(),
                            startPosition,
                            alignmentLength,
                            alignmentEntry.getMatchingReverseStrand()));
                }
            }
        }
        writer.close();
        reader.close();
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
        new AlignmentToTextMode().configure(args).execute();
    }
}
