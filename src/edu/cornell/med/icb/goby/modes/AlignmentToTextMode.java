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
import org.apache.commons.io.IOUtils;

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
     * The output file.
     */
    private String outputFile;

    /**
     * The basename of the compact alignment.
     */
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
     * @throws IOException error parsing
     * @throws JSAPException error parsing
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args)
            throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        final String inputFile = jsapResult.getString("input");
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
        AlignmentReader reader = null;
        PrintWriter writer = null;
        try {
            reader = new AlignmentReader(basename);
            reader.readHeader();
            final int numberOfReferences = reader.getNumberOfTargets();
            final DoubleIndexedIdentifier referenceIds =
                    new DoubleIndexedIdentifier(reader.getTargetIdentifiers());
            final int[] referenceLengths = reader.getTargetLength();
            System.out.println("Alignment contains " + numberOfReferences + " reference sequences");

            // create count writers, one for each reference sequence in the alignment:
            final IntSet referencesToProcess = new IntOpenHashSet();
            for (int referenceIndex = 0; referenceIndex < numberOfReferences; referenceIndex++) {
                referencesToProcess.add(referenceIndex);
            }

            // read the alignment:
            System.out.println("Converting the alignment..");

            writer = outputFile == null ? new PrintWriter(System.out) : new PrintWriter(outputFile);

            final boolean hasReadIds = reader.getQueryIdentifiers().size() > 0;
            final DoubleIndexedIdentifier readIds =
                    new DoubleIndexedIdentifier(reader.getQueryIdentifiers());

            for (final Alignments.AlignmentEntry alignmentEntry : reader) {
                final int referenceIndex = alignmentEntry.getTargetIndex();
                final String referenceName = referenceIds.getId(referenceIndex).toString();

                if (referencesToProcess.contains(referenceIndex)) {
                    final int startPosition = alignmentEntry.getPosition();
                    final int alignmentLength = alignmentEntry.getQueryAlignedLength();
                    for (int i = 0; i < alignmentEntry.getMultiplicity(); ++i) {
                        final int queryIndex = alignmentEntry.getQueryIndex();

                        // TODO - reference length?
                        writer.write(String.format("%s\t%s\t%d\t%d\t%g\t%d\t%d\t%b%n",
                                hasReadIds ? readIds.getId(queryIndex) : queryIndex,
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
        } finally {
            IOUtils.closeQuietly(writer);
            if (reader != null) {
                reader.close();
            }
        }
    }

    /**
     * Main method.
     *
     * @param args command line args.
     * @throws JSAPException error parsing
     * @throws IOException error parsing or executing.
     */
    public static void main(final String[] args) throws JSAPException, IOException {
        new AlignmentToTextMode().configure(args).execute();
    }
}
