/*
 * Copyright (C) 2009 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 * WEILL MEDICAL COLLEGE OF CORNELL UNIVERSITY MAKES NO REPRESENTATIONS
 * ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY PURPOSE. IT IS PROVIDED
 * "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY. THE WEILL MEDICAL COLLEGE
 * OF CORNELL UNIVERSITY SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY
 * THE USERS OF THIS SOFTWARE.
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
public class Alignment2TextMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    public static final String MODE_NAME = "alignment-to-text";
    public static final String MODE_DESCRIPTION = "Converts a compact alignment to plain text.";

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

        PrintWriter writer = outputFile == null ? new PrintWriter(System.out) :
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
                    int queryIndex = alignmentEntry.getQueryIndex();
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
        new Alignment2TextMode().configure(args).execute();
    }
}
