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

import cern.colt.Timer;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.algorithmic.algorithm.ComputeCount;
import edu.cornell.med.icb.goby.algorithmic.algorithm.ComputeStartCount;
import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.counts.CountsArchiveWriter;
import edu.cornell.med.icb.goby.counts.CountsWriter;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.apache.commons.io.IOUtils;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

/**
 * Converts a compact alignment to a compressed count archive.
 *
 * @author Fabien Campagne
 */
public class CompactAlignmentToCountsMode extends AbstractGobyMode {

    /**
     * The mode name.
     */
    private static final String MODE_NAME = "alignment-to-counts";
    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION =
            "Converts a compact alignment to counts.";

    /**
     * Default counts archive extension.
     */
    public static final String COUNT_ARCHIVE_MODIFIER_DEFAULT = "counts";

    /**
     * The output file.
     */
    private String outputFile;

    private boolean filterByReferenceNames;
    private ObjectSet<String> includeReferenceNames = new ObjectOpenHashSet<String>();
    private boolean fullGenomeAlignment;
    private String[] basenames;
    private String optionalOutputFile;
    private boolean accumulatePeakHistogram;
    private int focusOnStrand;

    private String countArchiveModifier = COUNT_ARCHIVE_MODIFIER_DEFAULT;

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

        final String[] inputFiles = jsapResult.getStringArray("input");
        final ObjectSet<String> basenameSet = new ObjectOpenHashSet<String>();

        for (final String inputFile : inputFiles) {
            basenameSet.add(AlignmentReader.getBasename(inputFile));
        }

        basenames = basenameSet.toArray(new String[basenameSet.size()]);

        optionalOutputFile = jsapResult.getString("output");
        fullGenomeAlignment = jsapResult.getBoolean("full-genome");

        final String includeReferenceNameCommas = jsapResult.getString("include-reference-names");
        if (includeReferenceNameCommas != null) {
            includeReferenceNames = new ObjectOpenHashSet<String>();
            includeReferenceNames.addAll(Arrays.asList(includeReferenceNameCommas.split("[,]")));
            System.out.println("Will write counts for the following sequences:");
            for (final String name : includeReferenceNames) {
                System.out.println(name);
            }
            filterByReferenceNames = true;
        }

        accumulatePeakHistogram = !jsapResult.getBoolean("start-only", false);
        if (!accumulatePeakHistogram) {
            final String strandChoice = jsapResult.getString("strand-choice");
            if ("POSITIVE".equalsIgnoreCase(strandChoice) || "FORWARD".equalsIgnoreCase(strandChoice)) {
                focusOnStrand = ComputeStartCount.POSITIVE_STRAND_ONLY;
                countArchiveModifier = "forward-strand-starts";
            } else if ("NEGATIVE".equalsIgnoreCase(strandChoice) || "REVERSE".equalsIgnoreCase(strandChoice)) {
                focusOnStrand = ComputeStartCount.REVERSE_STRAND_ONLY;
                countArchiveModifier = "reverse-strand-starts";
            } else if ("BOTH".equalsIgnoreCase(strandChoice) || "EITHER".equalsIgnoreCase(strandChoice)) {
                focusOnStrand = ComputeStartCount.BOTH_STRAND;
                countArchiveModifier = "either-strand-starts";
            } else {
                System.err.println("strand choice must be one of {forward, positive, reverse, negative, both, either}.");
                System.exit(2);
            }
        }
        return this;
    }

    /**
     * Run the map2text mode.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        for (final String basename : basenames) {
            if (optionalOutputFile == null) {
                outputFile = basename;
            }

            if (fullGenomeAlignment) {
                processFullGenomeAlignment(basename);
            } else {
                outputFile += "-transcript-counts.txt";
                processTranscriptAlignment(basename);
            }
        }
    }

    private void processTranscriptAlignment(final String basename) throws IOException {
        final AlignmentReader reader = new AlignmentReader(basename);
        PrintWriter outputWriter = null;
        try {
            outputWriter = new PrintWriter(new FileWriter(outputFile));
            // outputWriter.write("# One line per reference id. Count indicates the number of times a query \n" +
            //         "# partially overlaps a target, given the various quality filters used to create the alignment.\n");
            outputWriter.write("sampleId\treferenceId\tcount\tlog10(count+1)\tcumulativeBasesAligned\n");
            reader.readHeader();

            final int numberOfReferences = reader.getNumberOfTargets();
            final int[] numberOfReadsPerReference = new int[numberOfReferences];
            final int[] cumulativeBasesPerReference = new int[numberOfReferences];

            System.out.printf("Scanning alignment %s..%n", basename);
            for (final Alignments.AlignmentEntry alignmentEntry : reader) {
                ++numberOfReadsPerReference[alignmentEntry.getTargetIndex()];

                cumulativeBasesPerReference[alignmentEntry.getTargetIndex()] +=
                        Math.min(alignmentEntry.getQueryAlignedLength(),
                                alignmentEntry.getTargetAlignedLength());
            }
            final IndexedIdentifier targetIds = reader.getTargetIdentifiers();
            final DoubleIndexedIdentifier targetIdBackward = new DoubleIndexedIdentifier(targetIds);

            for (int referenceIndex = 0; referenceIndex < numberOfReferences; ++referenceIndex) {
                outputWriter.printf("%s\t%s\t%d\t%g\t%d%n", basename,
                        targetIdBackward.getId(referenceIndex),
                        numberOfReadsPerReference[referenceIndex],
                        Math.log10(numberOfReadsPerReference[referenceIndex] + 1),
                        cumulativeBasesPerReference[referenceIndex]);
                referenceIndex++;
            }
            outputWriter.flush();
        } finally {
            IOUtils.closeQuietly(outputWriter);
            reader.close();
        }
    }

    private void processFullGenomeAlignment(final String basename) throws IOException {
        final AlignmentReader reader = new AlignmentReader(basename);
        reader.readHeader();
        final int numberOfReferences = reader.getNumberOfTargets();

        final DoubleIndexedIdentifier referenceIds = new DoubleIndexedIdentifier(reader.getTargetIdentifiers());
        reader.close();
        System.out.println(String.format("Alignment contains %d reference sequences", numberOfReferences));
        final ComputeCount[] algs = new ComputeCount[numberOfReferences];
        final CountsArchiveWriter countArchive;
        countArchive = new CountsArchiveWriter(basename, countArchiveModifier);
        //  CountsWriter writers[] = new CountsWriter[numberOfReferences];
        final IntSet referencesToProcess = new IntOpenHashSet();

        // create count writers, one for each reference sequence in the alignment:
        for (int referenceIndex = 0; referenceIndex < numberOfReferences; referenceIndex++) {


            final String referenceName = referenceIds.getId(referenceIndex).toString();
            if (filterByReferenceNames) {
                if (includeReferenceNames.contains(referenceName)) {
                    // subset of reference names selected by the command line:
                    referencesToProcess.add(referenceIndex);
                }
            } else {
                // process each sequence:
                referencesToProcess.add(referenceIndex);
            }

            if (referencesToProcess.contains(referenceIndex)) {
                if (accumulatePeakHistogram) {
                    algs[referenceIndex] = new ComputeCount();
                } else {
                    algs[referenceIndex] = new ComputeStartCount(focusOnStrand);
                }
                algs[referenceIndex].startPopulating();
            }
        }

        final AlignmentReader referenceReader = new AlignmentReader(basename);
        referenceReader.readHeader();

        // read the alignment:
        System.out.println("Loading the alignment..");
        for (final Alignments.AlignmentEntry alignmentEntry : referenceReader) {
            final int referenceIndex = alignmentEntry.getTargetIndex();
            if (referencesToProcess.contains(referenceIndex)) {
                final int startPosition = alignmentEntry.getPosition();

                final int alignmentLength = alignmentEntry.getQueryAlignedLength();
                for (int i = 0; i < alignmentEntry.getMultiplicity(); ++i) {

                    algs[referenceIndex].populate(startPosition, startPosition + alignmentLength,
                            !alignmentEntry.getMatchingReverseStrand());
                }
            }
        }

        reader.close();
        final Timer timer = new Timer();
        timer.start();
        for (final int referenceIndex : referencesToProcess) {
            final String chromosomeName = referenceIds.getId(referenceIndex).toString();

            System.out.println("Writing counts for reference " + chromosomeName);

            algs[referenceIndex].accumulate();
            final CountsWriter countsWriter = countArchive.newCountWriter(referenceIndex, chromosomeName);
            algs[referenceIndex].baseCount(countsWriter);
            countArchive.returnWriter(countsWriter);
            algs[referenceIndex] = null;
            // Runtime.getRuntime().gc();

        }
        countArchive.close();
        timer.stop();
        System.out.println(String.format("time spent  %d ms %g secs %g mins",
                timer.millis(), timer.seconds(),
                timer.minutes()));
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
        new CompactAlignmentToCountsMode().configure(args).execute();
    }
}
