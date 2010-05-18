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
import edu.cornell.med.icb.goby.stats.DifferentialExpressionAnalysis;
import edu.cornell.med.icb.goby.stats.DifferentialExpressionCalculator;
import edu.cornell.med.icb.goby.stats.DifferentialExpressionResults;
import edu.cornell.med.icb.goby.stats.NormalizationMethod;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.floats.FloatArrayList;
import it.unimi.dsi.fastutil.io.BinIO;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Reads a compact alignment outputs read counts that overlap with transcript annotation segments.
 * <p/>
 * User: nyasha
 * Date: Mar 31, 2010
 * Time: 11:05:47 AM
 */
public class CompactAlignmentToTranscriptCountsMode extends AbstractGobyMode {

    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(CompactAlignmentToAnnotationCountsMode.class);

    /**
     * The mode name.
     */
    private static final String MODE_NAME = "alignment-to-transcript-counts";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Converts alignment to counts for transcripts";
    private String[] inputFiles;
    private String[] basenames;
    /**
     * The output file for transcript counts.
     */
    private String outputFile;
    /**
     * The output file giving summary statistics for differential expression between groups.
     */
    private String statsFilename;
    private boolean doComparison;

    private final DifferentialExpressionCalculator deCalculator = new DifferentialExpressionCalculator();
    private final DifferentialExpressionAnalysis deAnalyzer = new DifferentialExpressionAnalysis();
    /* The set of normalization methods to use for the comparison.
     */
    private ObjectArraySet<NormalizationMethod> normalizationMethods;
    private String weightsFilename;

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
        inputFiles = jsapResult.getStringArray("input");
        final ObjectSet<String> basenameSet = new ObjectOpenHashSet<String>();
        for (final String inputFile : inputFiles) {
            basenameSet.add(AlignmentReader.getBasename(inputFile));
        }
        basenames = basenameSet.toArray(new String[basenameSet.size()]);
        statsFilename = jsapResult.getString("stats");
        outputFile = jsapResult.getString("output");
        weightsFilename = jsapResult.getString("weights");
        final String groupsDefinition = jsapResult.getString("groups");

        deAnalyzer.parseGroupsDefinition(groupsDefinition, deCalculator, inputFiles);

        final String compare = jsapResult.getString("compare");
        if (compare == null) {
            doComparison = false;
        } else {
            doComparison = true;
        }
        if (doComparison) {
            deAnalyzer.parseCompare(compare);
        }
        normalizationMethods = deAnalyzer.parseNormalization(jsapResult);
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
            if (outputFile == null) {
                outputFile = basename;
            }
            outputFile += "-transcript-counts.txt";
            processTranscriptAlignment(basename);
            outputFile = null;
        }


        if (doComparison) {
            PrintWriter statsOutput = null;
            try {
                statsOutput = new PrintWriter(statsFilename);
                final DifferentialExpressionResults results =
                        deAnalyzer.evaluateDifferentialExpressionStatistics(deCalculator, doComparison, normalizationMethods);
                results.write(statsOutput, '\t', deCalculator);
            } finally {
                IOUtils.closeQuietly(statsOutput);
            }
        }
    }

    private void processTranscriptAlignment(final String basename) throws IOException {
        final AlignmentReader reader = new AlignmentReader(basename);
        PrintWriter outputWriter = null;
        try {
            FloatArrayList weights = null;
            if (weightsFilename != null) {
                weights = (FloatArrayList) BinIO.loadObject(weightsFilename);
            } else {
                System.err.println("Weights have not been provided, the 'reweightedCounts' column will contain the same value as count.");
            }
            outputWriter = new PrintWriter(new FileWriter(outputFile));

            // outputWriter.write("# One line per reference id. Count indicates the number of times a query \n" +
            //         "# partially overlaps a target, given the various quality filters used to create the alignment.\n");
            outputWriter.write("sampleId\treferenceId\tcount\tlog10(count+1)\tcumulativeBasesAligned\treweightedCounts\n");

            reader.readHeader();

            final int numberOfReferences = reader.getNumberOfTargets();
            final int[] numberOfReadsPerReference = new int[numberOfReferences];
            final double[] reweightedCounts = new double[numberOfReferences];
            final int[] cumulativeBasesPerReference = new int[numberOfReferences];


            System.out.printf("Scanning alignment %s%n", basename);
            for (final Alignments.AlignmentEntry alignmentEntry : reader) {
                final int referenceIndex = alignmentEntry.getTargetIndex();
                ++numberOfReadsPerReference[referenceIndex];

                reweightedCounts[referenceIndex] += (weights != null ? weights.getFloat(alignmentEntry.getQueryIndex()) : 1);

                cumulativeBasesPerReference[referenceIndex] +=
                        Math.min(alignmentEntry.getQueryAlignedLength(),
                                alignmentEntry.getTargetAlignedLength());
            }
            final IndexedIdentifier targetIds = reader.getTargetIdentifiers();

            final DoubleIndexedIdentifier targetIdBackward = new DoubleIndexedIdentifier(targetIds);

            final String sampleId = FilenameUtils.getBaseName(basename);
            deCalculator.reserve(numberOfReferences, inputFiles.length);
            int numAlignedReadsInSample = 0;
            // define elements that will be tested for differential expression:
            for (int referenceIndex = 0; referenceIndex < numberOfReferences; ++referenceIndex) {

                final String transcriptId = targetIdBackward.getId(referenceIndex).toString();
                final int index = deCalculator.defineElement(transcriptId, DifferentialExpressionCalculator.ElementType.TRANSCRIPT);

                deCalculator.defineElementLength(index, reader.getTargetLength(referenceIndex));
            }

            // observe elements:
            for (int referenceIndex = 0; referenceIndex < numberOfReferences; ++referenceIndex) {

                outputWriter.printf("%s\t%s\t%d\t%g\t%d\t%g%n", basename,
                        targetIdBackward.getId(referenceIndex),
                        numberOfReadsPerReference[referenceIndex],
                        Math.log10(numberOfReadsPerReference[referenceIndex] + 1),
                        cumulativeBasesPerReference[referenceIndex],
                        reweightedCounts[referenceIndex]);

                final String transcriptId = targetIdBackward.getId(referenceIndex).toString();

                deCalculator.observe(sampleId, transcriptId, numberOfReadsPerReference[referenceIndex]);
                numAlignedReadsInSample += numberOfReadsPerReference[referenceIndex];
            }
            deCalculator.setNumAlignedInSample(sampleId, numAlignedReadsInSample);
            outputWriter.flush();

        } catch (ClassNotFoundException e) {
            LOG.error("Cannot load serialized weights from filename " + weightsFilename);
            System.exit(1);
        } finally {
            IOUtils.closeQuietly(outputWriter);
            reader.close();
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
        new CompactAlignmentToTranscriptCountsMode().configure(args).execute();
    }

}
