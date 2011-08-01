/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
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
import edu.cornell.med.icb.goby.algorithmic.algorithm.*;
import edu.cornell.med.icb.goby.algorithmic.data.WeightsInfo;
import edu.cornell.med.icb.goby.alignments.*;
import edu.cornell.med.icb.goby.counts.CountWriterHelper;
import edu.cornell.med.icb.goby.counts.CountsArchiveWriter;
import edu.cornell.med.icb.goby.counts.CountsWriter;
import edu.cornell.med.icb.goby.util.Timer;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import java.io.IOException;
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
     * The overridden short mode name.
     */
    private static final String SHORT_MODE_NAME = "altc";

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

    private WeightParameters weightParams;
    private static final Logger LOG = Logger.getLogger(CompactAlignmentToCountsMode.class);
    private boolean verbose;
    private boolean filterAmbiguousReads;

    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    @Override
    public String getShortModeName() {
        return SHORT_MODE_NAME;
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
            basenameSet.add(AlignmentReaderImpl.getBasename(inputFile));
        }

        basenames = basenameSet.toArray(new String[basenameSet.size()]);
        optionalOutputFile = jsapResult.getString("output");

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

        weightParams = CompactAlignmentToAnnotationCountsMode.configureWeights(jsapResult);


        filterAmbiguousReads = jsapResult.getBoolean("filter-ambiguous-reads");
        if (filterAmbiguousReads) {
            System.out.println("Ambiguous reads will not be considered when estimating count statistics.");
        }
        return this;
    }

    /**
     * Run the mode.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        final ProgressLogger progress = new ProgressLogger(LOG);
        progress.priority = Level.INFO;
        progress.expectedUpdates = basenames.length;
        progress.itemsName = "basenames";
        progress.start();
        if (basenames.length > 1) {
            verbose = true;
        }
        for (final String basename : basenames) {
            if (optionalOutputFile == null) {
                outputFile = basename;
            }

            processFullGenomeAlignment(basename);
            progress.info = "last processed: " + basename;
            progress.lightUpdate();

        }
        progress.done();
    }

    private void processFullGenomeAlignment(final String basename) throws IOException {
        final AlignmentReaderFactory factory = filterAmbiguousReads ? new NonAmbiguousAlignmentReaderFactory() :
                new DefaultAlignmentReaderFactory();

        AlignmentReader reader = factory.createReader(basename);

        reader.readHeader();
        if (reader.isSorted() && reader.getNumberOfAlignedReads() > 50000000) {

            // if the alignment is sorted and has more than 50 million aligned reads, we switch to the more scalable
            // count production method. The method scales to any size alignment, but is about 3-4 times slower than
            // directly loading the entire alignment in memory..
            processSortedAlignmentFullGenome(basename, reader);
            return;
        }

        final int numberOfReferences = reader.getNumberOfTargets();

        final DoubleIndexedIdentifier referenceIds = new DoubleIndexedIdentifier(reader.getTargetIdentifiers());
        reader.close();
        //System.out.println(String.format("Alignment contains %d reference sequences", numberOfReferences));
        final ComputeCountInterface[] algs = new ComputeCountInterface[numberOfReferences];
        final CountsArchiveWriter countArchive;
        countArchive = new CountsArchiveWriter(basename, countArchiveModifier);
        // be verbose if only one basename to process.

        countArchive.setVerbose(verbose);
        //  CountsWriter writers[] = new CountsWriter[numberOfReferences];
        final IntSet referencesToProcess = new IntOpenHashSet();
        WeightsInfo weights = null;
        if (weightParams.useWeights) {
            weights = CompactAlignmentToAnnotationCountsMode.loadWeights(basename, weightParams.useWeights,
                    weightParams.weightId);
            if (weights != null) {
                System.err.println("Weights have been provided and loaded and will be used to reweight counts.");
            }
        }

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
                    final ComputeCountInterface algo = new ComputeCount();
                    algs[referenceIndex] = chooseAlgorithm(weightParams, weights, algo);
                } else {
                    algs[referenceIndex] = new ComputeStartCount(focusOnStrand);
                }
                algs[referenceIndex].startPopulating();
            }
        }

        final AlignmentReader referenceReader = new AlignmentReaderImpl(basename);
        referenceReader.readHeader();

        // read the alignment:
        //   System.out.println("Loading the alignment..");
        for (final Alignments.AlignmentEntry alignmentEntry : referenceReader) {
            final int referenceIndex = alignmentEntry.getTargetIndex();
            if (referencesToProcess.contains(referenceIndex)) {
                final int startPosition = alignmentEntry.getPosition();

                final int alignmentLength = alignmentEntry.getTargetAlignedLength();
                for (int i = 0; i < alignmentEntry.getMultiplicity(); ++i) {

                    algs[referenceIndex].populate(startPosition, startPosition + alignmentLength,
                            !alignmentEntry.getMatchingReverseStrand(), alignmentEntry.getQueryIndex());
                }
            }
        }

        reader.close();
        final Timer timer = new Timer();
        timer.start();
        for (final int referenceIndex : referencesToProcess) {
            final String chromosomeName = referenceIds.getId(referenceIndex).toString();

            if (verbose) {
                System.out.println("Writing counts for reference " + chromosomeName);
            }

            algs[referenceIndex].accumulate();
            final CountsWriter countsWriter = countArchive.newCountWriter(referenceIndex, chromosomeName);
            algs[referenceIndex].baseCount(countsWriter);
            countArchive.returnWriter(countsWriter);
            algs[referenceIndex] = null;
            // Runtime.getRuntime().gc();

        }
        countArchive.close();
        timer.stop();
        System.out.println(timer);
    }

    private void processSortedAlignmentFullGenome(final String basename, AlignmentReader reader) throws IOException {

        final CountsArchiveWriter countArchive = new CountsArchiveWriter(basename, countArchiveModifier);
        final IterateForCounts sortedPositionIterator = new IterateForCounts(countArchive);
        sortedPositionIterator.iterate(basename);
        sortedPositionIterator.finishWriter();
    }

    private class IterateForCounts extends IterateSortedAlignmentsListImpl {
        CountsWriter writer;
        CountsArchiveWriter archiveWriter;
        CountWriterHelper helper;

        private IterateForCounts(CountsArchiveWriter archiveWriter) {
            this.archiveWriter = archiveWriter;

        }

        int lastReferenceIndex = -1;

        @Override
        public void processPositions(int referenceIndex, int position, DiscoverVariantPositionData positionBaseInfos) {
            try {
                if (referenceIndex != lastReferenceIndex) {
                    if (writer != null) {
                        finishWriter();
                    }
                    writer = archiveWriter.newCountWriter(referenceIndex, getReferenceId(referenceIndex).toString());
                    helper = new CountWriterHelper(writer);
                    lastReferenceIndex = referenceIndex;
                }
                helper.appendCountAtPosition(positionBaseInfos.size(), position);
            } catch (IOException e) {
                LOG.error("cannot return counts writer to archive", e);

            }
        }

        public void finishWriter() throws IOException {
            helper.close();
            archiveWriter.returnWriter(writer);

            writer = null;
        }

    }

    private ComputeCountInterface chooseAlgorithm(final WeightParameters weightParams, final WeightsInfo weights, ComputeCountInterface algo) {

        if (weightParams.useWeights) {
            if (!weightParams.adjustGcBias) {
                // weights only:
                algo = new ComputeWeightCount(weights);
            } else {
                // use weights to reweight with formula:
                final FormulaWeightCount algo1 = new FormulaWeightCount(weights);

                algo1.setFormulaChoice(FormulaWeightAnnotationCount.FormulaChoice.valueOf(weightParams.formulaChoice));
                algo = algo1;
            }

        }
        return algo;

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

