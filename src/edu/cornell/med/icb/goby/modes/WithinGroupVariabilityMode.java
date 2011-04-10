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
import edu.cornell.med.icb.goby.algorithmic.algorithm.AnnotationCount;
import edu.cornell.med.icb.goby.algorithmic.data.Annotation;
import edu.cornell.med.icb.goby.algorithmic.data.Segment;
import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
import edu.cornell.med.icb.goby.stats.AverageFisherRCalculator;
import edu.cornell.med.icb.goby.stats.BullardUpperQuartileNormalization;
import edu.cornell.med.icb.goby.stats.DifferentialExpressionCalculator;
import edu.cornell.med.icb.goby.stats.DifferentialExpressionResults;
import edu.cornell.med.icb.goby.stats.FisherExactRCalculator;
import edu.cornell.med.icb.goby.stats.NormalizationMethod;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.Object2ObjectMap;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

/**
 * Reads a compact alignment and determines Fisher P-values comparing pairs of samples within a group.
 * <p/>
 * @author Fabien Campagne
 * Date: May 4, 2010
 * Time: 11:05:47 AM
 */
public class WithinGroupVariabilityMode extends AbstractGobyMode {

    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(WithinGroupVariabilityMode.class);

    /**
     * The mode name.
     */
    private static final String MODE_NAME = "within-group-variability";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Estimate variability within groups for each element of differential expression.";
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

    /* The set of normalization methods to use for the comparison.
     */
    private ObjectArraySet<NormalizationMethod> normalizationMethods;
    private ObjectOpenHashSet<String> includeAnnotationTypes;
    private final boolean transcripts = false;

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
            basenameSet.add(AlignmentReaderImpl.getBasename(inputFile));
        }
        basenames = basenameSet.toArray(new String[basenameSet.size()]);
        statsFilename = jsapResult.getString("stats");
        outputFile = jsapResult.getString("output");
        final String includeReferenceNameCommas = jsapResult.getString("include-reference-names");
        includeReferenceNames = new ObjectOpenHashSet<String>();
        if (includeReferenceNameCommas != null) {
            includeReferenceNames.addAll(Arrays.asList(includeReferenceNameCommas.split("[,]")));
            System.out.println("Will write counts for the following sequences:");
            for (final String name : includeReferenceNames) {
                System.out.println(name);
            }
            filterByReferenceNames = true;
        }

        parseAnnotations(jsapResult);

        return this;
    }

    private void parseAnnotations(final JSAPResult jsapResult) {
        annotationFile = jsapResult.getString("annotation");
        final String includeAnnotationTypeCommas = jsapResult.getString("include-annotation-types");
        includeAnnotationTypes = new ObjectOpenHashSet<String>();
        if (includeAnnotationTypeCommas != null) {
            includeAnnotationTypes.addAll(Arrays.asList(includeAnnotationTypeCommas.split("[,]")));
            for (final String name : includeAnnotationTypes) {
                if (name.equals("gene") || name.equals("other") || name.equals("exon")) {
                    continue;
                } else {
                    System.out.println("Please enter a valid annotation type. "
                            + "Valid annotation types include gene, exon, and or other.");
                    System.exit(1);
                }
                System.out.println("Will write counts for the following annotation types:");
                System.out.println(name);
            }
        }
    }

    /**
     * The annotation file.
     */
    private String annotationFile;
    private final DifferentialExpressionCalculator deCalculator = new DifferentialExpressionCalculator();

    /**
     * Run the map2text mode.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {


        // each sample belongs to its own group:
        for (final String sample : basenames) {
            final String basename = AlignmentReaderImpl.getBasename(FilenameUtils.getBaseName(sample));
            deCalculator.defineGroup(basename);
            deCalculator.associateSampleToGroup(basename, basename);


        }
        System.out.println("Reading annotations from " + annotationFile);
        final Object2ObjectMap<String, ObjectList<Annotation>> allAnnots = CompactAlignmentToAnnotationCountsMode.readAnnotations(annotationFile);
        for (final String basename : basenames) {
            if (outputFile == null) {
                outputFile = basename;
            }

            if (transcripts) {
                processTranscriptAlignment(basename);
            } else {
                processOneBasename(allAnnots, basename);
            }
            outputFile = null;
        }


        PrintWriter statsOutput = null;
        try {
            statsOutput = new PrintWriter(statsFilename);
            DifferentialExpressionResults results = null;

            results = null;

            final NormalizationMethod method = new BullardUpperQuartileNormalization();


            // evaluate differences between samples:
            int i = 0;
            int j = 0;
            final boolean[][] done = new boolean[basenames.length][basenames.length];

            for (i=0;i<basenames.length; i++) {
                for (j=0; j<basenames.length; j++) {
                    if (i != j && !done[i][j]) {
                        final String sample1 = basenames[i];
                        final String sample2 = basenames[j];
                        // consider each pair of samples only once:
                        final String basename1 = AlignmentReaderImpl.getBasename(FilenameUtils.getBaseName(sample1));
                        final String basename2 = AlignmentReaderImpl.getBasename(FilenameUtils.getBaseName(sample2));
                        method.normalize(deCalculator, basename1, basename2);
                        results = deCalculator.compare(results, method, new FisherExactRCalculator(), basename1, basename2);
                        done[i][j] = true;
                        done[j][i] = true;
                    }

                }

            }
            results = deCalculator.compare(results, method, new AverageFisherRCalculator());
            if (results != null) {
                results.write(statsOutput, '\t', deCalculator);
            } else {
                System.out.println("No results were produced. Make sure R is configured properly.");
            }

            IOUtils.closeQuietly(statsOutput);
        } finally {
            IOUtils.closeQuietly(statsOutput);

        }

    }

    private boolean filterByReferenceNames;
    private ObjectSet<String> includeReferenceNames;

    private void processOneBasename(final Object2ObjectMap<String, ObjectList<Annotation>> allAnnots, final String inputBasename) throws IOException {
        final AlignmentReaderImpl reader = new AlignmentReaderImpl(inputBasename);
        reader.readHeader();
        final int numberOfReferences = reader.getNumberOfTargets();

        final DoubleIndexedIdentifier referenceIds = new DoubleIndexedIdentifier(reader.getTargetIdentifiers());
        reader.close();
        System.out.println(String.format("Alignment contains %d reference sequences", numberOfReferences));
        final AnnotationCount[] algs = new AnnotationCount[numberOfReferences];
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
                algs[referenceIndex] = new AnnotationCount();
                algs[referenceIndex].getBaseCounter().startPopulating();
            }
        }


        final AlignmentReader referenceReader = new AlignmentReaderImpl(inputBasename);
        referenceReader.readHeader();

        // read the alignment:
        System.out.println("Loading alignment " + inputBasename + "..");
        int numAlignedReadsInSample = 0;
        for (final Alignments.AlignmentEntry alignmentEntry : referenceReader) {
            final int referenceIndex = alignmentEntry.getTargetIndex();
            if (referencesToProcess.contains(referenceIndex)) {
                final int startPosition = alignmentEntry.getPosition();

                final int alignmentLength = alignmentEntry.getQueryAlignedLength();
                //shifted the ends populating by 1
                for (int i = 0; i < alignmentEntry.getMultiplicity(); ++i) {
                    algs[referenceIndex].populate(startPosition, startPosition + alignmentLength);
                    ++numAlignedReadsInSample;
                }
            }
        }

        reader.close();
        final String sampleId = FilenameUtils.getBaseName(inputBasename);

        deCalculator.setNumAlignedInSample(sampleId, numAlignedReadsInSample);
        observeCounts(allAnnots, inputBasename, referenceIds, algs, referencesToProcess);

    }

    private void observeCounts(final Object2ObjectMap<String, ObjectList<Annotation>> allAnnots,
                               final String inputBasename,
                               final DoubleIndexedIdentifier referenceIds, final AnnotationCount[] algs,
                               final IntSet referencesToProcess) throws IOException {

        // collect all element ids:

        int numberOfElements = 0;
        int numberOfGenes = 0;
        int numberOfExons = 0;
        int numberOfIntrons = 0;
        for (final int referenceIndex : referencesToProcess) {
            final String chromosomeName = referenceIds.getId(referenceIndex).toString();

            if (!allAnnots.containsKey(chromosomeName)) {
                continue;
            }
            final ObjectList<Annotation> annots = allAnnots.get(chromosomeName);

            for (final Annotation annot : annots) {
                final String geneID = annot.getId();
                final int numExons = annot.getSegments().size();
                final int numberIntrons = numExons - 1;

                if (includeAnnotationTypes.contains("gene")) {
                    final int index = deCalculator.defineElement(geneID, DifferentialExpressionCalculator.ElementType.GENE);
                    deCalculator.defineElementLength(index, annot.getLength());
                    numberOfGenes++;
                    numberOfElements++;
                }

                if (includeAnnotationTypes.contains("exon")) {

                    for (int i = 0; i < numExons; i++) {
                        final Segment exonSegment = annot.getSegments().get(i);
                        final String exonID = exonSegment.getId();
                        final int index = deCalculator.defineElement(exonID, DifferentialExpressionCalculator.ElementType.EXON);
                        deCalculator.defineElementLength(index, annot.getLength());
                        numberOfExons++;
                        numberOfElements++;
                    }
                }

                if (includeAnnotationTypes.contains("other")) {
                    for (int i = 0; i < numExons; i++) {
                        if (i < numberIntrons) {
                            final Segment segment = annot.getSegments().get(i);
                            final int intronStart = segment.getEnd() + 1;
                            final Segment intronSegment = annot.getSegments().get(i + 1);
                            final int intronEnd = intronSegment.getStart() - 1;
                            final int intronLength = intronEnd - intronStart + 1;
                            final String intronID = segment.getId() + "-" + intronSegment.getId();
                            final int index = deCalculator.defineElement(intronID, DifferentialExpressionCalculator.ElementType.OTHER);
                            deCalculator.defineElementLength(index, intronLength);
                            numberOfIntrons++;
                            numberOfElements++;
                        }
                    }
                }
            }
        }
        LOG.info(String.format("%d Genes %d exons %d other total %d ", numberOfGenes, numberOfExons, numberOfIntrons, numberOfElements));
        deCalculator.reserve(numberOfElements, basenames.length);


        int numberOfAnottationCountsWritten = 0;
        for (final int referenceIndex : referencesToProcess) {
            final String chromosomeName = referenceIds.getId(referenceIndex).toString();
            System.out.println("Observing counts for reference " + chromosomeName);

            if (!allAnnots.containsKey(chromosomeName)) {
                continue;
            }

            final ObjectList<Annotation> annots = allAnnots.get(chromosomeName);
            algs[referenceIndex].sortReads();
            algs[referenceIndex].getBaseCounter().accumulate();
            algs[referenceIndex].getBaseCounter().baseCount();
            if (doComparison) {
                for (final Annotation annot : annots) {
                    final String geneID = annot.getId();
                    deCalculator.defineElement(geneID);
                }
            }
            final String basename = FilenameUtils.getBaseName(inputBasename);
            final String sampleId = basename;
            for (final Annotation annot : annots) {
                final String geneID = annot.getId();

                if (includeAnnotationTypes.contains("gene")) {
                    final int geneStart = annot.getStart();
                    final int geneEnd = annot.getEnd();
                    final int geneLength = geneEnd - geneStart + 1;
                    final float geneDepth = algs[referenceIndex].averageReadsPerPosition(geneStart, geneEnd);
                    final double geneOverlapReads = algs[referenceIndex].countReadsPartiallyOverlappingWithInterval(geneStart, geneEnd);
                    final double geneInsideReads = algs[referenceIndex].countReadsStriclyWithinInterval(geneStart, geneEnd);
                    final double geneExpression = algs[referenceIndex].geneExpressionCount(annot);
                    final int numExons = annot.getSegments().size();


                    final double geneRPKM = deCalculator.calculateNormalized(geneOverlapReads, annot.getLength(), deCalculator.getNumAlignedInSample(sampleId));

                    numberOfAnottationCountsWritten++;

                    deCalculator.observe(basename, geneID, geneExpression);

                }
                final int numberExons = annot.getSegments().size();
                final int numberIntrons = numberExons - 1;
                // skip unnecessary computation if we don't need exon or intron info:
                if (includeAnnotationTypes.contains("exon") || includeAnnotationTypes.contains("other")) {
                    for (int i = 0; i < numberExons; i++) {
                        final Segment segment = annot.getSegments().get(i);
                        final int exonStart = segment.getStart();
                        final int exonEnd = segment.getEnd();

                        final String exonStrand = segment.getStrand();
                        final int exonLength = segment.getLength();
                        final String exonID = segment.getId();
                        final float exonDepth = algs[referenceIndex].averageReadsPerPosition(exonStart, exonEnd);
                        final double exonOverlapReads = algs[referenceIndex].countReadsPartiallyOverlappingWithInterval(exonStart, exonEnd);
                        final double exonInsideReads = algs[referenceIndex].countReadsStriclyWithinInterval(exonStart, exonEnd);
                        final double exonRPKM = deCalculator.calculateNormalized(exonOverlapReads, segment.getLength(), deCalculator.getNumAlignedInSample(sampleId));
                        if (includeAnnotationTypes.contains("exon")) {

                            numberOfAnottationCountsWritten++;
                            if (includeAnnotationTypes.contains("exon")) {
                                deCalculator.observe(basename, exonID, exonOverlapReads);
                            }
                        }
                        if (i < numberIntrons) {
                            final int intronStart = segment.getEnd() + 1;
                            final Segment intronSegment = annot.getSegments().get(i + 1);
                            final int intronEnd = intronSegment.getStart() - 1;
                            final int intronLength = intronEnd - intronStart + 1;
                            final String intronID = segment.getId() + "-" + intronSegment.getId();
                            final float intronDepth = algs[referenceIndex].averageReadsPerPosition(intronStart, intronEnd);
                            final double intronOverlapReads = algs[referenceIndex].countReadsPartiallyOverlappingWithInterval(intronStart, intronEnd);
                            final double intronInsideReads = algs[referenceIndex].countReadsStriclyWithinInterval(intronStart, intronEnd);
                            final double intronRPKM = deCalculator.calculateNormalized(intronOverlapReads, intronSegment.getLength(), deCalculator.getNumAlignedInSample(sampleId));
                            if (intronLength > 0) {
                                if (includeAnnotationTypes.contains("intron")) {

                                    numberOfAnottationCountsWritten++;
                                    if (includeAnnotationTypes.contains("other")) {
                                        deCalculator.observe(basename, intronID, intronOverlapReads);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            algs[referenceIndex] = null;
        }

        LOG.info("Wrote " + numberOfAnottationCountsWritten + " entries");
        if (numberOfAnottationCountsWritten == 0) {
            LOG.warn("No entries were written.  This may be due to the fact that names "
                    + "in the reference dataset used do not match those in the annotation file.  "
                    + "For example, ENSEMBL names chromosomes \"1\",\"2\",\"3\" whereas UCSC "
                    + "names the same chromosomes \"chr1\",\"chr2\",\"chr3\". In these "
                    + "cases you will need to adjust the names in the annotation file being used "
                    + "so they match the names used in the reference dataset.");
        }
    }

    private void processTranscriptAlignment
            (
                    final String basename) throws IOException {
        final AlignmentReaderImpl reader = new AlignmentReaderImpl(basename);
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

            System.out.printf("Scanning alignment %s%n", basename);
            for (final Alignments.AlignmentEntry alignmentEntry : reader) {
                ++numberOfReadsPerReference[alignmentEntry.getTargetIndex()];

                cumulativeBasesPerReference[alignmentEntry.getTargetIndex()] +=
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

                outputWriter.printf("%s\t%s\t%d\t%g\t%d%n", basename,
                        targetIdBackward.getId(referenceIndex),
                        numberOfReadsPerReference[referenceIndex],
                        Math.log10(numberOfReadsPerReference[referenceIndex] + 1),
                        cumulativeBasesPerReference[referenceIndex]);

                final String transcriptId = targetIdBackward.getId(referenceIndex).toString();

                deCalculator.observe(sampleId, transcriptId, numberOfReadsPerReference[referenceIndex]);
                numAlignedReadsInSample += numberOfReadsPerReference[referenceIndex];
            }
            deCalculator.setNumAlignedInSample(sampleId, numAlignedReadsInSample);
            outputWriter.flush();

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
    public static void main
            (
                    final String[] args) throws JSAPException, IOException {
        new WithinGroupVariabilityMode().configure(args).execute();
    }

}
