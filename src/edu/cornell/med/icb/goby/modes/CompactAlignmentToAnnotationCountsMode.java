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
import edu.cornell.med.icb.goby.algorithmic.algorithm.AnnotationCountInterface;
import edu.cornell.med.icb.goby.algorithmic.algorithm.AnnotationCountIterateAlignments;
import edu.cornell.med.icb.goby.algorithmic.data.Annotation;
import edu.cornell.med.icb.goby.algorithmic.data.GroupComparison;
import edu.cornell.med.icb.goby.algorithmic.data.Segment;
import edu.cornell.med.icb.goby.algorithmic.data.WeightsInfo;
import edu.cornell.med.icb.goby.algorithmic.data.xml.AnnotationLength;
import edu.cornell.med.icb.goby.algorithmic.data.xml.InfoOutput;
import edu.cornell.med.icb.goby.algorithmic.data.xml.SampleTotalCount;
import edu.cornell.med.icb.goby.alignments.*;
import edu.cornell.med.icb.goby.exception.GobyRuntimeException;
import edu.cornell.med.icb.goby.stats.DifferentialExpressionAnalysis;
import edu.cornell.med.icb.goby.stats.DifferentialExpressionCalculator;
import edu.cornell.med.icb.goby.stats.DifferentialExpressionResults;
import edu.cornell.med.icb.goby.stats.NormalizationMethod;
import edu.cornell.med.icb.goby.util.Timer;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.ints.IntSortedSet;
import it.unimi.dsi.fastutil.objects.*;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;

/**
 * Reads Goby alignments and genome annotations and output read counts that overlap with
 * the annotation segments.
 *
 * @author Fabien Campagne
 */
public class CompactAlignmentToAnnotationCountsMode extends AbstractGobyMode {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(CompactAlignmentToAnnotationCountsMode.class);

    /**
     * The mode name.
     */
    private static final String MODE_NAME = "alignment-to-annotation-counts";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Converts alignment to counts for "
            + "annotations (e.g., gene annotations or exons).";

    /**
     * The natural log of the number two.
     */
    private static final double LOG_2 = Math.log(2);

    /**
     * The output file.
     */
    private String outputFilename;
    /**
     * The annotation file.
     */
    private String annotationFile;

    private boolean filterByReferenceNames;
    private ObjectSet<String> includeReferenceNames;
    private ObjectOpenHashSet<String> includeAnnotationTypes;
    private String[] inputFilenames;

    private boolean writeAnnotationCounts = true;
    private boolean omitNonInformativeColumns;
    private String statsFilename;
    private ParallelTeam team;
    private boolean parallel;
    private boolean doComparison;
    private final DifferentialExpressionCalculator deCalculator = new DifferentialExpressionCalculator();
    private final DifferentialExpressionAnalysis deAnalyzer = new DifferentialExpressionAnalysis();
    /**
     * The set of normalization methods to use for the comparison.
     */

    private ObjectArraySet<NormalizationMethod> normalizationMethods;

    WeightParameters weightParams;
    private String includeReferenceNameCommas;
    private boolean filterAmbiguousReads;
    /**
     * When not null, filename for an xml formatted output file. This file will contain information
     * needed by stats mode.
     */
    private String infoOutputFilename;
    private boolean removeSharedSegments;
    /**
     * List of comparisons to perform.
     */
    private ArrayList<GroupComparison> groupComparisonsList;


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
     * @throws IOException   error parsing
     * @throws JSAPException error parsing
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args) throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);
        parseGenomicRange(jsapResult);

        parallel = jsapResult.getBoolean("parallel", false);
        writeAnnotationCounts = jsapResult.getBoolean("write-annotation-counts");
        omitNonInformativeColumns = jsapResult.getBoolean("omit-non-informative-columns");
        inputFilenames = jsapResult.getStringArray("input");
        outputFilename = jsapResult.getString("output");
        filterAmbiguousReads = jsapResult.getBoolean("filter-ambiguous-reads");
        removeSharedSegments = jsapResult.getBoolean("remove-shared-segments");

        if (filterAmbiguousReads) {
            System.out.println("Ambiguous reads will not be considered when estimating count statistics.");
        }
        statsFilename = jsapResult.getString("stats");
        final String groupsDefinition = jsapResult.getString("groups");
        deAnalyzer.parseGroupsDefinition(groupsDefinition, deCalculator, inputFilenames);
        final String compare = jsapResult.getString("compare");
        if (compare == null) {
            doComparison = false;
        } else {
            doComparison = true;
        }
        if (doComparison) {
            groupComparisonsList = deAnalyzer.parseCompare(compare);
        }
        deAnalyzer.setRunInParallel(parallel);
        includeReferenceNameCommas = jsapResult.getString("include-reference-names");
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
        normalizationMethods = deAnalyzer.parseNormalization(jsapResult);
        parseEval(jsapResult, deAnalyzer);
        infoOutputFilename = jsapResult.getString("info-output");

        weightParams = configureWeights(jsapResult);
        return this;
    }

    private void parseGenomicRange(JSAPResult jsapResult) {
        String startOffsetArgument = jsapResult.getString("start-position");
        String endOffsetArgument = jsapResult.getString("end-position");
        if (startOffsetArgument != null && endOffsetArgument == null ||
                endOffsetArgument != null && startOffsetArgument == null) {
            System.err.println("Start (-s/--start-position) and end offset (-e/--end-position) arguments must be specified together or not at all.");
            System.exit(1);
        }
        if (startOffsetArgument != null) {
            genomicRange = new GenomicRange();
            final String[] startTokens = startOffsetArgument.split("[,]");
            final String[] endTokens = endOffsetArgument.split("[,]");
            genomicRange.startPosition = Integer.parseInt(startTokens[1]);
            genomicRange.endPosition = Integer.parseInt(endTokens[1]);

            genomicRange.startChromosome = startTokens[0];
            genomicRange.endChromosome = endTokens[0];
        }

    }

    protected static WeightParameters configureWeights(final JSAPResult jsapResult) {
        final WeightParameters params = new WeightParameters();
        params.weightId = jsapResult.getString("use-weights");
        if (params.weightId == null || params.weightId.equals("false")) {
            params.useWeights = false;
        } else {
            params.useWeights = true;
        }

        params.formulaChoice = jsapResult.getString("adjust-gc-bias");
        params.adjustGcBias = !(params.formulaChoice == null || "false".equals(params.formulaChoice));
        if (params.formulaChoice != null) {
            params.formulaChoice = params.formulaChoice.toUpperCase();

        }
        if (params.adjustGcBias && !params.useWeights) {
            System.err.println("Cannot adjust bias when use-weights is false");
            System.exit(1);
        }
        if (params.useWeights) {
            if (!params.adjustGcBias) {
                System.out.println("Estimating expression as sum of weights");

            } else {
                System.out.println("Estimating expression with gc bias adjustment formula=" + params.formulaChoice);

            }
        }
        return params;
    }

    public static void parseEval(final JSAPResult jsapResult, final DifferentialExpressionAnalysis deAnalyzer) {
        final String evalString = jsapResult.getString("eval");

        final String[] evalArray = evalString.split(",");
        final ObjectSet<String> evalSet = new ObjectOpenHashSet<String>();
        for (final String evalName : evalArray) {
            evalSet.add(evalName.trim().toLowerCase().intern());
        }
        deAnalyzer.setEvalNames(evalSet);
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


    class BasenameParallelRegion extends ParallelRegion {
        private final Object2ObjectMap<String, ObjectList<Annotation>> allAnnots;
        private final String[] inputFiles;
        private final BufferedWriter writer;
        private final GenomicRange genomicRange;

        BasenameParallelRegion(final Object2ObjectMap<String, ObjectList<Annotation>> allAnnots,
                               final String[] inputFiles, final BufferedWriter writer, final GenomicRange genomicRange) {
            super();
            this.allAnnots = allAnnots;
            this.inputFiles = inputFiles;
            this.writer = writer;
            this.genomicRange = genomicRange;
        }

        @Override
        public void run() throws Exception {
            execute(0, inputFiles.length - 1 /* end index must be inclusive. This is counter-intuitive */, new IntegerForLoop() {

                @Override
                public void run(final int startIndex, final int endIndex) {
                    //   System.out.println(String.format("executing start= %d end=%d ",startIndex, endIndex));
                    for (int i = startIndex; i <= endIndex; ++i) {
                        if (i >= 0 && i < inputFilenames.length) {
                            final String inputBasename = AlignmentReaderImpl.getBasename(inputFiles[i]);
                            try {
                                processOneBasename(allAnnots, writer, inputFiles[i], inputBasename, genomicRange);
                                Runtime.getRuntime().gc();
                            } catch (IOException e) {
                                throw new GobyRuntimeException(e);
                            } finally {
                                IOUtils.closeQuietly(writer);
                            }
                        }
                    }
                }
            });
        }


    }

    protected synchronized ParallelTeam getParallelTeam() {
        if (team == null) {
            if (!parallel) {
                // 1 thread only is sequential
                team = new ParallelTeam(1);
            } else {
                // as many threads as configured with -Dpj.nt or default.
                team = new ParallelTeam();
            }
        }
        LOG.info("Executing on " + team.getThreadCount() + " threads.");
        return team;
    }

    GenomicRange genomicRange = null;

    /**
     * Run the mode.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        System.out.println("Reading annotations from " + annotationFile);
        if (genomicRange != null) {
            // resolve chromsome ids against the alignment headers:

            ConcatAlignmentReader concat = new ConcatAlignmentReader(AlignmentReaderImpl.getBasenames(inputFilenames));
            concat.readHeader();

            genomicRange.setTargetIds(concat.getTargetIdentifiers());
            concat.close();
        }
        final Object2ObjectMap<String, ObjectList<Annotation>> allAnnots = filterAnnotations(removeNonConstitutiveSegments(
                readAnnotations(annotationFile)), genomicRange);

        final Timer timer = new Timer();
        timer.start();
        BufferedWriter writer = null;
        try {

            if (outputFilename != null) {
                writer = new BufferedWriter(new FileWriter(outputFilename));
                writer.write("basename\tmain-id\tsecondary-id\ttype\tchro\tstrand\tlength\tstart\tend\tin-count\tover-count\tRPKM\tlog2(RPKM+1)\texpression\tnum-exons\n");
            }

            final BasenameParallelRegion region = new BasenameParallelRegion(allAnnots, inputFilenames, writer, genomicRange);

            try {
                getParallelTeam().execute(region);
            } catch (Exception e) {
                LOG.error("An exception occurred.", e);
            }
            Runtime.getRuntime().gc();
            Runtime.getRuntime().gc();
            if (doComparison) {
                final DifferentialExpressionResults results =
                        deAnalyzer.evaluateDifferentialExpressionStatistics(deCalculator, doComparison, normalizationMethods);
                PrintWriter statsWriter = null;
                try {
                    statsWriter = new PrintWriter(statsFilename);
                    results.write(statsWriter, '\t', deCalculator);
                } finally {
                    IOUtils.closeQuietly(statsWriter);
                }
            }


            writeInfoOutput();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        } finally {
            IOUtils.closeQuietly(writer);
        }
        timer.stop();
        System.out.println("time spent  " + timer.toString());

    }

    /**
     * Remove segments of an annotation when they overlap with other annotations. This keep only segments that uniquely
     * tag a gene/transcript. Single base overlaps are sufficient to trigger the exclusion of an entire segment (secondary id).
     *
     * @param annotations annotations read from disk
     * @return
     */
    private Object2ObjectMap<String, ObjectList<Annotation>> removeNonConstitutiveSegments(final Object2ObjectMap<String,
            ObjectList<Annotation>> annotations) {
        if (!removeSharedSegments) {
            return annotations;
        }
        final Int2ObjectMap<ObjectSet<String>> positionMap = new Int2ObjectOpenHashMap<ObjectSet<String>>();
        positionMap.defaultReturnValue(new ObjectOpenHashSet<String>());

        for (final String chromosome : annotations.keySet()) {
            positionMap.clear();
            for (final Annotation annotation : annotations.get(chromosome)) {
                for (final Segment element : annotation.getSegments()) {
                    for (int i = element.getStart(); i < element.getEnd(); i++) {
                        final ObjectSet<String> setOfGenes = positionMap.get(i);
                        setOfGenes.add(annotation.getId());
                        if (setOfGenes.size() > 1) {
                            annotation.remove(element);
                            if (annotation.getSegments().isEmpty()) {
                                // no more segments. Remove the annotation entirely.
                                annotations.get(chromosome).remove(annotation);
                            }
                            break;
                        }
                    }
                }
            }
        }
        return annotations;
    }

    InfoOutput infoOutInstance = new InfoOutput();

    private void writeInfoOutput() throws JAXBException, IOException {
        if (infoOutputFilename != null) {


            final JAXBContext jc = JAXBContext.newInstance(InfoOutput.class);

            final Marshaller m = jc.createMarshaller();

            for (String sampleId : deCalculator.samples()) {
                SampleTotalCount stc = new SampleTotalCount();
                stc.sampleId = sampleId;
                stc.totalCount = deCalculator.getNumAlignedInSample(sampleId);
                System.out.println(stc);
                infoOutInstance.totalCounts.add(stc);
            }
            // lengths have been added when filtering annotations.
            /*for (MutableString elementId : deCalculator.getElementIds()) {
                AnnotationLength ae = new AnnotationLength();
                ae.length = deCalculator.getElementLength(elementId);
                ae.id = elementId.toString();
                infoOutInstance.lengths.add(ae);
            } */
            FileWriter fileWriter = new FileWriter(infoOutputFilename);
            m.marshal(infoOutInstance, fileWriter);
            fileWriter.close();


        }
    }

    // Remove annotations that do not fully map within the genomic range.
    private Object2ObjectMap<String, ObjectList<Annotation>> filterAnnotations(Object2ObjectMap<String, ObjectList<Annotation>> map, GenomicRange genomicRange) {
        if (genomicRange == null) {
            return map;
        }
        Object2ObjectMap<String, ObjectList<Annotation>> filtered = new Object2ObjectArrayMap<String, ObjectList<Annotation>>();
        for (Map.Entry<String, ObjectList<Annotation>> entry : map.entrySet()) {
            String key = entry.getKey();
            String chromosome = key;

            for (final Annotation value : entry.getValue()) {

                chromosome = value.getChromosome();
                // convert to zero-based coordinates:
                int segmentStart = value.getStart() - 1;
                int segmentEnd = value.getEnd() - 1;
                if (genomicRange.fullyContains(chromosome, segmentStart, segmentEnd)) {

                    ObjectList<Annotation> chromosomeList = filtered.get(key);
                    if (chromosomeList == null) {
                        chromosomeList = new ObjectArrayList<Annotation>();
                        filtered.put(key, chromosomeList);
                    }
                    chromosomeList.add(value);
                }
                // Even when the element is out of GR, we keep  the element length, used to write complete info-output:
                AnnotationLength annotationLength = new AnnotationLength();
                annotationLength.id = value.getId();
                annotationLength.length = value.getLength();
                infoOutInstance.lengths.add(annotationLength);
            }
        }
        return filtered;
    }

    private void processOneBasename(final Object2ObjectMap<String, ObjectList<Annotation>> allAnnots,
                                    BufferedWriter writer, final String inputFile, final String inputBasename, GenomicRange range) throws IOException {

        WeightsInfo weights = null;
        if (weightParams.useWeights) {
            weights = loadWeights(inputBasename, weightParams.useWeights, weightParams.weightId);
            if (weights != null) {
                System.err.println("Weights have been provided and loaded and will be used to reweight transcript counts.");
            }
        }
        final AlignmentReaderFactory factory = filterAmbiguousReads ? new NonAmbiguousAlignmentReaderFactory() :
                new DefaultAlignmentReaderFactory();
        final int numberOfReferences;
        final DoubleIndexedIdentifier referenceIds;
        {
            final AlignmentReader reader = factory.createReader(inputBasename);
            reader.readHeader();
            numberOfReferences = reader.getNumberOfTargets();

            referenceIds = new DoubleIndexedIdentifier(reader.getTargetIdentifiers());
            reader.close();
        }
        System.out.println(String.format("Alignment contains %d reference sequences", numberOfReferences));
        if (genomicRange != null) {
            genomicRange.resolveChromosomeIndices(referenceIds);
        }
        final AnnotationCountIterateAlignments iterateAlignment = new AnnotationCountIterateAlignments();
        iterateAlignment.setWeightInfo(weightParams, weights);
        iterateAlignment.parseIncludeReferenceArgument(includeReferenceNameCommas);

        // Iterate through the alignment and retrieve algs:
        System.out.println("Loading alignment " + inputBasename + "..");
        iterateAlignment.setAlignmentReaderFactory(factory);
        iterateAlignment.iterate(genomicRange, inputBasename);

        final long numAlignedReadsInSample = iterateAlignment.getNumAlignedReadsInSample();
        final AnnotationCountInterface[] algs = iterateAlignment.getAlgs();
        final IntSortedSet referencesToProcess = iterateAlignment.getReferencesSelected();

        final String sampleId = FilenameUtils.getBaseName(inputBasename);

        deCalculator.setNumAlignedInSample(sampleId, numAlignedReadsInSample);

        if (outputFilename == null) {
            // output filename was not provided on the command line. We make one output per input basename
            if (writeAnnotationCounts) {
                final String outputFileTmp = FilenameUtils.removeExtension(inputFile) + ".ann-counts.tsv";
                writer = new BufferedWriter(new FileWriter(outputFileTmp));
                writer.write("basename\tmain-id\tsecondary-id\ttype\tchro\tstrand\tlength\tstart\tend\tin-count\tover-count\tRPKM\tlog2(RPKM+1)\texpression\tnum-exons\n");
            }
        }

        writeAnnotationCounts(allAnnots, writer, inputBasename, referenceIds, algs, referencesToProcess);

        if (outputFilename == null) {
            // output filename was not provided on the command line. We close each basename output.
            IOUtils.closeQuietly(writer);
        }

    }


    public static WeightsInfo loadWeights(final String inputBasename, final boolean useWeights, final String id) {
        WeightsInfo weights = null;
        if (useWeights) {
            try {
                weights = WeightsInfo.loadForBasename(inputBasename, id);
            } catch (Exception e) {
                LOG.warn("Cannot load weights file for " + inputBasename, e);
                LOG.warn("Using weight=1.0 for each read.");
                return new WeightsInfo() {
                    @Override
                    public float getWeight(final int readIndex) {
                        return 1f;
                    }
                };
            }
        }
        return weights;
    }


    private void writeAnnotationCounts(final Object2ObjectMap<String, ObjectList<Annotation>> allAnnots,
                                       final BufferedWriter writer, final String inputBasename,
                                       final DoubleIndexedIdentifier referenceIds, final AnnotationCountInterface[] algs,
                                       final IntSet referencesToProcess) throws IOException {

        // collect all element ids:
        if (doComparison) {
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
            deCalculator.reserve(numberOfElements, inputFilenames.length);
        }

        int numberOfAnnotationCountsWritten = 0;
        for (final int referenceIndex : referencesToProcess) {
            final String chromosomeName = referenceIds.getId(referenceIndex).toString();
            System.out.println("Writing annotation counts for reference " + chromosomeName);

            if (!allAnnots.containsKey(chromosomeName)) {
                continue;
            }

            final ObjectList<Annotation> annots = allAnnots.get(chromosomeName);
            algs[referenceIndex].sortReads();
            algs[referenceIndex].accumulate();
            algs[referenceIndex].baseCount();
            if (doComparison) {
                for (final Annotation annot : annots) {
                    final String geneID = annot.getId();
                    deCalculator.defineElement(geneID);
                }
            }

            // get just the filename (strip the path, not the extension)
            final String basename = FilenameUtils.getName(inputBasename);
            final String sampleId = inputBasename;
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

                    final double geneRPKM = deCalculator.calculateNormalized(geneOverlapReads, annot.getLength(),
                            deCalculator.getNumAlignedInSample(sampleId));
                    if (writeAnnotationCounts) {
                        writer.write(String.format("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%g\t%g\t%g\t%g\t%g\t%d%n",
                                basename,
                                geneID,
                                "",
                                "gene",
                                annot.getChromosome(),
                                annot.getStrand(),
                                geneLength,
                                geneStart,
                                geneEnd,
                                geneInsideReads,
                                geneOverlapReads,
                                geneRPKM,
                                log2(geneRPKM),
                                geneExpression,
                                numExons));

                    }
                    numberOfAnnotationCountsWritten++;
                    if (doComparison) {
                        deCalculator.observe(basename, geneID, geneExpression);
                    }
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
                            if (writeAnnotationCounts) {
                                writer.write(String.format("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%g\t%g\t%g\t%g%n",
                                        basename,
                                        geneID,
                                        exonID,
                                        "exon",
                                        annot.getChromosome(),
                                        exonStrand,
                                        exonLength,
                                        exonStart,
                                        exonEnd,
                                        exonInsideReads,
                                        exonOverlapReads,
                                        exonRPKM,
                                        log2(exonRPKM)));
                            }
                            numberOfAnnotationCountsWritten++;
                            if (doComparison && includeAnnotationTypes.contains("exon")) {
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
                                    if (writeAnnotationCounts) {
                                        writer.write(String.format("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%g\t%g\t%e\t%g%n",
                                                basename,
                                                geneID,
                                                intronID,
                                                "other",
                                                annot.getChromosome(),
                                                exonStrand,
                                                exonLength,
                                                exonStart,
                                                exonEnd,
                                                intronInsideReads,
                                                intronOverlapReads,
                                                intronRPKM,
                                                log2(intronRPKM)));
                                    }
                                    numberOfAnnotationCountsWritten++;
                                    if (doComparison && includeAnnotationTypes.contains("other")) {
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

        LOG.info("Wrote " + numberOfAnnotationCountsWritten + " entries");
        if (numberOfAnnotationCountsWritten == 0) {
            LOG.warn("No entries were written.  This may be due to the fact that names "
                    + "in the reference dataset used do not match those in the annotation file.  "
                    + "For example, ENSEMBL names chromosomes \"1\",\"2\",\"3\" whereas UCSC "
                    + "names the same chromosomes \"chr1\",\"chr2\",\"chr3\". In these "
                    + "cases you will need to adjust the names in the annotation file being used "
                    + "so they match the names used in the reference dataset.");
        }
    }

    /**
     * Calculate the log2 of x+1.
     *
     * @param x
     * @return log2(x+1)=Math.log1p(x)/Math.log(2)
     */
    private double log2(final double x) {
        return Math.log1p(x) / LOG_2;
    }


    /**
     * Read tab delimited annotation file including 6 columns : chromosome, strand, transcriptID,
     * segmentID, start, end each row is a segment, read each chromosome at a time.
     *
     * @param annotFile
     * @return
     * @throws IOException
     */
    public static Object2ObjectMap<String, ObjectList<Annotation>> readAnnotations(final String annotFile) throws IOException {
        return readAnnotations(new FileReader(annotFile));
    }

    /**
     * Read tab delimited annotation file including 6 columns : chromosome, strand, transcriptID,
     * segmentID, start, end each row is a segment, read each chromosome at a time.
     *
     * @param annotReader reader over an annotation file.
     * @return
     * @throws IOException
     */
    public static Object2ObjectMap<String, ObjectList<Annotation>> readAnnotations(final Reader annotReader) throws IOException {

        BufferedReader reader = null;
        final Object2ObjectMap<String, Annotation> annots = new Object2ObjectAVLTreeMap<String, Annotation>();
        try {
            reader = new BufferedReader(annotReader);
            String line;
            final String header = reader.readLine();

            while ((line = reader.readLine()) != null) {
                if (!line.startsWith("#")) {
                    final String[] linearray = line.trim().split("\t");
                    if (linearray.length < 6) {
                        LOG.warn("Annotation file, encountered truncated line, ignoring: " + line);
                        continue;
                    }
                    final String chromosome = linearray[0];
                    //           if(!chromosome.equalsIgnoreCase(chroName)) continue;
                    final String strand = linearray[1];
                    final String transcriptID = linearray[2];
                    final String exonID = linearray[3];
                    final int segmentStart = Integer.parseInt(linearray[4]);
                    final int segmentEnd = Integer.parseInt(linearray[5]);
                    final Segment segment = new Segment(segmentStart, segmentEnd, exonID, strand);
                    boolean ignoreElement = false;


                    // annotation segment must be within the genomic range to be considered. Consider all segments when
                    // we are processing the entire alignment (genomicRange==null)
                    if (annots.containsKey(transcriptID)) {
                        annots.get(transcriptID).addSegment(segment);
                    } else {
                        final Annotation annot = new Annotation(transcriptID, chromosome, strand);
                        annot.addSegment(segment);
                        annots.put(transcriptID, annot);
                    }

                }
            }
        } finally {
            IOUtils.closeQuietly(reader);
        }

        //Group annotations by chromosome (key of the map returned)
        final Object2ObjectMap<String, ObjectList<Annotation>> allAnnots
                = new Object2ObjectOpenHashMap<String, ObjectList<Annotation>>();
        for (final Object2ObjectMap.Entry<String, Annotation> entry : annots.object2ObjectEntrySet()) {
            final Annotation annotation = entry.getValue();
            annotation.sortSegments();
            final String chromosome = annotation.getChromosome();
            if (allAnnots.containsKey(chromosome)) {
                allAnnots.get(chromosome).add(annotation);
            } else {
                final ObjectList<Annotation> annotations = new ObjectArrayList<Annotation>();
                annotations.add(annotation);
                allAnnots.put(chromosome, annotations);
            }
        }
        return allAnnots;
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
        new CompactAlignmentToAnnotationCountsMode().configure(args).execute();
    }

}
