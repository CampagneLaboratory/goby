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
import edu.cornell.med.icb.goby.algorithmic.data.Segment;
import edu.cornell.med.icb.goby.algorithmic.data.WeightsInfo;
import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
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
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.Object2ObjectMap;
import it.unimi.dsi.fastutil.objects.Object2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

/**
 * Reads a compact alignment and genome annotations and output read counts that overlap with
 * the annotation segments.
 *
 * @author Xutao Deng
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

        parallel = jsapResult.getBoolean("parallel", false);
        writeAnnotationCounts = jsapResult.getBoolean("write-annotation-counts");
        omitNonInformativeColumns = jsapResult.getBoolean("omit-non-informative-columns");
        inputFilenames = jsapResult.getStringArray("input");
        outputFilename = jsapResult.getString("output");


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
            deAnalyzer.parseCompare(compare);
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


        weightParams = configureWeights(jsapResult);
        return this;
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

        BasenameParallelRegion(final Object2ObjectMap<String, ObjectList<Annotation>> allAnnots,
                               final String[] inputFiles, final BufferedWriter writer) {
            super();
            this.allAnnots = allAnnots;
            this.inputFiles = inputFiles;
            this.writer = writer;
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
                                processOneBasename(allAnnots, writer, inputFiles[i], inputBasename);
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


    /**
     * Run the mode.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        System.out.println("Reading annotations from " + annotationFile);
        final Object2ObjectMap<String, ObjectList<Annotation>> allAnnots = readAnnotations(annotationFile);
        final Timer timer = new Timer();
        timer.start();
        BufferedWriter writer = null;
        try {

            if (outputFilename != null) {
                writer = new BufferedWriter(new FileWriter(outputFilename));
                writer.write("basename\tmain-id\tsecondary-id\ttype\tchro\tstrand\tlength\tstart\tend\tin-count\tover-count\tRPKM\tlog2(RPKM+1)\texpression\tnum-exons\n");
            }
            final BasenameParallelRegion region = new BasenameParallelRegion(allAnnots, inputFilenames, writer);

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
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        } finally {
            IOUtils.closeQuietly(writer);
        }
        timer.stop();
        System.out.println("time spent  " + timer.toString());
    }

    private void processOneBasename(final Object2ObjectMap<String, ObjectList<Annotation>> allAnnots,
                                    BufferedWriter writer, final String inputFile, final String inputBasename) throws IOException {

        WeightsInfo weights = null;
        if (weightParams.useWeights) {
            weights = loadWeights(inputBasename, weightParams.useWeights, weightParams.weightId);
            if (weights != null) {
                System.err.println("Weights have been provided and loaded and will be used to reweight transcript counts.");
            }
        }

        final AlignmentReaderImpl reader = new AlignmentReaderImpl(inputBasename);
        reader.readHeader();
        final int numberOfReferences = reader.getNumberOfTargets();

        final DoubleIndexedIdentifier referenceIds = new DoubleIndexedIdentifier(reader.getTargetIdentifiers());
        reader.close();
        System.out.println(String.format("Alignment contains %d reference sequences", numberOfReferences));

        final AnnotationCountIterateAlignments iterateAlignment = new AnnotationCountIterateAlignments();
        iterateAlignment.setWeightInfo(weightParams, weights);
        iterateAlignment.parseIncludeReferenceArgument(includeReferenceNameCommas);


        final AlignmentReader referenceReader = new AlignmentReaderImpl(inputBasename);
        referenceReader.readHeader();

        // ITerate through the alignment and retrieve algs:
        System.out.println("Loading alignment " + inputBasename + "..");
        iterateAlignment.iterate(inputBasename);

        final int numAlignedReadsInSample = iterateAlignment.getNumAlignedReadsInSample();
        final AnnotationCountInterface[] algs = iterateAlignment.getAlgs();
        final IntSet referencesToProcess = iterateAlignment.getReferencesSelected();

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

        int numberOfAnottationCountsWritten = 0;
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
                    numberOfAnottationCountsWritten++;
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
                            numberOfAnottationCountsWritten++;
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
                                    numberOfAnottationCountsWritten++;
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

        BufferedReader reader = null;
        final Object2ObjectMap<String, Annotation> annots = new Object2ObjectOpenHashMap<String, Annotation>();
        try {
            reader = new BufferedReader(new FileReader(annotFile));
            String line;
            final String header = reader.readLine();

            while ((line = reader.readLine()) != null) {
                if (!line.startsWith("#")) {
                    final String[] linearray = line.trim().split("\t");
                    final String chromosome = linearray[0];
                    //           if(!chromosome.equalsIgnoreCase(chroName)) continue;
                    final String strand = linearray[1];
                    final String transcriptID = linearray[2];
                    final String exonID = linearray[3];
                    final int segmentStart = Integer.parseInt(linearray[4]);
                    final int segmentEnd = Integer.parseInt(linearray[5]);
                    final Segment segment = new Segment(segmentStart, segmentEnd, exonID, strand);
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

        //organize the Annotations to chromosome
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
