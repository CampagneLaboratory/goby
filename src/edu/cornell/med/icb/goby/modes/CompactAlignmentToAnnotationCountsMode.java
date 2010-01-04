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

import cern.colt.Timer;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.algorithmic.algorithm.AnnotationCount;
import edu.cornell.med.icb.goby.algorithmic.data.Annotation;
import edu.cornell.med.icb.goby.algorithmic.data.Segment;
import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.Object2ObjectMap;
import it.unimi.dsi.fastutil.objects.Object2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.apache.commons.io.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

/**
 * Receive reads file and annotation files and output counts for gene or exons.
 *
 * @author Xutao Deng
 */
public class CompactAlignmentToAnnotationCountsMode extends AbstractGobyMode {

    /**
     * The mode name.
     */
    public static final String MODE_NAME = "alignment-to-annotation-counts";
    public static final String MODE_DESCRIPTION = "Converts alignment to counts for annotations (e.g., gene transcript annotations or exons).";


    /**
     * The input file.
     */
    private String inputFile;

    /**
     * The output file.
     */
    private String outputFile;

    private String annotationFile;

    boolean filterByReferenceNames;
    private ObjectSet<String> includeReferenceNames;
    private ObjectOpenHashSet<String> includeAnnotationTypes;

    public String getModeName() {
        return MODE_NAME;
    }

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

        inputFile = AlignmentReader.getBasename(jsapResult.getString("input"));
        outputFile = jsapResult.getString("output");
        annotationFile = jsapResult.getString("annotation");
        final String includeReferenceNameComas = jsapResult.getString("include-reference-names");
        includeReferenceNames = new ObjectOpenHashSet<String>();
        if (includeReferenceNameComas != null) {

            includeReferenceNames.addAll(Arrays.asList(includeReferenceNameComas.split("[,]")));
            System.out.println("Will write counts for the following sequences:");
            for (final String name : includeReferenceNames) {
                System.out.println(name);
            }
            filterByReferenceNames = true;
        }
        final String includeAnnotationTypeComas = jsapResult.getString("include-annotation-types");
        includeAnnotationTypes = new ObjectOpenHashSet<String>();
        if (includeAnnotationTypeComas != null) {

            includeAnnotationTypes.addAll(Arrays.asList(includeAnnotationTypeComas.split("[,]")));
            System.out.println("Will write counts for the following annotation types:");
            for (final String name : includeAnnotationTypes) {
                System.out.println(name);
            }
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
        final Object2ObjectMap<String, ObjectList<Annotation>> allAnnots = readAnnotations(annotationFile);
        final AlignmentReader reader = new AlignmentReader(inputFile);
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
                algs[referenceIndex].baseCounter.startPopulating();
            }
        }


        final AlignmentReader referenceReader = new AlignmentReader(inputFile);
        referenceReader.readHeader();

        // read the alignment:
        System.out.println("Loading the alignment..");
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

        final Timer timer = new Timer();
        timer.start();
        BufferedWriter writer = null;
        try {
            writer = new BufferedWriter(new FileWriter(outputFile));
            writer.write("basename\tmain-id\tsecondary-id\ttype\tchro\tstrand\tlength\tstart\tend\tin-count\tover-count\tRPKM\tlog2(PRKM+1)\texpression\tnum-exons\n");

            //       System.out.println("id\ttype\tchro\tstart\tend\tin_count\tover_count\tdepth\texpression");
            for (final int referenceIndex : referencesToProcess) {

                final String chromosomeName = referenceIds.getId(referenceIndex).toString();

                System.out.println("Writing annotation counts for reference " + chromosomeName);

                if (!allAnnots.containsKey(chromosomeName)) {
                    continue;
                }
                final ObjectList<Annotation> annots = allAnnots.get(chromosomeName);
                algs[referenceIndex].sortReads();
                algs[referenceIndex].baseCounter.accumulate();
                algs[referenceIndex].baseCounter.baseCount();

                for (final Annotation annot : annots) {
                    final String geneID = annot.id;
                    String basename = inputFile;
                    if (includeAnnotationTypes.contains("gene")) {

                        final int geneStart = annot.getStart();
                        final int geneEnd = annot.getEnd();
                        final int geneLength = geneEnd - geneStart + 1;
                        final float geneDepth = algs[referenceIndex].averageReadsPerPosition(geneStart, geneEnd);
                        final int geneOverlapReads = algs[referenceIndex].readsOverlapSegmentCount(geneStart, geneEnd);
                        final int geneInsideReads = algs[referenceIndex].readsInSegmentCount(geneStart, geneEnd);
                        final int geneExpression = algs[referenceIndex].geneExpressionCount(annot);
                        final int numExons = annot.segments.size();


                        double geneRPKM = calculateRPKM(geneOverlapReads, annot.getLength(), numAlignedReadsInSample);
                        writer.write(String.format("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%g\t%g\t%d\t%d\t%n",
                                basename,
                                geneID,
                                "",
                                "gene",
                                annot.chromosome,
                                annot.strand,
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
                    final int numberExons = annot.segments.size();
                    final int numberIntrons = numberExons - 1;

                    for (int i = 0; i < numberExons; i++) {
                        final Segment segment = annot.segments.get(i);
                        final int exonStart = segment.start;
                        final int exonEnd = segment.end;

                        final String exonStrand = segment.strand;
                        final int exonLength = segment.getLength();
                        final String exonID = segment.id;
                        final float exonDepth = algs[referenceIndex].averageReadsPerPosition(exonStart, exonEnd);
                        final int exonOverlapReads = algs[referenceIndex].readsOverlapSegmentCount(exonStart, exonEnd);
                        final int exonInsideReads = algs[referenceIndex].readsInSegmentCount(exonStart, exonEnd);
                        double exonRPKM = calculateRPKM(exonOverlapReads, segment.getLength(), numAlignedReadsInSample);
                        if (includeAnnotationTypes.contains("exon")) {
                            writer.write(String.format("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%g\t%g\t\t%n",
                                    basename,
                                    geneID,
                                    exonID,
                                    "exon",
                                    annot.chromosome,
                                    exonStrand,
                                    exonLength,
                                    exonStart,
                                    exonEnd,
                                    exonInsideReads,
                                    exonOverlapReads,
                                    exonRPKM,
                                    log2(exonRPKM)));
                        }
                        if (i < numberIntrons) {
                            final int intronStart = segment.end + 1;
                            final Segment intronSegment = annot.segments.get(i + 1);
                            final int intronEnd = intronSegment.start - 1;
                            final int intronLength = intronEnd - intronStart + 1;
                            final String intronID = segment.id + "-" + intronSegment.id;
                            final float intronDepth = algs[referenceIndex].averageReadsPerPosition(intronStart, intronEnd);
                            final int intronOverlapReads = algs[referenceIndex].readsOverlapSegmentCount(intronStart, intronEnd);
                            final int intronInsideReads = algs[referenceIndex].readsInSegmentCount(intronStart, intronEnd);
                            double intronRPKM = calculateRPKM(intronOverlapReads, intronSegment.getLength(), numAlignedReadsInSample);
                            if (intronLength > 0) {
                                if (includeAnnotationTypes.contains("intron"))
                                    writer.write(String.format("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%e\t%g\t\t%n",
                                            basename,
                                            geneID,
                                            intronID,
                                            "intron",
                                            annot.chromosome,
                                            exonStrand,
                                            exonLength,
                                            exonStart,
                                            exonEnd,
                                            intronInsideReads,
                                            intronOverlapReads,
                                            intronRPKM,
                                            log2(intronRPKM)));
                            }
                        }
                    }
                }
                algs[referenceIndex] = null;
            }
        } finally {
            IOUtils.closeQuietly(writer);
        }

        timer.stop();
        System.out.println("time spent  " + timer.millis());
    }

    final double LOG_2 = Math.log(2);

    /**
     * Calculate the log2 of x +1.
     * @param x
     * @return log2(x+1)=Math.log1p(x)/Math.log(2)
     */
    private double log2(double x) {

        return Math.log1p(x) / LOG_2;
    }

    private double calculateRPKM(int readCountInt, int annotLength, int geneOverlapReads) {
        double readCount = readCountInt;
        double length = annotLength; // in bases
        double sampleReadCount = geneOverlapReads; // in read
        return readCount / (length / 1000.0f) / (sampleReadCount / 1E6f);
    }

    /**
     * Read tab delimited annotation file including 6 columns : chromosome, strand, transcriptID,
     * segmentID, start, end each row is a segment, read each chromosome at a time.
     *
     * @param annotFile
     * @return
     * @throws IOException
     */
    public Object2ObjectMap<String, ObjectList<Annotation>> readAnnotations(final String annotFile) throws IOException {

        BufferedReader reader = null;
        final Object2ObjectMap<String, Annotation> annots = new Object2ObjectOpenHashMap<String, Annotation>();
        try {
            reader = new BufferedReader(new FileReader(annotFile));
            String line;
            final String header = reader.readLine();
            System.out.println(header);
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
                        final Annotation annot = new Annotation(transcriptID, chromosome);
                        annot.strand = strand;
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
            entry.getValue().sortSegments();
            final String chromosome = entry.getValue().chromosome;
            if (allAnnots.containsKey(chromosome)) {
                allAnnots.get(chromosome).add(entry.getValue());
            } else {
                final ObjectList<Annotation> annotations = new ObjectArrayList<Annotation>();
                annotations.add(entry.getValue());
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
