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

import cern.colt.Timer;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.algorithmic.algorithm.AnnotationCount;
import edu.cornell.med.icb.goby.algorithmic.data.Annotation;
import edu.cornell.med.icb.goby.algorithmic.data.AnnotationRPKM;
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
import it.unimi.dsi.fastutil.objects.ObjectListIterator;
import org.apache.commons.io.IOUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Set;

/**
 * Aggregates peaks according to their interpeak RPKM differences given a threshold.
 *
 * @author Jaaved Mohammed
 */
public class AggregatePeaksByRPKMDifferenceMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "aggregate-by-peak-rpkm-difference";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Aggregates peaks according to their interpeak RPKM differences given a threshold.";

    /**
     * The input files of compact reads.
     */
    private String[] inputFilenames;

    /**
     * The output file.
     */
    private String outputFile;

    /**
     * The annotation-like file containing the union of all peaks accross all samples.
     */
    private String unionAnnotationFile;
    private int rpkmDiffThreshold;
    private String rpkmFileName;


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

        inputFilenames = jsapResult.getStringArray("input");
        outputFile = jsapResult.getString("output");
        unionAnnotationFile = jsapResult.getString("peak-union");
        rpkmDiffThreshold = jsapResult.getInt("threshold");
        rpkmFileName = jsapResult.getString("rpkm-file");

        return this;
    }

    /**
     * Run the mode.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        System.out.println("Reading annotations from " + unionAnnotationFile);
        final Object2ObjectMap<String, ObjectList<AnnotationRPKM>> allAnnots = readAnnotations(unionAnnotationFile);
        final Timer timer = new Timer();
        timer.start();

        //For all input files, and given the union peaks, compute the average RPKM for each peak accross all input samples
        for (String inputFile : inputFilenames) {
            System.out.println("Reading alignment file: " + inputFile);
            String inputBasename = AlignmentReader.getBasename(inputFile);
            final AlignmentReader reader = new AlignmentReader(inputBasename);
            reader.readHeader();
            final int numberOfReferences = reader.getNumberOfTargets();

            final DoubleIndexedIdentifier referenceIds = new DoubleIndexedIdentifier(reader.getTargetIdentifiers());
            reader.close();
            final AnnotationCount[] algs = new AnnotationCount[numberOfReferences];
            final IntSet referencesToProcess = new IntOpenHashSet();

            Set<String> chromosomeNames = allAnnots.keySet();

            // create count writers, one for each reference sequence in the alignment:
            for (int referenceIndex = 0; referenceIndex < numberOfReferences; referenceIndex++) {
                final String referenceName = referenceIds.getId(referenceIndex).toString();
                if (chromosomeNames.contains(referenceName)) {
                    // subset of reference names selected by the command line:
                    referencesToProcess.add(referenceIndex);
                    algs[referenceIndex] = new AnnotationCount();
                    algs[referenceIndex].baseCounter.startPopulating();
                }
            }


            final AlignmentReader referenceReader = new AlignmentReader(inputBasename);
            referenceReader.readHeader();

            // read the alignment:
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
            computeAnnotationRPKM(allAnnots, referenceIds, algs, referencesToProcess, numAlignedReadsInSample);
        }

        //Now that we've finished processing all samples, compute the average of the RPKMs accross all samples
        //Also if a rpkm file is specified, write the averaged RPKM values to this file.
        int numSamples = inputFilenames.length;
        Set<String> references = allAnnots.keySet();
        if (rpkmFileName != null) {
            System.out.println("Writing the RPKM values to file: " + rpkmFileName);
        }

        for (final String reference : references) {
            System.out.println("Averaging RPKM for reference " + reference);
            final ObjectList<AnnotationRPKM> annots = allAnnots.get(reference);

            for (final AnnotationRPKM annot : annots) {
                annot.RPKM /= numSamples;
            }

            if (rpkmFileName != null) {
                writeRPKM(rpkmFileName, annots, true /*append*/);
            }
        }

        //FINALLY, we can perform the aggregation
        for (final String reference : references) {
            System.out.println("Aggregating annotations by RPKM for reference " + reference);
            final ObjectList<AnnotationRPKM> annots = allAnnots.get(reference);
            ObjectList<Annotation> mergedAnnots = new ObjectArrayList<Annotation>();
            ObjectListIterator<AnnotationRPKM> annotsIterator = annots.listIterator();

            AnnotationRPKM annot = null, nextAnnot = null;
            if (annotsIterator.hasNext()) {
                annot = annotsIterator.next();
            }
            while (annot != null) {
                final Annotation newAnnot = new Annotation(annot.chromosome + "." + annot.getStart(), annot.chromosome);
                newAnnot.segments.addAll(annot.segments);
                if (annotsIterator.hasNext()) {
                    nextAnnot = annotsIterator.next();
                    while ((nextAnnot != null) && Math.abs(annot.RPKM - nextAnnot.RPKM) <= rpkmDiffThreshold) {
                        //Merge the two annotations
                        newAnnot.segments.addAll(nextAnnot.segments);
                        if (annotsIterator.hasNext()) {
                            nextAnnot = annotsIterator.next();
                        } else {
                            nextAnnot = null;
                        }
                    }
                }
                newAnnot.sortSegments();
                mergedAnnots.add(newAnnot);
                annot = nextAnnot;
                nextAnnot = null;
            }

            //clean out the old list to preserve memory explosion
            annots.clear();

            //Finally write the new annotations to the output file
            writeAnnotations(outputFile, mergedAnnots, true /*append*/);
        }

        //DONE
        System.out.println("time spent  " + timer.toString());
    }


    private void computeAnnotationRPKM(Object2ObjectMap<String, ObjectList<AnnotationRPKM>> allAnnots, DoubleIndexedIdentifier referenceIds, AnnotationCount[] algs, IntSet referencesToProcess, int numAlignedReadsInSample) throws IOException {
        for (final int referenceIndex : referencesToProcess) {

            final String chromosomeName = referenceIds.getId(referenceIndex).toString();

            if (!allAnnots.containsKey(chromosomeName)) {
                continue;
            }
            final ObjectList<AnnotationRPKM> annots = allAnnots.get(chromosomeName);
            algs[referenceIndex].sortReads();
            algs[referenceIndex].baseCounter.accumulate();
            algs[referenceIndex].baseCounter.baseCount();

            for (final AnnotationRPKM annot : annots) {
                final int start = annot.getStart();
                final int end = annot.getEnd();
                final int overlapReads = algs[referenceIndex].countReadsPartiallyOverlappingWithInterval(start, end);

                final double RPKM = calculateRPKM(overlapReads, annot.getLength(), numAlignedReadsInSample);
                annot.RPKM += RPKM;
            }
        }
    }

    private float calculateRPKM(final int readCountInt, final int annotLength, final int geneOverlapReads) {
        final float readCount = readCountInt;
        final float length = annotLength; // in bases
        final float sampleReadCount = geneOverlapReads; // in read
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
    public Object2ObjectMap<String, ObjectList<AnnotationRPKM>> readAnnotations(final String annotFile) throws IOException {

        BufferedReader reader = null;
        final Object2ObjectMap<String, AnnotationRPKM> annots = new Object2ObjectOpenHashMap<String, AnnotationRPKM>();
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
                        final AnnotationRPKM annot = new AnnotationRPKM(transcriptID, chromosome, 0d);
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
        final Object2ObjectMap<String, ObjectList<AnnotationRPKM>> allAnnots
                = new Object2ObjectOpenHashMap<String, ObjectList<AnnotationRPKM>>();
        for (final Object2ObjectMap.Entry<String, AnnotationRPKM> entry : annots.object2ObjectEntrySet()) {
            entry.getValue().sortSegments();
            final String chromosome = entry.getValue().chromosome;
            if (allAnnots.containsKey(chromosome)) {
                allAnnots.get(chromosome).add(entry.getValue());
            } else {
                final ObjectList<AnnotationRPKM> annotations = new ObjectArrayList<AnnotationRPKM>();
                annotations.add(entry.getValue());
                allAnnots.put(chromosome, annotations);
            }
        }
        return allAnnots;
    }

    public static void writeRPKM(String outputFileName, ObjectList<AnnotationRPKM> annotationList, boolean append) {
        PrintWriter writer = null;
        File outputFile = new File(outputFileName);

        try {
            if (!outputFile.exists()) {
                writer = new PrintWriter(outputFile);

                //Write the file header
                writer.write("Chromosome_Name\tStrand\tPrimary_ID\tSecondary_ID\tTranscript_Start\tTranscript_End\tAverage_RPKM\n");

            } else {
                writer = new PrintWriter(new FileOutputStream(outputFile, append));
            }

            if (writer != null) {
                ObjectListIterator<AnnotationRPKM> annotIterator = annotationList.listIterator();
                while (annotIterator.hasNext()) {
                    final AnnotationRPKM annotation = annotIterator.next();
                    annotation.write(writer);
                }
            } else {
                System.err.println("Cannot write annotations to file: " + outputFileName);
                System.err.println("The writer failed to initialize.");
                System.exit(1);
            }
        }
        catch (FileNotFoundException fnfe) {
            System.err.println("Caught exception in writeAnnotations: " + fnfe.getMessage());
            System.exit(1);
        }
        finally {
            IOUtils.closeQuietly(writer);
        }
    }

    public static void writeAnnotations(String outputFileName, ObjectList<Annotation> annotationList, boolean append) {
        PrintWriter writer = null;
        File outputFile = new File(outputFileName);

        try {
            if (!outputFile.exists()) {
                writer = new PrintWriter(outputFile);

                //Write the file header

                writer.write("Chromosome_Name\tStrand\tPrimary_ID\tSecondary_ID\tTranscript_Start\tTranscript_End\n");

            } else {
                writer = new PrintWriter(new FileOutputStream(outputFile, append));
            }

            if (writer != null) {
                ObjectListIterator<Annotation> annotIterator = annotationList.listIterator();
                while (annotIterator.hasNext()) {

                    final Annotation annotation = annotIterator.next();
                    annotation.write(writer);
                }
            } else {
                System.err.println("Cannot write annotations to file: " + outputFileName);
                System.err.println("The writer failed to initialize.");
                System.exit(1);
            }
        }
        catch (FileNotFoundException fnfe) {
            System.err.println("Caught exception in writeAnnotations: " + fnfe.getMessage());
            System.exit(1);
        }
        finally {
            IOUtils.closeQuietly(writer);
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
        new AggregatePeaksByRPKMDifferenceMode().configure(args).execute();
    }
}
