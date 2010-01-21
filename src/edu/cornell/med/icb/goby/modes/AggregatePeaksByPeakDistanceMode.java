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
import edu.cornell.med.icb.goby.algorithmic.data.Segment;
import edu.cornell.med.icb.goby.algorithmic.data.Annotation;
import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.ConcatAlignmentReader;
import edu.cornell.med.icb.goby.counts.AnyTransitionCountsIterator;
import edu.cornell.med.icb.goby.counts.CountsArchiveReader;
import edu.cornell.med.icb.goby.counts.CountsReaderI;
import edu.cornell.med.icb.goby.counts.Peak;
import edu.cornell.med.icb.goby.counts.PeakAggregator;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.*;

import java.io.*;
import java.util.Arrays;
import java.util.Set;
import java.util.Collections;
import java.util.Iterator;

import org.apache.commons.io.IOUtils;


/**
 * Write annotations corresponding to consensus peaks found in each sequence of count archives.
 *
 * @author Jaaved Mohammed
 */
public class AggregatePeaksByPeakDistanceMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    public static final String MODE_NAME = "aggregate-by-peak-distance";
    public static final String MODE_DESCRIPTION = "Write annotations corresponding to consensus peaks found in each sequence of count archives.";


    /**
     * The input files. Must reduce to alignment basenames.
     */
    private String inputFile;

    /**
     * The output file.
     */
    private String outputFile;

    /**
     * Other program arguments.
     */
    private int peakDistThreshold;

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
        outputFile = jsapResult.getString("output");
        peakDistThreshold = jsapResult.getInt("threshold");

        return this;
    }

    /**
     * Run the mode.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        System.out.println("Reading official proposal annotations from " + inputFile);
        final Object2ObjectMap<String, ObjectList<Annotation>> unionAnnots = readAnnotations(inputFile);

        //In the union file, all peaks should be non-overlapping 
        System.out.println("Aggregating peaks by distance.");
        final Set<String> chromsomes = unionAnnots.keySet();
        ObjectList<Annotation> annotationList = null;
        for (final String chromsome : chromsomes) {
            annotationList = new ObjectArrayList<Annotation>();

            //Get the list of annotations and rip out the segments
            final ObjectList<Segment> segmentsList = new ObjectArrayList<Segment>();
            final ObjectList<Annotation> unionAnnotList = unionAnnots.get(chromsome);
            System.out.println(String.format("Reference %s: Starting with %d annotations.", chromsome, unionAnnotList.size()));
            ObjectListIterator<Annotation> unionAnnotListIt = unionAnnotList.listIterator();
            while (unionAnnotListIt.hasNext()) {
                segmentsList.addAll(unionAnnotListIt.next().segments);
            }
            unionAnnotList.clear();
            Collections.sort(segmentsList);
            ObjectListIterator<Segment> segIterator = segmentsList.listIterator();

            Segment segment = null, nextSeg = null;
            if (segIterator.hasNext()) {
                segment = segIterator.next();
            }
            while (segment != null) {
                final Annotation annotation = new Annotation(chromsome + "." + segment.start + "." + segment.getLength(), chromsome);
                annotation.addSegment(segment);
                if (segIterator.hasNext()) {
                    nextSeg = segIterator.next();
                    while ((nextSeg != null) && segment.distance(nextSeg) <= peakDistThreshold) {
                        //Add the nextSegment to the annotation
                        annotation.addSegment(nextSeg);
                        if (segIterator.hasNext()) {
                            nextSeg = segIterator.next();
                        } else {
                            nextSeg = null;
                        }
                    }
                }
                annotationList.add(annotation);
                segment = nextSeg;
                nextSeg = null;
            }
            System.out.println(String.format("Reference %s: Finished with %d annotations.", chromsome, annotationList.size()));
            writeAnnotations(outputFile, annotationList, true /* append */);
        }
    }

    public static void writeAnnotations(String outputFileName, ObjectList<Annotation> annotationList, boolean append) {
        PrintWriter writer = null;
        File outputFile = new File(outputFileName);

        try {
            if (!outputFile.exists()) {
                writer = new PrintWriter(outputFile);
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
            if (writer != null) {
                writer.close();
            }
        }

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
        new AggregatePeaksByPeakDistanceMode().configure(args).execute();
    }
}