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

import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.algorithmic.data.Segment;
import edu.cornell.med.icb.goby.algorithmic.data.Annotation;

import it.unimi.dsi.fastutil.objects.*;

import java.io.*;
import java.util.*;

import org.apache.commons.io.IOUtils;


/**
 * Write annotations corresponding to consensus peaks found in each sequence of count archives.
 *
 * @author Jaaved Mohammed
 */
public class AnnotationPenaltyMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    public static final String MODE_NAME = "annotation-penalty";
    public static final String MODE_DESCRIPTION = "Calculates the distances between two collections of annotations.";

     /**
     * Constants
      * These penalty scores are designed to be small enough for whole genome sequences (mouse)
      * They may need to be further shrunk to reduce overflow error for larger genomes.
     */
    public static final double NON_OVERLAP_COST = 0.02;
    public static final double OVERLAP_COST = 0.01;


    /**
     * The input files. Must reduce to alignment basenames.
     */
    private String[] inputFiles;

    /**
     * The output file.
     */
    private String outputFile;
    private String geneAnnotFile;
    private String proposalAnnotFile;

    private static double score; //The total penalty score, resets to zero at each execution of the mode
    private static int lastPos; //Keeps track of the last position processed until to prevent duplication.

    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }

    public double getScore() {
        return score;
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

        outputFile = jsapResult.getString("output");
        geneAnnotFile = jsapResult.getString("gene-annotation");
        proposalAnnotFile = jsapResult.getString("proposal-annotation");

        return this;
    }

    /**
     * Run the mode.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        score = 0d;
        lastPos = -1;


        System.out.println("Reading official gene annotations from " + geneAnnotFile);
        final Object2ObjectMap<String, ObjectList<Annotation>> allGeneAnnots = readAnnotations(geneAnnotFile);
        System.out.println("Reading official proposal annotations from " + proposalAnnotFile);
        final Object2ObjectMap<String, ObjectList<Annotation>> allProposalAnnots = readAnnotations(proposalAnnotFile);

        //Find the references (chromosomes) common in both annotation files
        Set<String> allGeneChroms = allGeneAnnots.keySet();
        Set<String> allProposalChroms = allProposalAnnots.keySet();
        Set<String> intersect = allGeneChroms;
        intersect.retainAll(allProposalChroms);
        Iterator<String> commonChromsIt = intersect.iterator();

        //Find the overlapping differences between the official gene annotation and the proposed one.
        while (commonChromsIt.hasNext()) {
            lastPos = 0;
            final String chromosome = commonChromsIt.next();

            final ObjectList<Annotation> geneAnnotations = allGeneAnnots.get(chromosome);
            final ObjectList<Annotation> proposalAnnotations = allProposalAnnots.get(chromosome);
            //Sort to allow fast comparison
            Collections.sort(geneAnnotations);
            Collections.sort(proposalAnnotations);

            final ObjectIterator<Annotation> geneAnnotIt = geneAnnotations.iterator();
            final ObjectIterator<Annotation> propAnnotIt = proposalAnnotations.iterator();

            Annotation geneAnnot = getNextAnnotation(geneAnnotIt);
            Annotation propAnnot = getNextAnnotation(propAnnotIt);

            while ((geneAnnot != null) || (propAnnot != null)) {
                if ((geneAnnot != null) && (propAnnot != null)) {
                    //check to see if they overlap
                    if (!geneAnnot.overlap(propAnnot)) {
                        //System.out.println("\tDon't overlap: geneAnnot = " + geneAnnot.getStart() + "-" + geneAnnot.getEnd() + ", propAnnot = " + propAnnot.getStart() + "-" + propAnnot.getEnd());
                        //determine which iterator to increment
                        if (geneAnnot.getStart() > propAnnot.getStart()) {
                            addPenalty(propAnnot);
                            propAnnot = getNextAnnotation(propAnnotIt);
                        } else {
                            addPenalty(geneAnnot);
                            geneAnnot = getNextAnnotation(geneAnnotIt);
                        }
                    } else {
                        //They do overlap
                        //Need to determine which annotation is the shorter sequence and advance it's iterator
                        //System.out.println("\t\tDO overlap: geneAnnot = " + geneAnnot.getStart() + "-" + geneAnnot.getEnd() + ", propAnnot = " + propAnnot.getStart() + "-" + propAnnot.getEnd());
                        addPenalty(geneAnnot, propAnnot);
                        if (geneAnnot.getEnd() < propAnnot.getEnd()) {
                            geneAnnot = getNextAnnotation(geneAnnotIt);
                        } else {
                            propAnnot = getNextAnnotation(propAnnotIt);
                        }
                    }
                } else if (geneAnnot == null) {
                    //System.out.println("No more geneAnnot remaining");
                    //simply add the cost of all remaining propAnnot
                    addPenalty(propAnnot);
                    propAnnot = getNextAnnotation(propAnnotIt);
                } else if (propAnnot == null) {
                    //System.out.println("No more propAnnot remaining");
                    //simply add the cost of all remaining geneAnnot
                    addPenalty(geneAnnot);
                    geneAnnot = getNextAnnotation(geneAnnotIt);
                }
            }
        }
        System.out.println("Score = " + score);
    }

    private static void addPenalty(Annotation annot) {
        score += (annot.getLength() * NON_OVERLAP_COST);
        lastPos = annot.getEnd();
    }

    private static void addPenalty(Segment seg) {
        if(seg.end < lastPos)
            return;

        score += (seg.getLength() * NON_OVERLAP_COST);
        lastPos = seg.end;
    }

    private static void addPenalty(Segment geneSeg, Segment propSeg) {
        int stop = Math.min(geneSeg.end, propSeg.end);

        if(stop < lastPos)
            return;

        int overlapLen = 0, nonOverlapLen = 0;
        for(int i = lastPos; i < stop; i++)
        {
            if((geneSeg.start <= i && i <= geneSeg.end) && (propSeg.start <= i && i <= propSeg.end))
            {
                overlapLen++;
            }
            else
            {
                nonOverlapLen++;
            }
        }

        score += ((nonOverlapLen * NON_OVERLAP_COST) + (overlapLen * OVERLAP_COST));
        lastPos = stop;
    }

    private static void addPenalty(Annotation geneAnnot, Annotation propAnnot) {

        final ObjectIterator<Segment> geneSegIt = geneAnnot.segments.iterator();
        final ObjectIterator<Segment> propSegIt = propAnnot.segments.iterator();
        Segment geneSeg = getNextSegment(geneSegIt);
        Segment propSeg = getNextSegment(propSegIt);

            while ((geneSeg != null) || (propSeg != null)) {
                if ((geneSeg != null) && (propSeg != null)) {
                    //check to see if they overlap
                    if (!geneSeg.overlap(propSeg)) {
                        //System.out.println("\t\t\tDon't overlap: geneSeg = " + geneSeg.start + "-" + geneSeg.end + ", propSeg = " + propSeg.start + "-" + propSeg.end);
                        //determine which iterator to increment
                        if (geneSeg.start > propSeg.start) {
                            addPenalty(propSeg);
                            propSeg = getNextSegment(propSegIt);
                        } else {
                            addPenalty(geneSeg);
                            geneSeg = getNextSegment(geneSegIt);
                        }
                    } else {
                        //They do overlap
                        //Need to determine which annotation is the shorter sequence and advance it's iterator
                        //System.out.println("\t\t\t\tDO overlap: geneSeg = " + geneSeg.start + "-" + geneSeg.end + ", propSeg = " + propSeg.start + "-" + propSeg.end);
                        addPenalty(geneSeg, propSeg);
                        if (geneSeg.end < propSeg.end) {
                            geneSeg = getNextSegment(geneSegIt);
                        } else {
                            propSeg = getNextSegment(propSegIt);
                        }
                    }
                } else if (geneSeg == null) {
                    //System.out.println("\t\t\tNo more geneSeg remaining");
                    //simply add the cost of all remaining propSeg
                    addPenalty(propSeg);
                    propSeg = getNextSegment(propSegIt);
                } else if (propSeg == null) {
                    //System.out.println("\t\t\tNo more propSeg remaining");
                    //simply add the cost of all remaining geneSeg
                    addPenalty(geneSeg);
                    geneSeg = getNextSegment(geneSegIt);
                }
            }
    }

    private static Annotation getNextAnnotation(ObjectIterator<Annotation> iterator) {
        Annotation annotation = null;
        if (iterator.hasNext()) {
            annotation = iterator.next();
            annotation.sortSegments();
        }
        return annotation;
    }

        private static Segment getNextSegment(ObjectIterator<Segment> iterator) {
        Segment segment = null;
        if (iterator.hasNext()) {
            segment = iterator.next();
        }
        return segment;
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
        new AnnotationPenaltyMode().configure(args).execute();
    }
}
