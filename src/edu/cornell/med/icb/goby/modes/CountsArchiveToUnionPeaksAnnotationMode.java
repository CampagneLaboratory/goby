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
import edu.cornell.med.icb.goby.algorithmic.data.Annotation;
import edu.cornell.med.icb.goby.algorithmic.data.Segment;
import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.ConcatAlignmentReader;
import edu.cornell.med.icb.goby.counts.AnyTransitionCountsIterator;
import edu.cornell.med.icb.goby.counts.CountsArchiveReader;
import edu.cornell.med.icb.goby.counts.CountsReaderI;
import edu.cornell.med.icb.goby.counts.Peak;
import edu.cornell.med.icb.goby.counts.PeakAggregator;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.objects.Object2ObjectMap;
import it.unimi.dsi.fastutil.objects.Object2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.fastutil.objects.ObjectListIterator;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.apache.commons.io.IOUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collections;
import java.util.Set;


/**
 * Write annotations corresponding to a consolidation of peaks accross each sequence in
 * the count archives.
 *
 * @author Jaaved Mohammed
 */
public class CountsArchiveToUnionPeaksAnnotationMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "count-archive-to-peak-union-annotations";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Write annotations corresponding "
            + "to a consolidation of peaks accross each sequence in the count archives.";

    /**
     * The input files. Must reduce to alignment basenames.
     */
    private String[] inputFiles;

    /**
     * The output file.
     */
    private String outputFile;

    private boolean filterByReferenceNames;
    private ObjectSet<String> includeReferenceNames;
    private String alternativeCountsName;
    private int detectionThreshold;

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
     * @throws IOException error parsing
     * @throws JSAPException error parsing
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args) throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        inputFiles = AlignmentReader.getBasenames(jsapResult.getStringArray("input"));
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
        alternativeCountsName = jsapResult.getString("alternative-count-archive");
        detectionThreshold = jsapResult.getInt("threshold");
        System.out.println("Peak Detection Threshold is " + detectionThreshold);

        return this;
    }

    /**
     * Run the mode.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        final String[] basenames = AlignmentReader.getBasenames(inputFiles);
        final ConcatAlignmentReader reader = new ConcatAlignmentReader(basenames);
        reader.readHeader();
        final int numberOfReferences = reader.getNumberOfTargets();
        final IndexedIdentifier referenceIds = reader.getTargetIdentifiers();
        final DoubleIndexedIdentifier backwards = new DoubleIndexedIdentifier(referenceIds);
        reader.close();

        // read counts archive:
        final CountsArchiveReader[] countArchiveReaders = new CountsArchiveReader[inputFiles.length];
        int i = 0;
        for (final String inputFile : basenames) {
            countArchiveReaders[i++] = new CountsArchiveReader(inputFile, alternativeCountsName);
        }

        // for each reference sequence, generate annotations for the union of peaks across all
        //  the input samples. More precisely, we generate an annotation across the widest peak
        // that can be called with the sum of counts across all the input samples.
        final Object2ObjectMap<String, ObjectList<Segment>> allSegments = new Object2ObjectOpenHashMap<String, ObjectList<Segment>>();
        for (int referenceIndex = 0; referenceIndex < numberOfReferences; referenceIndex++) {
            final String referenceId = backwards.getId(referenceIndex).toString();
            boolean processThisSequence = true;

            if (filterByReferenceNames && !includeReferenceNames.contains(referenceId)) {
                processThisSequence = false;
            }

            if (processThisSequence) {
                System.out.println("Processing reference " + referenceId);
                final CountsReaderI[] readers = new CountsReaderI[countArchiveReaders.length];
                int readerIndex = 0;
                for (final CountsArchiveReader archive : countArchiveReaders) {
                    if (archive.getIdentifier(referenceIndex) != null) {
                        readers[readerIndex++] = archive.getCountReader(referenceIndex);

                        if (allSegments.get(referenceId) == null) {
                            allSegments.put(referenceId, new ObjectArrayList<Segment>());
                        }
                    }
                }

                // Reads in all the files and defines one iterator over all input files
                final AnyTransitionCountsIterator iterator = new AnyTransitionCountsIterator(readers);
                // Given all input files and one iterator over them, start collecting (possibly overlapping)
                // peaks across all input files
                final PeakAggregator peakAggregator = new PeakAggregator(iterator);
                peakAggregator.setPeakDetectionThreshold(detectionThreshold);
                while (peakAggregator.hasNext()) {
                    final Peak peak = peakAggregator.next();

                    final ObjectList<Segment> segmentsList = allSegments.get(referenceId);
                    segmentsList.add(new Segment(peak.start, peak.start + peak.length, "id", "either"));
                }
            }
        }

        // Consolidate all overlapping peaks across all input files by taking the
        // union of overlapping peaks
        System.out.println("Consolidating overlapping peaks over all counts files.");
        final Set<String> chromsomes = allSegments.keySet();
        ObjectList<Segment> consensusList = null;
        for (final String chromsome : chromsomes) {
            consensusList = new ObjectArrayList<Segment>();
            final ObjectList<Segment> segmentsList = allSegments.get(chromsome);

            // Sort the segments to make consolidation by union easier.
            Collections.sort(segmentsList);
            System.out.println(String.format("Reference %s: Starting with %d segments.", chromsome, segmentsList.size()));
            ObjectListIterator<Segment> segIterator = segmentsList.listIterator();

            Segment segment = null;
            Segment nextSeg = null;
            Segment consolidatedSeg = null;
            if (segIterator.hasNext()) {
                segment = segIterator.next();
            }
            while (segment != null) {
                if (segIterator.hasNext()) {
                    nextSeg = segIterator.next();
                    while ((nextSeg != null) && segment.overlap(nextSeg)) {
                        // merge the two segments
                        consolidatedSeg = segment.merge(nextSeg);
                        segment = consolidatedSeg;
                        if (segIterator.hasNext()) {
                            nextSeg = segIterator.next();
                        } else {
                            nextSeg = null;
                        }
                    }
                }
                consensusList.add(segment);
                segment = nextSeg;
                nextSeg = null;
            }
            Collections.sort(consensusList);
            segmentsList.clear();
            segmentsList.addAll(consensusList);

            System.out.println(String.format("Reference %s: Finished with %d segments.", chromsome, segmentsList.size()));

            // Transform the segments back to annotations for easy writing to file
            final ObjectList<Annotation> annotationList = new ObjectArrayList<Annotation>();
            segIterator = segmentsList.listIterator();
            while (segIterator.hasNext()) {
                segment = segIterator.next();
                final Annotation annot = new Annotation(chromsome + "." + segment.getStart() + "." + segment.getLength(), chromsome);
                annot.addSegment(segment);
                annotationList.add(annot);
            }

            // Write the list of annotations to file
            writeAnnotations(outputFile, annotationList, true);
        }
    }

    public static void writeAnnotations(final String outputFileName, final ObjectList<Annotation> annotationList, final boolean append) {
        final File outputFile = new File(outputFileName);
        PrintWriter writer = null;
        try {
            if (!outputFile.exists()) {
                writer = new PrintWriter(outputFile);
                writer.write("Chromosome_Name\tStrand\tPrimary_ID\tSecondary_ID\tTranscript_Start\tTranscript_End\n");
            } else {
                writer = new PrintWriter(new FileOutputStream(outputFile, append));
            }

            final ObjectListIterator<Annotation> annotIterator = annotationList.listIterator();
            while (annotIterator.hasNext()) {
                final Annotation annotation = annotIterator.next();
                annotation.write(writer);
            }
        } catch (FileNotFoundException fnfe) {
            System.err.println("Caught exception in writeAnnotations: " + fnfe.getMessage());
            System.exit(1);
        } finally {
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
        new CountsArchiveToUnionPeaksAnnotationMode().configure(args).execute();
    }
}
