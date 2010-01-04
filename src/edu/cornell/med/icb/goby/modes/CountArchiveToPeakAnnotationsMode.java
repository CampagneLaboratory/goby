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

import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.alignments.AlignmentReader;
import edu.cornell.med.icb.alignments.ConcatAlignmentReader;
import edu.cornell.med.icb.counts.AnyTransitionCountsIterator;
import edu.cornell.med.icb.counts.CountsArchiveReader;
import edu.cornell.med.icb.counts.CountsReaderI;
import edu.cornell.med.icb.counts.Peak;
import edu.cornell.med.icb.counts.PeakAggregator;
import edu.cornell.med.icb.goby.algorithmic.data.Annotation;
import edu.cornell.med.icb.goby.algorithmic.data.Segment;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.Object2ObjectMap;
import it.unimi.dsi.fastutil.objects.Object2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;


/**
 * Write annotations corresponding to consensus peaks found in each sequence of count archives.
 *
 * @author Fabien Campagne
 */
public class CountArchiveToPeakAnnotationsMode extends AbstractGobyMode {

    /**
     * The mode name.
     */
    public static final String MODE_NAME = "count-archive-to-peak-annotations";
    public static final String MODE_DESCRIPTION = "Write annotations corresponding to consensus peaks found in each sequence of count archives.";


    /**
     * The input files. Must reduce to alignment basenames.
     */
    private String[] inputFiles;

    /**
     * The output file.
     */
    private String outputFile;

    boolean filterByReferenceNames;
    private ObjectSet<String> includeReferenceNames;
    private String alternativeCountsName;
    private int detectionThreshold;

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

        inputFiles = AlignmentReader.getBasenames(jsapResult.getStringArray("input"));
        outputFile = jsapResult.getString("output");
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
        alternativeCountsName = jsapResult.getString("alternative-count-archive");
        detectionThreshold=jsapResult.getInt("threshold");
        return this;
    }

    /**
     * Run the mode.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        final IntSet referencesToProcess = new IntOpenHashSet();

        final String[] basenames = AlignmentReader.getBasenames(inputFiles);
        final ConcatAlignmentReader reader = new ConcatAlignmentReader(basenames);
        reader.readHeader();
        final int numberOfReferences = reader.getNumberOfTargets();
        final IndexedIdentifier referenceIds = reader.getTargetIdentifiers();
        final DoubleIndexedIdentifier backwards = new DoubleIndexedIdentifier(referenceIds);
        reader.close();

        // read counts archive:
        CountsArchiveReader countArchiveReaders[] = new CountsArchiveReader[inputFiles.length];
        int i = 0;
        for (String inputFile : basenames) {
            countArchiveReaders[i++] = new CountsArchiveReader(inputFile, alternativeCountsName);
        }

        // for each reference sequence, generate annotations for the union of peaks across all the input samples.
        // More precisely, we generate an annotation across the widest peak that can be called with the sum of counts
        // across all the input samples.
        final Object2ObjectMap<String, ObjectList<Annotation>> allAnnots = new Object2ObjectOpenHashMap<String, ObjectList<Annotation>>();
        PrintWriter annotationWriter = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));
        for (int referenceIndex = 0; referenceIndex < numberOfReferences; referenceIndex++) {
            String referenceId = backwards.getId(referenceIndex).toString();
            System.out.println("Processing reference " + referenceId);
            CountsReaderI[] readers = new CountsReaderI[countArchiveReaders.length];
            int readerIndex = 0;
            for (CountsArchiveReader archive : countArchiveReaders) {

                if (archive.getIdentifier(referenceIndex) != null) {

                    readers[readerIndex++] = archive.getCountReader(referenceIndex);

                    if (allAnnots.get(referenceId) == null) {
                        allAnnots.put(referenceId, new ObjectArrayList<Annotation>());
                    }
                }
            }
            AnyTransitionCountsIterator iterator = new AnyTransitionCountsIterator(readers);

            PeakAggregator peakAggregator = new PeakAggregator(iterator);
            peakAggregator.setPeakDetectionThreshold(detectionThreshold);
            while (peakAggregator.hasNext()) {
                Peak peak = peakAggregator.next();

                ObjectList<Annotation> annotationList = allAnnots.get(referenceId);
                final Annotation annotation = new Annotation(referenceId + "." + peak.start + "." + peak.length, referenceId);
                annotation.strand = "either";
                annotation.addSegment(new Segment(peak.start, peak.start + peak.length, "id", "either"));
                annotationList.add(annotation);
                annotation.write(annotationWriter);


            }
        }
        annotationWriter.close();
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
        new CountArchiveToPeakAnnotationsMode().configure(args).execute();
    }
}
