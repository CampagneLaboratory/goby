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
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.apache.commons.io.IOUtils;

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
    private static final String MODE_NAME = "count-archive-to-peak-annotations";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Write annotations corresponding to consensus "
            + "peaks found in each sequence of count archives.";

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
     * @throws java.io.IOException error parsing
     * @throws com.martiansoftware.jsap.JSAPException
     *                             error parsing
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
        /**
         * TODO: Determine of adjustQueryIndices should be the default of true.
         */
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

        // for each reference sequence, generate annotations for the union of peaks across all the input samples.
        // More precisely, we generate an annotation across the widest peak that can be called with the sum of counts
        // across all the input samples.
        final Object2ObjectMap<String, ObjectList<Annotation>> allAnnots = new Object2ObjectOpenHashMap<String, ObjectList<Annotation>>();
        PrintWriter annotationWriter = null;
        try {
            annotationWriter = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));
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
                            if (allAnnots.get(referenceId) == null) {
                                allAnnots.put(referenceId, new ObjectArrayList<Annotation>());
                            }
                        }
                    }
                    final AnyTransitionCountsIterator iterator = new AnyTransitionCountsIterator(readers);

                    final PeakAggregator peakAggregator = new PeakAggregator(iterator);
                    peakAggregator.setPeakDetectionThreshold(detectionThreshold);
                    while (peakAggregator.hasNext()) {
                        final Peak peak = peakAggregator.next();

                        final ObjectList<Annotation> annotationList = allAnnots.get(referenceId);
                        final String id = referenceId + "." + peak.start + "." + peak.length;
                        final Annotation annotation = new Annotation(id, referenceId, "either");
                        annotation.addSegment(new Segment(peak.start, peak.start + peak.length, "id", "either"));
                        annotationList.add(annotation);
                        annotation.write(annotationWriter);
                    }
                }
            }
        } finally {
            IOUtils.closeQuietly(annotationWriter);
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
        new CountArchiveToPeakAnnotationsMode().configure(args).execute();
    }
}
