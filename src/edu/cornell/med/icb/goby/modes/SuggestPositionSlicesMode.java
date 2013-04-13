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
import edu.cornell.med.icb.goby.algorithmic.data.Annotation;
import edu.cornell.med.icb.goby.algorithmic.data.ranges.Range;
import edu.cornell.med.icb.goby.algorithmic.data.ranges.Ranges;
import edu.cornell.med.icb.goby.alignments.*;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import it.unimi.dsi.fastutil.objects.*;
import org.apache.commons.io.IOUtils;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Map;

/**
 * Converts a compact alignment to plain text.
 *
 * @author Fabien Campagne
 */
public class SuggestPositionSlicesMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "suggest-position-slices";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Suggest how to slice an alignment by position to yield roughly equally sized slices. ";

    /**
     * The output file.
     */
    private String outputFilename;

    /**
     * The basename of the compact alignment.
     */
    private String[] basenames;
    private int modulo;
    private int numberOfSlices;
    private int bufferLength;
    /**
     * Optional: file name for an annotation file.
     */
    private String annotationFilename;
    private int numBytesPerSlice = -1;
    private boolean useModulo;
    /**
     * When this switch is active, slices that would span two chromosomes are truncated to the end of the first
     * spanned chromosome. The next slice starts at the beginning of the next chromosome.
     */
    private boolean restrictPerChromosome;


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
    public AbstractCommandLineMode configure(final String[] args)
            throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        final String[] inputFiles = jsapResult.getStringArray("input");
        basenames = AlignmentReaderImpl.getBasenames(inputFiles);
        annotationFilename = jsapResult.getString("annotations");
        outputFilename = jsapResult.getString("output");
        modulo = jsapResult.getInt("modulo");
        numberOfSlices = jsapResult.getInt("number-of-slices");
        numBytesPerSlice = jsapResult.getInt("number-of-bytes");
        if (!(jsapResult.userSpecified("number-of-slices") ||
                jsapResult.userSpecified("number-of-bytes"))) {
            System.err.println("You must specify either --number-of-bytes or --number-of-slices");
            System.exit(1);
        }
        restrictPerChromosome = jsapResult.getBoolean("restrict-per-chromosome");
        if (jsapResult.userSpecified("number-of-slices") && jsapResult.userSpecified("number-of-bytes")) {
            System.err.println("You must select either number-of-slices or number-of-bytes, but not both. ");
            System.exit(1);
        }
        useModulo = jsapResult.userSpecified("number-of-slices");
        if (!useModulo) System.err.println("Splitting with " + numBytesPerSlice);
        return this;
    }


    /**
     * Suggests slices to process a large alignment file in parallel.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        PrintStream stream = null;

        try {

            stream = outputFilename == null ? System.out
                    : new PrintStream(new FileOutputStream(outputFilename));


            ConcatSortedAlignmentReader input = new ConcatSortedAlignmentReader(basenames);
            input.readHeader();

            DoubleIndexedIdentifier ids = new DoubleIndexedIdentifier(input.getTargetIdentifiers());

            Ranges ranges = null;
            if (annotationFilename != null) {
                ranges = convertAnnotationsToRanges(annotationFilename, ids, input);
            }

            ReferenceLocation[] breakpoints = useModulo ? getReferenceLocationsWithModulo(stream, input, ids) :
                    getReferenceLocationsWithBytes(stream, input, ids);

            if (ranges != null) {
                adjustBreakpointsWithAnnotations(breakpoints, ranges);
            }
            if (restrictPerChromosome) {
                breakpoints = restrictPerChromosome(breakpoints, input);
            }
            for (int i = 0; i < numberOfSlices; i++) {
                if (!restrictPerChromosome ||
                        (restrictPerChromosome && breakpoints[i].targetIndex == breakpoints[i + 1].targetIndex))
                stream.printf(String.format("%s\t%d\t%s,%d\t%s\t%d\t%s,%d%n",
                        ids.getId(breakpoints[i].targetIndex),
                        breakpoints[i].position,
                        ids.getId(breakpoints[i].targetIndex),
                        breakpoints[i].position,
                        ids.getId(breakpoints[i + 1].targetIndex),
                        breakpoints[i + 1].position,
                        ids.getId(breakpoints[i + 1].targetIndex),
                        breakpoints[i + 1].position));

            }

        } finally {
            if (stream != System.out) {
                IOUtils.closeQuietly(stream);
            }
        }
    }

    private ReferenceLocation[] restrictPerChromosome(ReferenceLocation[] breakpoints, ConcatSortedAlignmentReader reader) {
        ObjectArrayList<ReferenceLocation> result = new ObjectArrayList<ReferenceLocation>();
        int lastTargetIndex = -1;
        int index = 0;
        for (ReferenceLocation breakpoint : breakpoints) {
            if (breakpoint.targetIndex != lastTargetIndex) {
                // we switch to a new chromosome, introduce a new breakpoint at the end of the previous chromosome:
                if (lastTargetIndex != -1) {
                    result.add(new ReferenceLocation(lastTargetIndex, reader.getTargetLength(lastTargetIndex)));
                    result.add(new ReferenceLocation(breakpoint.targetIndex, 0));
                    System.out.println("Adding breakpoint at end of "+lastTargetIndex);
                }
                lastTargetIndex = breakpoint.targetIndex;
            }
            result.add(breakpoint);
            index++;
        }
        return result.toArray(new ReferenceLocation[result.size()]);
    }

    private ReferenceLocation[] getReferenceLocationsWithBytes(PrintStream stream, ConcatSortedAlignmentReader input, DoubleIndexedIdentifier ids) throws IOException {
        ObjectList<ReferenceLocation> locations = input.getLocationsByBytes(numBytesPerSlice);
        numberOfSlices = locations.size();
        return prepareBreakpoints(stream, input, ids, locations);
    }

    private ReferenceLocation[] getReferenceLocationsWithModulo(PrintStream stream, ConcatSortedAlignmentReader input, DoubleIndexedIdentifier ids) throws IOException {
        ObjectList<ReferenceLocation> locations = input.getLocations(modulo);

        if (locations.size() < numberOfSlices) {
            numberOfSlices = locations.size();
        }

        return prepareBreakpoints(stream, input, ids, locations);
    }

    private ReferenceLocation[] prepareBreakpoints(PrintStream stream, ConcatSortedAlignmentReader input, DoubleIndexedIdentifier ids, ObjectList<ReferenceLocation> locations) {
        ReferenceLocation first;
        ReferenceLocation[] breakpoints = new ReferenceLocation[numberOfSlices + 1];
        breakpoints[0] = first = locations.get(0);
        locations.remove(first);
        stream.println("targetIdStart\t%positionStart\tstart:(ref,pos)\ttargetIdEnd\t%positionEnd\tend:(ref,pos)");
        for (int i = 0; i < numberOfSlices - 1; i++) {
            breakpoints[i + 1] = locations.get(locations.size() / (numberOfSlices - 1) * i);
        }

        // largest position in the last reference sequence:
        final int lastTargetIndex = ids.size() - 1;
        breakpoints[breakpoints.length - 1] = new ReferenceLocation(lastTargetIndex, input.getTargetLength(lastTargetIndex));
        return breakpoints;
    }

    private void adjustBreakpointsWithAnnotations(final ReferenceLocation[] breakpoints, final Ranges ranges) {
        ranges.order();
        for (final ReferenceLocation breakpoint : breakpoints) {
            int referenceIndex = breakpoint.targetIndex;
            int position = breakpoint.position;
            Range nonOverlapping = ranges.findNonOverlappingRange(referenceIndex, position);

            // change the breakpoint to a position in the middle of the non-overlapping segment:
            breakpoint.position = (nonOverlapping.min + nonOverlapping.max) / 2;

        }
    }

    private Ranges convertAnnotationsToRanges(String annotationFilename, DoubleIndexedIdentifier ids, ConcatSortedAlignmentReader input) throws IOException {
        Ranges ranges = new Ranges();
        if (annotationFilename == null) return null;
        Object2ObjectMap<String, ObjectList<Annotation>> annotations = CompactAlignmentToAnnotationCountsMode.readAnnotations(annotationFilename);
        Object2IntMap<String> starts = new Object2IntLinkedOpenHashMap<String>();
        Object2IntMap<String> ends = new Object2IntLinkedOpenHashMap<String>();
        Object2IntMap<String> refIndices = new Object2IntLinkedOpenHashMap<String>();
        starts.defaultReturnValue(Integer.MAX_VALUE);
        ends.defaultReturnValue(Integer.MIN_VALUE);

        for (ObjectList<Annotation> value : annotations.values()) {

            for (Annotation ann : value) {

                int min = Integer.MAX_VALUE;
                int max = Integer.MIN_VALUE;
                String chromosome = ann.getChromosome();
                String id = ann.getId();
                min = starts.getInt(id);
                max = ends.getInt(id);

                starts.put(id, Math.min(min, ann.getStart()));
                ends.put(id, Math.max(max, ann.getEnd()));

                assert chromosome != null : " annotation must contain non null chromosome.";

            }
        }
        for (final String id : starts.keySet()) {
            final Range range = new Range();
            final int referenceIndex = refIndices.getInt(id);
            range.min = starts.getInt(id);
            range.max = ends.getInt(id);
            ranges.add(range, referenceIndex);
        }

        ranges.order();
        return ranges;
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
        new SuggestPositionSlicesMode().configure(args).execute();
    }
}