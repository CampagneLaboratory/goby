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
import edu.cornell.med.icb.goby.alignments.*;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;

/**
 * Sort an alignment by reference and reference position.
 *
 * @author Fabien Campagne
 */
public class SortMode extends AbstractGobyMode {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(CompactAlignmentToAnnotationCountsMode.class);

    /**
     * The mode name.
     */
    private static final String MODE_NAME = "sort";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Sort a compact alignment by reference position. The output alignment is sorted and indexed.";

    /**
     * The output file.
     */
    private String outputFilename;

    /**
     * The basename of the compact alignment.
     */
    private String basename;
    private SortIterateAlignments alignmentIterator;
    private int[] targetLengths;
    private boolean hasStartOrEndPosition;
    private long startPosition;
    private long endPosition;


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

        final String inputFile = jsapResult.getString("input");
        basename = AlignmentReader.getBasename(inputFile);
        outputFilename = jsapResult.getString("output");
        alignmentIterator = new SortIterateAlignments();
        alignmentIterator.parseIncludeReferenceArgument(jsapResult);
        if (jsapResult.contains("start-position") || jsapResult.contains("end-position")) {
            hasStartOrEndPosition = true;
            startPosition = jsapResult.getLong("start-position", 0L);
            endPosition = jsapResult.getLong("end-position", Long.MAX_VALUE);
        }

        if (startPosition < 0L) {
            throw new JSAPException("Start position must not be less than zero");
        }
        if (endPosition < 0L) {
            throw new JSAPException("End position must not be less than zero");
        }
        if (startPosition > endPosition) {
            throw new JSAPException("Start position must not be greater than the end position");
        }

        return this;
    }


    private static final Comparator<Alignments.AlignmentEntry> POSITION_ENTRY_COMPARATOR =
            new AlignmentPositionComparator();

    /**
     * Sort the alignment.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {


        AlignmentWriter writer = new AlignmentWriter(outputFilename);

        // Iterate through each alignment and write sequence variations to output file:
        LOG.info("Loading entries..");
        alignmentIterator.iterate(startPosition, endPosition, basename);
        LOG.info("Sorting..");
        alignmentIterator.sort();
        LOG.info("Writing sorted alignment..");
        alignmentIterator.write(writer);

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
        new SortMode().configure(args).execute();
    }

    private class SortIterateAlignments extends IterateAlignments {


        private SortIterateAlignments() {

            entries = new ObjectArrayList<Alignments.AlignmentEntry>();
        }


        ObjectArrayList<Alignments.AlignmentEntry> entries;

        @Override
        public void prepareDataStructuresForReference(AlignmentReader alignmentReader, int referenceIndex) {
            if (this.alignmentReader == null) {
                if (alignmentReader.isSorted()) {
                    LOG.warn("Warning: An input alignment is already sorted.");
                }
                this.alignmentReader = alignmentReader;
            }

        }

        @Override
        public void processAlignmentEntry(final AlignmentReader alignmentReader,
                                          final Alignments.AlignmentEntry alignmentEntry) {

            entries.add(alignmentEntry);
        }

        AlignmentReader alignmentReader = null;

        public void sort() {
            Collections.sort(entries, POSITION_ENTRY_COMPARATOR);
        }

        public void write(AlignmentWriter writer) throws IOException {
            // too many hits is prepared as for Merge:
            try {
                Merge.prepareMergedTooManyHits(outputFilename, alignmentReader.getNumberOfQueries(), 0, basename);
                writer.setTargetIdentifiers(alignmentReader.getTargetIdentifiers());
                writer.setQueryIdentifiers(alignmentReader.getQueryIdentifiers());
                writer.setTargetLengths(alignmentReader.getTargetLength());
                writer.setSorted(true);

                for (Alignments.AlignmentEntry entry : entries) {
                    writer.appendEntry(entry);
                }
            } finally {
                writer.close();
            }

        }
    }
}