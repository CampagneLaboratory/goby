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
import org.apache.commons.io.IOUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.IOException;
import java.io.PrintStream;
import java.io.FileOutputStream;
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
    private static final String MODE_DESCRIPTION = "Sort a compact alignment by reference position.";

    /**
     * The output file.
     */
    private String outputFilename;

    /**
     * The basename of the compact alignment.
     */
    private String[] basenames;
    private MyIterateAlignments alignmentIterator;


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
        basenames = AlignmentReader.getBasenames(inputFiles);
        outputFilename = jsapResult.getString("output");
        alignmentIterator = new MyIterateAlignments();
        alignmentIterator.parseIncludeReferenceArgument(jsapResult);

        return this;
    }

    private class MyIterateAlignments extends IterateAlignments {


        private MyIterateAlignments() {

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
            Merge.prepareMergedTooManyHits(outputFilename, alignmentReader.getNumberOfQueries(), 0, basenames);
            writer.setTargetIdentifiers(alignmentReader.getTargetIdentifiers());
            writer.setQueryIdentifiers(alignmentReader.getQueryIdentifiers());


            for (Alignments.AlignmentEntry entry : entries) {
                writer.appendEntry(entry);
            }

            writer.close();
        }
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
        PrintStream stream = null;
        try {
            stream = outputFilename == null ? System.out
                    : new PrintStream(new FileOutputStream(outputFilename));
            AlignmentWriter writer = new AlignmentWriter(outputFilename);
            writer.setSorted(true);

            // Iterate through each alignment and write sequence variations to output file:
            LOG.info("Loading entries..");
            alignmentIterator.iterate(basenames);
            LOG.info("Sorting..");
            alignmentIterator.sort();
            LOG.info("Writing sorted alignment..");
            alignmentIterator.write(writer);

        } finally {
            if (stream != System.out) {
                IOUtils.closeQuietly(stream);
            }
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
        new SortMode().configure(args).execute();
    }

}