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
import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
import edu.cornell.med.icb.goby.alignments.AlignmentWriter;
import edu.cornell.med.icb.goby.alignments.SortIterateAlignments;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.IOException;

/**
 * Sort an alignment by reference and reference position.
 *
 * @author Fabien Campagne
 */
public class Sort1Mode extends AbstractGobyMode {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(Sort1Mode.class);

    /**
     * The mode name.
     */
    private static final String MODE_NAME = "sort-1";

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
    private long startPosition;
    private long endPosition = Long.MAX_VALUE;


    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }


    public void setInput(final String input) {
        this.basename = AlignmentReaderImpl.getBasename(input);
    }

    public void setOutput(final String output) {
        this.outputFilename = output;
    }

    public void setStartPosition(final long startPosition) {
        this.startPosition = startPosition;
    }

    public void setEndPosition(final long endPosition) {
        this.endPosition = endPosition;
    }

    public void setIncludeReferenceNames(final String includeReferenceNameCommas) {
        alignmentIterator = new SortIterateAlignments();
        alignmentIterator.parseIncludeReferenceArgument(includeReferenceNameCommas);
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
        basename = AlignmentReaderImpl.getBasename(inputFile);
        outputFilename = jsapResult.getString("output");
        alignmentIterator = new SortIterateAlignments();
        alignmentIterator.parseIncludeReferenceArgument(jsapResult);
        endPosition = Long.MAX_VALUE;
        if (jsapResult.contains("start-position") || jsapResult.contains("end-position")) {
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

    /**
     * Sort the alignment.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        final AlignmentReader reader = new AlignmentReaderImpl(basename);
        try {
            reader.readHeader();
            if (reader.isSorted()) {
                LOG.warn("Warning: The input alignment is already sorted.");
            }
        } finally {
            reader.close();
        }

        if (alignmentIterator == null) {
            // If executing from API, this may not have been created yet
            alignmentIterator = new SortIterateAlignments();
        }
        alignmentIterator.setBasename(basename);
        alignmentIterator.setOutputFilename(outputFilename);

        final AlignmentWriter writer = new AlignmentWriter(outputFilename);

        // Iterate through each alignment and write sequence variations to output file:
        LOG.info("Loading entries..");
        alignmentIterator.iterate(startPosition, endPosition, basename);
        LOG.info("Sorting..");
        alignmentIterator.sort();
        LOG.info("Writing sorted alignment..");
        alignmentIterator.writeTmh();
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
        new Sort1Mode().configure(args).execute();
    }
}
