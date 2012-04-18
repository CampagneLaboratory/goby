/*
 * Copyright (C) 2009-2012 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.alignments;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;

public class SortIterateAlignments extends IterateAlignments {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(SortIterateAlignments.class);

    ObjectArrayList<Alignments.AlignmentEntry> entries;

    String outputFilename;
    String basename;


    private static final Comparator<Alignments.AlignmentEntry> POSITION_ENTRY_COMPARATOR =
            new AlignmentPositionComparator();

    public SortIterateAlignments() {
        entries = new ObjectArrayList<Alignments.AlignmentEntry>();
    }

    public String getBasename() {
        return basename;
    }

    public void setBasename(final String basename) {
        this.basename = basename;
    }

    public String getOutputFilename() {
        return outputFilename;
    }

    public void setOutputFilename(final String outputFilename) {
        this.outputFilename = outputFilename;
    }

    @Override
    public void prepareDataStructuresForReference(final AlignmentReader alignmentReader, final int referenceIndex) {
        if (this.alignmentReader == null) {
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

    public void writeTmh() throws IOException {
        Merge.prepareMergedTooManyHits(outputFilename, alignmentReader.getNumberOfQueries(), 0, basename);
    }

    public void write(final AlignmentWriterImpl writer) throws IOException {
        // too many hits is prepared as for Merge:
        try {
            writer.setTargetIdentifiers(alignmentReader.getTargetIdentifiers());
            writer.setQueryIdentifiers(alignmentReader.getQueryIdentifiers());
            final int[] targetLengths = alignmentReader.getTargetLength();
            if (targetLengths != null) {
                writer.setTargetLengths(targetLengths);
            }
            writer.setLargestSplitQueryIndex(alignmentReader.getLargestSplitQueryIndex());
            writer.setSmallestSplitQueryIndex(alignmentReader.getSmallestSplitQueryIndex());
            writer.setSorted(true);
            writer.setAlignerName(alignmentReader.getAlignerName());
            writer.setAlignerVersion(alignmentReader.getAlignerVersion());

            // Propagate the statistics from the input, but update the basename
            writer.setStatistics(alignmentReader.getStatistics());
            writer.putStatistic("basename", FilenameUtils.getName(basename));
            writer.putStatistic("basename.full", basename);

            for (final Alignments.AlignmentEntry entry : entries) {
                writer.appendEntry(entry);
            }
        } finally {
            writer.close();
        }

    }
}
