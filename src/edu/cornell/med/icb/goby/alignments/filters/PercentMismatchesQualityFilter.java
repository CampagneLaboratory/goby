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

package edu.cornell.med.icb.goby.alignments.filters;

import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.mssm.crover.cli.CLI;
import org.apache.commons.lang.StringUtils;

import java.io.PrintStream;

/**
 * @author Fabien Campagne
 *         Date: May 6, 2009
 *         Time: 11:11:48 AM
 */
public class PercentMismatchesQualityFilter implements AlignmentQualityFilter {
    private double qualityThresholdPercent;

    /**
     * Reject entries that have more than threshold % differences with the query, in either mismatches or indels.
     *
     * @param header header of the alignment to which this entry belongs.
     * @param entry  The entry to inspect.
     * @return
     */
    public final boolean keepEntry(final Alignments.AlignmentHeader header,
                                   final Alignments.AlignmentEntry entry) {
        final int queryLength = header.getQueryLength(entry.getQueryIndex());
        return keepEntry(queryLength, entry);
    }

    public PercentMismatchesQualityFilter() {
        qualityThresholdPercent = 0.05;
    }

    public final boolean keepEntry(final int queryLength, final Alignments.AlignmentEntry entry) {
        final int numberOfMismatches = entry.getNumberOfMismatches();
        final int numberOfIndels = entry.getNumberOfIndels();
        final int sumDifferences = numberOfIndels + numberOfMismatches;
        return sumDifferences <= (int) Math.round(((float) queryLength) * qualityThresholdPercent);
    }

    public void setParameters(final String parameters) {
        final String[] args = StringUtils.defaultString(parameters).split("[',=]");
        qualityThresholdPercent = CLI.getDoubleOption(args, "threshold", 0.05d);
        System.out.println("Setting quality threshold to " + qualityThresholdPercent);
    }

    public void printUsage(final PrintStream out) {
        out.append("This quality filter rejects alignment entries that have more than a certain threshold " +
                "of differences with the target sequence. Base mismatches, as well as insertion or deletion differences " +
                "are counted towards the difference count. The threshold is set by default to 5% (0.05), " +
                "but can be changed with the threshold parameter. Use syntax threshold=value.\n");

    }
}
