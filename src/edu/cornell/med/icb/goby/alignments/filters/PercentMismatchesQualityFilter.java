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

package edu.cornell.med.icb.alignments.filters;

import edu.cornell.med.icb.alignments.Alignments;
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
