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

package edu.cornell.med.icb.goby.alignments.filters;

import edu.cornell.med.icb.goby.alignments.Alignments;

import java.io.PrintStream;

/**
 * Quality filters make decision about whether alignment entries are of sufficient quality
 * to be kept.
 *
 * @author Fabien Campagne
 *         Date: May 6, 2009
 *         Time: 11:08:48 AM
 */
public interface AlignmentQualityFilter {
    /**
     * Returns true if the entry passes the quality criteria of this filter.
     *
     * @param header header of the alignment to which this entry belongs.
     * @param entry  The entry to inspect.
     * @return True or false.
     */
    boolean keepEntry(Alignments.AlignmentHeader header, Alignments.AlignmentEntry entry);
     /**
     * Returns true if the entry passes the quality criteria of this filter.
     *
     * @param queryLength the length of the query sequence described by this entry.
     * @param entry  The entry to inspect.
     * @return True or false.
     */
    boolean keepEntry(int queryLength, Alignments.AlignmentEntry entry);

    /**
     * Set parameters for this filter.
     * @param parameters A string in the format described by printUsage.
     */
    public void setParameters(String parameters);

    /**
     * Print usage information for this filter. Describes any optional parameters to the end user.
     * @param out
     */
    public void printUsage(PrintStream out);

}
