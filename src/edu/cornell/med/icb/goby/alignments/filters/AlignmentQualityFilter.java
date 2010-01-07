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
