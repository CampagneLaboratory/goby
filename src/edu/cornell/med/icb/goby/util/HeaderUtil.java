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

package edu.cornell.med.icb.goby.util;

import edu.cornell.med.icb.goby.alignments.AlignmentWriter;
import edu.cornell.med.icb.goby.alignments.ConcatAlignmentReader;

/**
 * @author Fabien Campagne
 *         Date: 5/9/12
 *         Time: 12:01 PM
 */
public class HeaderUtil {

    /**
     * Copy header information from the source to the destination.
     *
     * @param source      AlignmentReader where to obtain reader information from.
     * @param destination AlignmentWriter where to write header information to.
     */
    public static void copyHeader(final ConcatAlignmentReader source, final AlignmentWriter destination) {
        destination.setSmallestSplitQueryIndex(source.getSmallestSplitQueryIndex());
        destination.setLargestSplitQueryIndex(source.getLargestSplitQueryIndex());

        if (source.getQueryIdentifiers() != null) {
            destination.setQueryIdentifiers(source.getQueryIdentifiers());
        }
        if (source.getTargetIdentifiers() != null) {
            destination.setTargetIdentifiers(source.getTargetIdentifiers());
        }
        if (source.getTargetLength() != null) {
           destination.setTargetLengths(source.getTargetLength());
        }
        destination.setNumQueries(source.getNumberOfQueries());
        destination.setNumTargets(source.getNumberOfTargets());

        destination.setStatistics(source.getStatistics());

    }
}
