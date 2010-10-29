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

package edu.cornell.med.icb.goby.alignments;

import java.util.Comparator;

/**
 * Compare alignment entries by reference/target sequence then by position.
 *
 * @author Fabien Campagne
 *         Date: Jun 21, 2010
 *         Time: 6:03:58 PM
 */
public class AlignmentPositionComparator implements Comparator<Alignments.AlignmentEntry> {


    public final int compare(final Alignments.AlignmentEntry a,
                             final Alignments.AlignmentEntry b) {

        final int sameTargets = a.getTargetIndex() - b.getTargetIndex();
        if (sameTargets == 0) {
            return a.getPosition() - b.getPosition();
        } else {
            return sameTargets;
        }

    }
}

