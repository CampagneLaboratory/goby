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

package edu.cornell.med.icb.goby.algorithmic.algorithm;

import edu.cornell.med.icb.goby.algorithmic.data.HeptamerInfo;
import it.unimi.dsi.lang.MutableString;

/**
 * @author Fabien Campagne
 *         Date: May 21, 2010
 *         Time: 5:07:17 PM
 */
public class HeptamerWeight implements WeightCalculator {
    private HeptamerInfo heptamers;

    public HeptamerWeight(final HeptamerInfo heptamers) {
        this.heptamers = heptamers;
    }

    public float weight(final MutableString sequence) {
        if (heptamers.colorSpace) {
            sequence.delete(0, 1);
        }
        // if (count++ > 1000000) break;
        final int item = 0;

        final int positionInRead = 1;


        final int end = positionInRead - 1 + heptamers.heptamerLength;
        final int start = positionInRead - 1;
        //     System.out.printf("%d %d %d %d%n", positionInRead,sequence.length(),start, end );
        final MutableString heptamer = sequence.substring(start, end);

        if (heptamer.indexOf('N') == -1) {
            // heptamers that include any number of Ns are ignored.
            final short heptamerIndex = (short) heptamers.heptamerToIndices.getInt(heptamer);

            final float weight = heptamerIndex == -1 ? 1 : heptamers.heptamerIndexToWeight.get(heptamerIndex);
            return weight;
        } else {
            return 1f;
        }

    }
     public String id() {
        return "heptamers";
    }
}
