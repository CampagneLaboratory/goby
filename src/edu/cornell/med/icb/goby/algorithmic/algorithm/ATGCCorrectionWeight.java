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

package edu.cornell.med.icb.goby.algorithmic.algorithm;

import it.unimi.dsi.lang.MutableString;

/**
 * @author Fabien Campagne
 *         Date: May 21, 2010
 *         Time: 5:19:54 PM
 */
public class ATGCCorrectionWeight implements WeightCalculator {

    public ATGCCorrectionWeight(final boolean colorSpace) {
        if (colorSpace) {
            throw new UnsupportedOperationException("ATGC content is not implemented for color-space reads");
        }
    }


    public float weight(final MutableString sequence) {

        float GC = 0;
        float AT = 0;

        for (int i = 0; i < sequence.length(); i++) {

            final char c = sequence.charAt(i);
            GC += (c == 'G' || c == 'C') ? 1 : 0;
            AT += (c == 'A' || c == 'T') ? 1 : 0;
        }

        final float normGC = GC / (AT + GC);
        return (10.4436f - 8.266972f * normGC);
    }

    public String id() {
        return "atgc";
    }
}
