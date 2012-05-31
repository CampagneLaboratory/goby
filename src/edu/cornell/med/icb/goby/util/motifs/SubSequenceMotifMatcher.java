/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.util.motifs;

import it.unimi.dsi.lang.MutableString;

/**
 * Motif matcher that checks if a subsequence of accepted characters exactly match a pre-determined motif.
 *
 * @author campagne
 *         Date: 9/28/11
 *         Time: 12:28 PM
 */
public class SubSequenceMotifMatcher implements MotifMatcher {
    MutableString observed = new MutableString();
    MutableString motif;
    private int motifLength;

    public SubSequenceMotifMatcher(final String motif) {
        this.motif = new MutableString(motif).compact();
        motifLength=motif.length();
    }


    @Override
    public void accept(char c) {
        if (observed.length() == motifLength) {
            observed.deleteCharAt(0);
        }
        observed.append(c);
    }

    @Override
    public void newSequence() {
        observed.setLength(0);
    }

    @Override
    public boolean matched() {
        return motif.equals(observed);
    }
}
