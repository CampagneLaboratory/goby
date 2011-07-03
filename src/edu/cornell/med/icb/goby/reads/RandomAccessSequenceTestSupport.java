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

package edu.cornell.med.icb.goby.reads;

import it.unimi.dsi.lang.MutableString;

/**
 * @author Fabien Campagne
 *         Date: May 1, 2011
 *         Time: 1:46:38 PM
 */
public class RandomAccessSequenceTestSupport implements RandomAccessSequenceInterface {
    private String[] referenceSequences;

    public RandomAccessSequenceTestSupport(String[] referenceSequences) {
        this.referenceSequences = referenceSequences;
    }


    public char get(int referenceIndex, int position) {
        return referenceSequences[referenceIndex].charAt(position);
    }

    public int getLength(int targetIndex) {
        return referenceSequences[targetIndex].length();
    }

    @Override
    public void getRange(int referenceIndex, int position, int length, MutableString bases) {
        bases.setLength(0);
        bases.append(referenceSequences[referenceIndex].substring(position, position + length));
    }

    /**
     * Override this method to return another integer than zero.
     *
     * @param referenceId an id
     * @return zero
     */
    @Override
    public int getReferenceIndex(String referenceId) {
        return 0;
    }

    @Override
    public int size() {
        return referenceSequences.length;
    }
}
