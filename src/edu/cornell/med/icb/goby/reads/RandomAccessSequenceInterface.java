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
 * An interface to obtain bases from a RandomAccessSequence cache implementation.
 *
 * @author Fabien Campagne
 *         Date: May 1, 2011
 *         Time: 1:36:27 PM
 */
public interface RandomAccessSequenceInterface {
    /**
     * Return the base in the specified reference at the given position.
     *
     * @param referenceIndex index of the reference sequence.
     * @param position       zero-based position in the reference for which the base is sought.
     * @return base at position in reference sequence.
     */
    char get(final int referenceIndex, final int position);

    /**
     * Returns the length of the target sequence at index targetIndex.
     *
     * @param targetIndex Index of the target sequence.
     * @return Length of the sequence, in bases.
     */
    int getLength(int targetIndex);

    /**
     * Return bases from a range of positions in the cache.
     *
     * @param referenceIndex Name of the reference sequence
     * @param position       position where the range starts
     * @param length         length of the range for which bases should be returned.
     * @param bases          where the bases will be written.
     */
    void getRange(final int referenceIndex, final int position, final int length, MutableString bases);

    int getReferenceIndex(String referenceId);
   /**
     * Return the reference name corresponding to this index.
     *
     * @param index for the specified reference
     * @return referenceName The name of the sequence to get the index for
     */

    String getReferenceName(final int index) ;
    /**
     * Return the number of sequences in the cache.
     *
     * @return size of the cache.
     */
    int size();
}
