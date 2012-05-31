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

package edu.cornell.med.icb.goby.alignments;

/**
 * Stores statistics about sequence variations in a specific sample.
 *
 * @author Fabien Campagne
*         Date: Mar 21, 2011
*         Time: 12:27:53 PM
*/
public class ReadIndexStats {
    /**
     * The basename of the sample from which stats were estimated.
     */
    public String basename;
    /**
     * The index of the alignment reader that is reading over this basename, will be populated when we know.
     */
    public int readerIndex = -1;
    /**
     * Indexed by readIndex
     */
    public long[] countVariationBases;
    /**
     * Indexed by readIndex
     */
    public long[] countReferenceBases;
}
