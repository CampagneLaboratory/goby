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

package edu.cornell.med.icb.goby.counts;

/**
 * @author Fabien Campagne
 *         Date: 5/27/11
 *         Time: 2:43 PM
 */
public interface CountsAggregatorI extends CountsReaderI {
     /**
     * Return the sum of counts over the readers that have non zero counts at the current position.
     */
    public int getCount() ;

    /**
     * Return the count for the reader at index readerIndex.
     * @param readerIndex Index of the reader for which count is sought.
     * @return count or the current transition for the specified reader.
     */
     public int getCount(final int readerIndex);
}
