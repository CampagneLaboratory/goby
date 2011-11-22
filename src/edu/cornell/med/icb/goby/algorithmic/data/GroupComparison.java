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

package edu.cornell.med.icb.goby.algorithmic.data;

/**
 * Describe two groups under comparison.
 * @author Fabien Campagne
 *         Date: 11/21/11
 *         Time: 3:33 PM
 */
public class GroupComparison {
    final public String nameGroup1;
    final public String nameGroup2;
    final public int indexGroup1;
    final public int indexGroup2;
    /**
     * Index of this comparison among all during study.
     */
    final public int index;
    public GroupComparison(String nameGroup1, String nameGroup2, int indexGroup1, int indexGroup2, int comparisonIndex) {
        this.nameGroup1 = nameGroup1;
        this.nameGroup2 = nameGroup2;
        this.indexGroup1 = indexGroup1;
        this.indexGroup2 = indexGroup2;
        this.index=comparisonIndex;
    }
}
