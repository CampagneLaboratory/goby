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

package edu.cornell.med.icb.goby.algorithmic.data;

/**
* @author Fabien Campagne
*         Date: 7/27/12
*         Time: 12:33 PM
*/
public class Interval {
    public int referenceIndex;
    public int start;
    public int end;
    public String id;

    @Override
    public String toString() {
        return "Interval{" +
                "end=" + end +
                ", referenceIndex=" + referenceIndex +
                ", start=" + start +
                ", id='" + id + '\'' +
                '}';
    }
}
