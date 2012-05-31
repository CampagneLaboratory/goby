/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 *  This file is part of the Goby IO API.
 *
 *     The Goby IO API is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     The Goby IO API is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with the Goby IO API.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.cornell.med.icb.goby.counts;

/**
 * Information about a peak of something (e.g., read counts) along a reference.
 *
 * @author Fabien Campagne
 *         Date: May 27, 2009
 *         Time: 6:44:06 PM
 */
public class Peak {
    public int start;
    public int count;
    public int length;

    /**
     * Create a copy of this peak information.
     *
     * @return a new instance of Peak with the same information.
     */
    public Peak copy() {
        final Peak result = new Peak();
        result.start = this.start;
        result.count = this.count;
        result.length = this.length;
        return result;
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder();
        sb.append(" peak :");
        sb.append(" start :");
        sb.append(start);
        sb.append(" count :");
        sb.append(count);
        sb.append(" length :");
        sb.append(length);
        return sb.toString();
    }
}
