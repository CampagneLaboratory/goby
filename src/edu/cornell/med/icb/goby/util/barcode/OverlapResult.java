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

package edu.cornell.med.icb.goby.util.barcode;

/**
 * Describe class here.
 *
 * @author Kevin Dorff
 */
public class OverlapResult {

    public int start;

    public int length;

    public OverlapResult() {
    }

    public OverlapResult(final int start, final int length) {
        this.start = start;
        this.length = length;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        final OverlapResult that = (OverlapResult) o;

        if (length != that.length) {
            return false;
        }
        if (start != that.start) {
            return false;
        }

        return true;
    }

    @Override
    public int hashCode() {
        int result = start;
        result = 31 * result + length;
        return result;
    }
}
