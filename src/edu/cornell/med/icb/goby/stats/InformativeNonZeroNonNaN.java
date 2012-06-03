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

package edu.cornell.med.icb.goby.stats;

/**
 * Class defining informative as being non-zero and non-NaN.
 *
 * @author Kevin Dorff
 */
public class InformativeNonZeroNonNaN implements InformativeDouble {
    /**
     * Check if a value is not NaN and not zero.
     * @param value the value to check
     * @return true if the value is informative (non-zero, non-NaN).
     */
    public boolean isInformative(final double value) {
        boolean informative = false;
        if (value != 0 && !Double.isNaN(value)) {
            // require something else than NaN or zero to be have an informative DE.
            informative = true;
        }
        return informative;
    }
}
