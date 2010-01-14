/*
 * Copyright (C) 2010 Institute for Computational Biomedicine,
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
 * Abstract class for all implementations of FDR correction methods.
 *  @author Fabien Campagne
 *         Date: Jan 12, 2010
 *         Time: 6:59:01 PM
 */
public abstract class FDRAdjustment {

    public DifferentialExpressionResults adjust(DifferentialExpressionResults list, String... statisticIds) {
        for (String statisticId : statisticIds) {
            adjust(list, statisticId);
        }
        return list;
    }

    abstract public DifferentialExpressionResults adjust(DifferentialExpressionResults list, String statisticId);
}
