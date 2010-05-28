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
 * @author Fabien Campagne
 *         Date: May 21, 2010
 *         Time: 4:04:04 PM
 */
public class Log2RPKMCalculator extends StatisticCalculator {
    public boolean canDo(String[] group) {
        return true;
    }

    public DifferentialExpressionInfo evaluate(DifferentialExpressionCalculator differentialExpressionCalculator, NormalizationMethod method, DifferentialExpressionResults results, DifferentialExpressionInfo info, String... group) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
