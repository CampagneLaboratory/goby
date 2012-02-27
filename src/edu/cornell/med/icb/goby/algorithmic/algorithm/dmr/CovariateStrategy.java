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

package edu.cornell.med.icb.goby.algorithmic.algorithm.dmr;

import java.io.Serializable;

/**
 * Strategies to determine how to map covariates to discrete distributions.
 * @author Fabien Campagne
 *         Date: 2/26/12
 *         Time: 12:43 PM
 */
public abstract class CovariateStrategy implements Serializable {
    /**
     * Take a set of covariates and return the index of the distribution that the covariates map to.
     * @param covariates list of covariates associated with an observation.
     * @return an index between zero and the number of discrete distributions.
     */
    public abstract int getIndex(int... covariates);
}
