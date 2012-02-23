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

package edu.cornell.med.icb.goby.algorithmic.algorithm;

/**
 * An interface for algorithms that combine evidence across tests.
 * @author Fabien Campagne
 *         Date: 2/23/12
 *         Time: 2:37 PM
 */
public interface EvidenceCombinator {
    /**
     * Observe a pValue from one test. This method can be called repetitively for multiple tests.
     * @param pValue the pValue obtained on one test.
     */
    void observe(double pValue);

    /**
     * Reset the combinator to construction state.
     */
    void reset();

    /**
     * Derive the adjusted p-value across the test(s) observed. By convention, when no p-values were observed, the adjust method
     * will return 1.0.
     * @return q-value or adjusted p-value.
     */
    double adjust();
}
