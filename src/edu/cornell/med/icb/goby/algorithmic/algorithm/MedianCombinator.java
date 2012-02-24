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

import it.unimi.dsi.fastutil.doubles.DoubleArrayList;

import java.util.Collections;

/**
 * A combinator that returns the median of a set of p-values.
 *
 * @author Fabien Campagne
 *         Date: 2/23/12
 *         Time: 2:49 PM
 */
public class MedianCombinator implements EvidenceCombinator {

    private DoubleArrayList list = new DoubleArrayList();

    @Override
    public void observe(double pValue) {
        list.add(pValue);

    }

    @Override
    public void reset() {
        list.clear();

    }

    @Override
    public double adjust() {
        Collections.sort(list);
        int size = list.size();
        if (size > 0) {
            return list.get(size / 2);
        } else {
            return 1.0;
        }
    }
}
