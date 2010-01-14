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

import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.fastutil.doubles.DoubleList;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;

import java.io.PrintWriter;

/**
 * @author Fabien Campagne
 *         Date: Jan 11, 2010
 *         Time: 6:57:49 PM
 */
public class DifferentialExpressionInfo {
    MutableString elementId;
    DoubleArrayList statistics = new DoubleArrayList();

    @Override
    public String toString() {
        return String.format("[ %s %s ]", elementId, statistics.toString());
    }

    public void write(PrintWriter printWriter, char delimiter) {
        printWriter.append(elementId);
        for (double value : statistics) {
            printWriter.append(delimiter);
            printWriter.append(String.format("%g", value));
        }
    }
}
