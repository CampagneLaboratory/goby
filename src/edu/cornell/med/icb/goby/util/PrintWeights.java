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

package edu.cornell.med.icb.goby.util;

import edu.cornell.med.icb.goby.algorithmic.data.WeightsInfo;

import java.io.IOException;

/**
 * Print data in a weight file.
 * Usage: java -cp goby.jar edu.cornell.med.icb.goby.util.PrintWeights DLTTEJH-Bullard-HBR-SRR037439.gc-weights
 *
 * @author Fabien Campagne
 *         Date: Jun 23, 2010
 *         Time: 11:12:36 AM
 */
public class PrintWeights {
    private PrintWeights() {
    }

    public static void main(final String[] args) throws IOException, ClassNotFoundException {
        if (args.length == 0) {
            System.err.println("This utility requires exactly one arguments: the filename of the weights file to print.");
        }
        final String filename = args[0];

        final WeightsInfo weights = WeightsInfo.load(filename);
        final int size = weights.size();
        for (int i = 0; i < size; i++) {
            System.out.printf("readIndex %d weight %g%n", i, weights.getWeight(i));
        }
        System.exit(0);
    }
}
