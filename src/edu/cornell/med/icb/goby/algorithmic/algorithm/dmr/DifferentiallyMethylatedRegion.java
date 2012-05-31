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

/**
 * Defines a differentially methylated region of the genome
 * @author  Nyasha Chambwe
 *          Date: 4/18/12
 *          Time: 11:02 AM
 */

public class DifferentiallyMethylatedRegion {
    private final int chromosome;//target index
    private final int start;
    private final int end;
    /* defines the number of cytosines in the DMR window*/
    private final int numCytosines;

    public DifferentiallyMethylatedRegion(int chromosome, int end, int start, int numCytosines) {
        this.chromosome = chromosome;
        this.end = end;
        this.start = start;
        this.numCytosines= numCytosines;
    }


}
