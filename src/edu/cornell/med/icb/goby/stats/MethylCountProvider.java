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

package edu.cornell.med.icb.goby.stats;

/**
 * @author Fabien Campagne
 * @author Nyasha Chambwe
 * Date: 1/27/12
 * Time: 4:06 PM
 */
public interface MethylCountProvider {
    /**
     * Get the chromosome for which the counts are returned.
     * @return reference sequence identifier.
     */
    public CharSequence getChromosome();

    /**
     * Get the position for which the counts are returned.
     * @return position one-based position in chromosome.
     */
    public int getPosition();

    /**
     * Get the sample ids for each index.
     * @return array of sample identifiers.
     */
    public String[] getSamples();

    /**
     * Get the count of unmethylated cytosines in a given sample.
     *
     * @param sampleIndex index of the sample.
     * @return number of unmethylated cytosine bases.
     */
    public int getC(int sampleIndex);

    /**
     * Get the count of methylated cytosines in a given sample.
     *
     * @param sampleIndex index of the sample.
     * @return number of methylated cytosine bases.
     */
    public int getCm(int sampleIndex);

    /**
     * Advance to the next position.
     */
    public void next();


}
