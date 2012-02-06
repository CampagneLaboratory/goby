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
 * User: nyasha
 * Date: 1/27/12
 * Time: 4:08 PM
 */
public class MethylCountProviderTestSupport implements MethylCountProvider {


    private String[] samples;
    private int[] positions;
    private String chromosome;
    private int[][] C;
    private int[][] Cm;
    private String[] groups;

    public int getIndex() {
        return index;
    }

    private int index;


    @Override
    public CharSequence getChromosome() {
        return chromosome;
    }

    @Override
    public int getPosition() {
        return positions[index];
    }

    @Override
    public String[] getSamples() {
        return samples;
    }

    @Override
    public String[] getGroups(){
        return groups;
    }

    public MethylCountProviderTestSupport(String samples[],
                                          int positions[], String chromosome,
                                          int[][] C, int[][] Cm) {

        this.samples = samples;

        this.positions = positions;
        this.chromosome = chromosome;
        this.C = C;
        this.Cm = Cm;
    }


    public MethylCountProviderTestSupport(String[]groups, String samples[],
                                          int positions[], String chromosome,
                                          int[][] C, int[][] Cm) {
        this(samples,positions, chromosome, C, Cm);
        this.groups= groups;
    }
    @Override
    public int getC(int sampleIndex) {
        return C[sampleIndex][index];
    }

    @Override
    public int getCm(int sampleIndex) {
        return Cm[sampleIndex][index];
    }


    @Override
    public void next() {
        index++;
    }


}
