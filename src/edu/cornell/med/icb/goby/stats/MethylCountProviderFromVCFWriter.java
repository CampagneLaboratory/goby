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
 * Time: 4:21 PM
 */
public class MethylCountProviderFromVCFWriter implements MethylCountProvider {
    private VCFWriter vcfWriter;
    private int[][] C;

    private int[][] Cm;
    private String[] samples;

    public MethylCountProviderFromVCFWriter(VCFWriter vcfWriter) {
        this.vcfWriter = vcfWriter;
    }

    @Override
    public CharSequence getChromosome() {
        return vcfWriter.getChromosome();
    }

    @Override
    public int getPosition() {
        return vcfWriter.getPosition();
    }

    @Override
    public String[] getSamples() {
        return vcfWriter.getSamples();
    }

    @Override
    public int getC(int sampleIndex, int positionIndex) {
        return C[sampleIndex][positionIndex];
    }

    @Override
    public int getCm(int sampleIndex, int positionIndex) {
              return Cm[sampleIndex][positionIndex];
    }

    @Override
    public void next() {
        vcfWriter.clear();
    }


    @Override
    public void feedFrom(VCFWriter vcfAveragingWriter) {
    // TODO complete this
        throw new InternalError("Complete this section.");
        /*
        this.samples = vcfWriter.getSamples();
        for (int sampleIndex = 0; sampleIndex < samples.length; sampleIndex++) {
            C[sampleIndex] = Integer.parseInt(vcfAveragingWriter.getSampleValue(fieldIndexC, sampleIndex));
            Cm[sampleIndex] = Integer.parseInt(vcfAveragingWriter.getSampleValue(fieldIndexCm, sampleIndex));
        }*/
    }

    @Override
    public int getIndex(){
        return 0;
    }
}
