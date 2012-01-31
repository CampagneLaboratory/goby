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

import edu.cornell.med.icb.goby.algorithmic.data.MethylCountInfo;
import edu.cornell.med.icb.goby.modes.MethylationRegionsOutputFormat;

/**
 * A provider that fetches methylation counts from a region output format.
 * @author Fabien Campagne
 *         Date: 1/31/12
 *         Time: 1:09 PM
 */
public class MethylCountProviderFromRegionsOutputFormat implements MethylCountProvider {
    private MethylCountInfo mci;

    public MethylCountProviderFromRegionsOutputFormat(MethylationRegionsOutputFormat regionFormat) {
        this.regionFormat = regionFormat;
        mci=regionFormat.getMci();
    }

    MethylationRegionsOutputFormat regionFormat;

    @Override
    public CharSequence getChromosome() {
        return regionFormat.getChromosome();
    }

    @Override
    public int getPosition() {
        return regionFormat.getPosition();
    }

    @Override
    public String[] getSamples() {
        return regionFormat.getSamples();
    }

    @Override
    public int getC(int sampleIndex) {
        return mci.unmethylatedCCountPerSample[sampleIndex];
    }

    @Override
    public int getCm(int sampleIndex) {
         return mci.methylatedCCountsPerSample[sampleIndex];
    }

    @Override
    public void next() {
       mci.reset();
    }



}
