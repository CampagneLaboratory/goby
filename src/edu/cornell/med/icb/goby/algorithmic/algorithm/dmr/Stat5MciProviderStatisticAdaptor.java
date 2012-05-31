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

import edu.cornell.med.icb.goby.algorithmic.data.MethylCountInfo;

/**
 * A stat5 implementation that takes Cm and C counts from a MethylCountInfo instance. Used to estimate empirical
 * P-values for individual sites.
 *
 * @author Fabien Campagne
 *         Date: 2/29/12
 *         Time: 6:36 PM
 */
public class Stat5MciProviderStatisticAdaptor extends Stat5StatisticAdaptor {
    private static final long serialVersionUID = -6171007111030942764L;

    public String statName() {
        return "stat5_mci";
    }

    @Override
    protected void checkProviderType(final Object dataProvider) {
        if (!(dataProvider instanceof MethylCountInfo)) {
            throw new InternalError();
        }

    }

    @Override
    protected int getC(final int sampleIndexA, final Object sampleDataPool, final int contextIndex) {
        final MethylCountInfo mci = (MethylCountInfo) sampleDataPool;
        return mci.unmethylatedCCountPerSample[sampleIndexA];
    }

    @Override
    protected int getCm(final int sampleIndexA, final Object sampleDataPool, final int contextIndex) {
        final MethylCountInfo mci = (MethylCountInfo) sampleDataPool;
        return mci.methylatedCCountPerSample[sampleIndexA];
    }
}
