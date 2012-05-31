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

import edu.cornell.med.icb.goby.algorithmic.algorithm.FenwickTree;
import edu.cornell.med.icb.goby.algorithmic.algorithm.dmr.BinningStrategy;
import edu.cornell.med.icb.goby.algorithmic.algorithm.dmr.EstimatedDistribution;
import edu.cornell.med.icb.goby.algorithmic.algorithm.dmr.StatisticAdaptor;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.ints.Int2ObjectArrayMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.io.BinIO;
import it.unimi.dsi.fastutil.objects.Object2IntArrayMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.lang.MutableString;

import java.io.IOException;
import java.io.Serializable;

/**
 * Provide access to estimated test distributions.
 * @author Fabien Campagne
 *         Date: 5/1/12
 *         Time: 2:56 PM
 */
public class EstimatedTestDistributions implements Serializable {

    private static final long serialVersionUID = -4803501043413548993L;
    private ObjectArrayList<EstimatedDistribution> testDistributions=new ObjectArrayList<EstimatedDistribution>();
    private Object2IntOpenHashMap<String> groupToIndex = new Object2IntOpenHashMap<String>();
    private int currentIndex;
    private int numberOfContexts;
    private BinningStrategy binningStrategy;
    private StatisticAdaptor statAdaptor;


    public EstimatedTestDistributions() {
        groupToIndex.defaultReturnValue(-1);
    }

    public EstimatedTestDistributions(int numberOfContexts, StatisticAdaptor statAdaptor) {
        this.numberOfContexts=numberOfContexts;
        groupToIndex.defaultReturnValue(-1);
        this.statAdaptor=statAdaptor;
    }

    public EstimatedDistribution getTestDistribution(String groupComparison) {
        int index = groupToIndex.getInt(groupComparison);
        if (index == -1) {
            final EstimatedDistribution distribution = new EstimatedDistribution(numberOfContexts,statAdaptor);
            distribution.setBinningStrategy(binningStrategy);

            testDistributions.add(distribution);
            groupToIndex.put(groupComparison, currentIndex);
            index = currentIndex;
            currentIndex++;
        }
        return testDistributions.get(index);
    }

    public FenwickTree getDensity(String groupComparison, int... covariates) {
        return getTestDistribution(groupComparison).getDensity(covariates);
    }

    public void setBinningStrategy(BinningStrategy binningStrategy) {
        this.binningStrategy=binningStrategy;
    }

    public static EstimatedTestDistributions load(String filename) throws ClassNotFoundException, IOException {
        return (EstimatedTestDistributions) BinIO.loadObject(filename);
    }

    public StatisticAdaptor getStatAdaptor() {
        return statAdaptor;
    }

    public static void store(EstimatedTestDistributions testDistributions, String testFilename) throws IOException {
        BinIO.storeObject(testDistributions,testFilename);
    }
}
