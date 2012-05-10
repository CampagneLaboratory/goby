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

package edu.cornell.med.icb.goby.alignments;

import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;

import java.io.Closeable;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Properties;

/**
 * An interface for alignment writers.
 * @author Fabien Campagne
 *         Date: 4/17/12
 *         Time: 5:28 PM
 */
public interface AlignmentWriter extends Closeable {
    void setPermutation(boolean state);

    void setSorted(boolean sortedState);

    void appendEntry(Alignments.AlignmentEntry entry) throws IOException;

    @Override
    void close() throws IOException;

    void setQueryIdentifiersArray(String[] queryIdentifiersArray);

    void setTargetIdentifiersArray(String[] targetIdentifiersArray);

    void setQueryIdentifiers(IndexedIdentifier queryIdentifiers);

    void setTargetIdentifiers(IndexedIdentifier targetIdentifiers);

    void setTargetLengths(int[] targetLengths);

    void setNumQueries(int numQueries);

    void setNumTargets(int numTargets);

    void putStatistic(String description, String value);

    void putStatistic(String description, double value);

    void putStatistic(String description, int value);

    void setAlignerVersion(String alignerVersion);

    void setAlignerName(String alignerName);

    void setReadOriginInfo(ObjectArrayList<Alignments.ReadOriginInfo.Builder> readOriginInfoBuilderList);
    void addReadOriginInfo(ObjectArrayList<Alignments.ReadOriginInfo.Builder> readOriginInfoBuilderList);

    void printStats(PrintStream out);

    void setStatistics(Properties statistics);

    void setSmallestSplitQueryIndex(int smallestSplitQueryIndex);

    void setLargestSplitQueryIndex(int largestSplitQueryIndex);
}
