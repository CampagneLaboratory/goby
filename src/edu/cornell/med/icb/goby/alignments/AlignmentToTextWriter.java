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
import it.unimi.dsi.lang.MutableString;

import java.io.IOException;
import java.io.PrintStream;
import java.util.Properties;

/**
 * @author Fabien Campagne
 *         Date: 4/17/12
 *         Time: 5:50 PM
 */
public class AlignmentToTextWriter implements AlignmentWriter {
    private MutableString textOutput = new MutableString();

    public MutableString getTextOutput() {
        return textOutput;
    }

    public AlignmentToTextWriter() {

    }

    @Override
    public void setPermutation(boolean state) {
        textOutput.append("Set permutation=");
        textOutput.append(state);
        textOutput.append("\n");
    }

    @Override
    public void setSorted(boolean sortedState) {
        textOutput.append("Set sortedState=");
        textOutput.append(sortedState);
        textOutput.append("\n");
    }

    @Override
    public void appendEntry(Alignments.AlignmentEntry entry) throws IOException {
        textOutput.append("{");
        textOutput.append(entry.toString());
        textOutput.append("}");
        textOutput.append("\n");
    }

    @Override
    public void close() throws IOException {
        textOutput.append("Closed");
        textOutput.append("\n");
    }

    @Override
    public void setQueryIdentifiersArray(String[] queryIdentifiersArray) {
        textOutput.append("Set queryIndentifiersArray, length=");
        textOutput.append(queryIdentifiersArray.length);
    }

    @Override
    public void setTargetIdentifiersArray(String[] targetIdentifiersArray) {
        textOutput.append("Set setTargetIdentifiersArray:");
        textOutput.append(ObjectArrayList.wrap(targetIdentifiersArray).toString());
        textOutput.append("\n");
    }

    @Override
    public void setQueryIdentifiers(IndexedIdentifier queryIdentifiers) {
        throw new UnsupportedOperationException("This method has not yet been implemented");
    }

    @Override
    public void setTargetIdentifiers(IndexedIdentifier targetIdentifiers) {
        textOutput.append("Set targetIdentifiers: {");
        for (MutableString key : targetIdentifiers.keySet()) {
            textOutput.append(targetIdentifiers.getInt(key));
            textOutput.append("->");
            textOutput.append(key);
            textOutput.append("\n");
        }
        textOutput.append("}\n");
    }

    @Override
    public void setTargetLengths(int[] targetLengths) {
        textOutput.append("Set targetLengths: {");
        int i = 0;
        for (int length : targetLengths) {

            textOutput.append(i++);
            textOutput.append("->");
            textOutput.append(length);
            textOutput.append("\n");
        }
        textOutput.append("}\n");
    }

    @Override
    public void setNumQueries(int numQueries) {
        textOutput.append("Set numQueries=");
        textOutput.append(numQueries);
        textOutput.append("\n");
    }

    @Override
    public void setNumTargets(int numTargets) {
        textOutput.append("Set numTargets=");
        textOutput.append(numTargets);
        textOutput.append("\n");
    }

    @Override
    public void putStatistic(String description, String value) {
        textOutput.append(String.format("Put statistic description=%s value=%s %n", description, value));

    }

    @Override
    public void putStatistic(String description, double value) {
        textOutput.append(String.format("Put statistic description=%s value=%g %n", description, value));
    }

    @Override
    public void putStatistic(String description, int value) {
        textOutput.append(String.format("Put statistic description=%s value=%d %n", description, value));
    }

    @Override
    public void setAlignerVersion(String alignerVersion) {
        textOutput.append(String.format("Set aligner version=%s %n", alignerVersion));
    }

    @Override
    public void setAlignerName(String alignerName) {
        textOutput.append(String.format("Set aligner version=%s %n", alignerName));
    }

    @Override
    public void setReadOriginInfo(ObjectArrayList<Alignments.ReadOriginInfo.Builder> readOriginInfoBuilderList) {
        throw new UnsupportedOperationException("This method has not yet been implemented");
    }

    @Override
    public void addReadOriginInfo(ObjectArrayList<Alignments.ReadOriginInfo.Builder> readOriginInfoBuilderList) {
        throw new UnsupportedOperationException("This method has not yet been implemented");
    }

    @Override
    public void printStats(PrintStream out) {
        throw new UnsupportedOperationException("This method has not yet been implemented");
    }

    @Override
    public void setStatistics(Properties statistics) {
        textOutput.append(String.format("Set statistics { %s }%n", statistics.toString()));
    }

    @Override
    public void setSmallestSplitQueryIndex(int smallestSplitQueryIndex) {
        textOutput.append(String.format("Set smallestSplitQueryIndex=%d%n", smallestSplitQueryIndex));
    }

    @Override
    public void setLargestSplitQueryIndex(int largestSplitQueryIndex) {
        textOutput.append(String.format("Set largestSplitQueryIndex=%d%n", largestSplitQueryIndex));
    }
}
