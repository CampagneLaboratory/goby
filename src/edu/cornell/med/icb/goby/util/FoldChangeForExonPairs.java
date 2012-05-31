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

import edu.cornell.med.icb.io.TSVReader;
import edu.mssm.crover.cli.CLI;

import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import it.unimi.dsi.fastutil.objects.Object2DoubleMap;
import it.unimi.dsi.fastutil.objects.Object2DoubleArrayMap;
import it.unimi.dsi.fastutil.objects.Object2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.io.LineIterator;
import it.unimi.dsi.io.FastBufferedReader;
import org.apache.commons.io.FilenameUtils;

/**
 * @author Fabien Campagne
 *         Date: Aug 17, 2010
 *         Time: 2:00:35 PM
 */
public class FoldChangeForExonPairs {
    public static void main(String args[]) throws FileNotFoundException {
        FoldChangeForExonPairs processor = new FoldChangeForExonPairs();
        processor.process(args);
    }

    private void process(String[] args) throws FileNotFoundException {
        String inputFilenames = CLI.getOption(args, "--data", "expression-data.tsv,expression-data2.tsv");
        String pairsFilename = CLI.getOption(args, "--pairs", "pairs.tsv");
        String log2AverageId = CLI.getOption(args, "--group", "average log2_RPKM group UHR(BUQ)");
        String outputFilename = CLI.getOption(args, "--output", "/data/gc-weights/exons/out.tsv");
        String averageCountId = CLI.getOption(args, "--threshold-id", "average count group Brain");
        int thresholdValue = CLI.getIntOption(args, "--threshold-value", 20);

        //  String groupB= CLI.getOption(args,"-2","average RPKM group UHR(BUQ)");
        int log2AverageColumnIndex = -1;
        Object2DoubleMap<MutableString> exonExpressionData;


        ObjectArrayList<Pair> pairs = loadPairs(pairsFilename);
        PrintWriter output = new PrintWriter(outputFilename);
        output.printf("experiment\texonId1\texonId2\tlog2RpkmExon1-log2RpkmExon2%n");
        for (String inputFilename : inputFilenames.split("[,]")) {

            System.out.println("Processing " + inputFilename);
            exonExpressionData = loadData(inputFilename, log2AverageId, averageCountId, thresholdValue);
            System.out.println("Size: " + exonExpressionData.size());
            System.out.println("Writing data..");
            for (Pair pair : pairs) {
                //calculate log2 RPKM exonId1 - exonId2:
                if (exonExpressionData.containsKey(pair.exonId1) && exonExpressionData.containsKey(pair.exonId2)) {
                    double log2FoldChangeExons = exonExpressionData.get(pair.exonId1) - exonExpressionData.get(pair.exonId2);
                    output.printf("%s\t%s\t%s\t%g%n", FilenameUtils.getBaseName(inputFilename),
                            pair.exonId1, pair.exonId2, log2FoldChangeExons);
                }

            }
        }

    }

    private ObjectArrayList<Pair> loadPairs(String pairsFilename) throws FileNotFoundException {
        ObjectArrayList<Pair> pairs = new ObjectArrayList<Pair>();
        LineIterator lineIt = new LineIterator(new FastBufferedReader(new FileReader(pairsFilename)));
        while (lineIt.hasNext()) {
            MutableString mutableString = lineIt.next();
            String[] exonIds = mutableString.toString().split("[ \t]");
            Pair pair = new Pair();
            pair.exonId1 = new MutableString(exonIds[0]);
            pair.exonId2 = new MutableString(exonIds[1]);
            pairs.add(pair);
        }
        return pairs;
    }

    private Object2DoubleMap<MutableString> loadData(String inputFilename, String log2AverageId, String thresholdId,
                                                     int thresholdValue) {
        int log2AverageColumnIndex = -1;
        int thresholdColumnIndex = -1;
        Object2DoubleMap<MutableString> exonExpressionData = null;
        try {
            // load log2 RPKM expression values from the data table:
            TSVReader reader = new TSVReader(new FileReader(inputFilename));
            int numCols = -1;
            if (reader.hasNext()) {
                reader.next();
                numCols = reader.numTokens();
                for (int i = 0; i < numCols; i++) {
                    final String columnId = reader.getString();

                    if (log2AverageId.equals(columnId)) {
                        log2AverageColumnIndex = i;
                    }
                    if (thresholdId.equals(columnId)) {
                        thresholdColumnIndex = i;
                    }
                }

                if (log2AverageColumnIndex == -1) {
                    System.err.printf("Could not find the group log2 average column. Aborting.", log2AverageId);
                    System.exit(1);
                }
                if (thresholdColumnIndex == -1) {
                    System.err.printf("Could not find the threshold column. Aborting.", thresholdId);
                    System.exit(1);
                }


            }
            System.out.println("Loading data..");
            exonExpressionData = new Object2DoubleOpenHashMap<MutableString>();

            while (reader.hasNext()) {
                reader.next();
                String elementId = reader.getString();
                double log2Rpkm = 0;
                double averageCount = 0;

                for (int i = 1; i < numCols; i++) {

                    if (i == log2AverageColumnIndex) {
                        log2Rpkm = reader.getDouble();

                    } else if (i == thresholdColumnIndex) {
                        averageCount = reader.getDouble();

                    } else {
                        reader.getString();
                    }

                }
                if (averageCount >= thresholdValue) {
                    exonExpressionData.put(new MutableString(elementId), log2Rpkm);
                }

            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return exonExpressionData;
    }

    private class Pair {
        MutableString exonId1;
        MutableString exonId2;
    }
}
