/*
 * Copyright (C) 2010 Institute for Computational Biomedicine,
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
import it.unimi.dsi.fastutil.floats.FloatArrayList;
import it.unimi.dsi.fastutil.floats.FloatList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.logging.ProgressLogger;

import java.io.FileReader;
import java.io.IOException;

/**
 *
 * @author Fabien Campagne
 *         Date: Oct 21, 2010
 *         Time: 3:54:56 PM
 */
public class MethylSimilarityScan {
    private int scoreComparator;

    public static void main(String args[]) throws IOException {
        MethylSimilarityScan scanner = new MethylSimilarityScan();
        scanner.process(args);

    }

    private class Query {
        int positionStart;
        int positionEnd;
        FloatList query;

    }

    private void process(String[] args) throws IOException {
        String inputFilename = CLI.getOption(args, "-i", "/data/lister/mc_h1.tsv");
        String chromosomeSelected = CLI.getOption(args, "-c", "1");
        char strandSelected = CLI.getOption(args, "-s", "-").charAt(0);
        int queryStartPosition = CLI.getIntOption(args, "-qs", -1);
        int queryEndPosition = CLI.getIntOption(args, "-qe", -1);

        FloatArrayList methylationRates = load(inputFilename, chromosomeSelected, strandSelected);
        Query query;

        if (queryStartPosition == -1 || queryEndPosition == -1) {

            int minSites = 5;
            int windowWidth = 20;
            System.out.printf("minSites=%d windowWidth=%d %n", minSites, windowWidth);
            query = findQuery(methylationRates, windowWidth, minSites);
            //  System.out.println("You must specify query start (-qs) and query end positions (-qe).");
            //  System.exit(1);
        } else {
            query = new Query();
            query.positionEnd = queryEndPosition;
            query.positionStart = queryStartPosition;
            query.query = methylationRates.subList(queryStartPosition, queryEndPosition);

        }


        for (int position = 0; position < query.query.size(); position++) {
            System.out.printf("position %d rate: %f %n", query.positionStart + position, query.query.get(position));

        }
        HitBoundedPriorityQueue results = bruteForceSearch(query, methylationRates, 100);

        printResults(chromosomeSelected, strandSelected, results, query, methylationRates);
    }

    private Query findQuery(FloatArrayList methylationRates, int windowWidth, int minSites) {
        int windowTally = 0;
        IntArrayList startTally = new IntArrayList();
        IntArrayList endTally = new IntArrayList();
        int startAdjustment = 0;
        int endAdjustment = 0;

        for (int position = 0; position < methylationRates.size(); position++) {
            final float rate = methylationRates.getFloat(position);
            if (rate == rate) {
                startAdjustment++;
            }
            final float previousRate = position > 1 ? methylationRates.get(position - 1) : Float.NaN;
            if (previousRate == previousRate) {
                endAdjustment++;
            }
            startTally.add(startAdjustment);
            endTally.add(endAdjustment);

        }
        for (int start = 0; start < methylationRates.size(); start++) {
            final int end = start + windowWidth;
            if (endTally.get(end) - startTally.get(start) >= minSites) {
                System.out.printf("queryStart=%d queryEnd=%d%n", start, end);
                // return the first window of size windowSize with at least minSites in the window
                Query query = new Query();
                query.positionEnd = end;
                query.positionStart = start;
                query.query = methylationRates.subList(start, end);
                return query;
            }
        }
        return null;
    }


    private void printResults(String chromosome, char strand, HitBoundedPriorityQueue results,
                              MethylSimilarityScan.Query query, FloatList target) {
        ObjectArrayList<MethylationSimilarityMatch> sortedHits = new ObjectArrayList();
        sortedHits.size(results.size());

        int queryLength = query.query.size();
        int i = sortedHits.size();
        while (!results.isEmpty()) {
            MethylationSimilarityMatch hit = results.dequeue();
            sortedHits.set(--i, hit);
        }
        for (MethylationSimilarityMatch hit : sortedHits) {
            System.out.printf("query %s %c %d-%d [ %d-%d %f]%n", chromosome, strand,
                    query.positionStart, query.positionEnd,
                    hit.targetPosition, hit.targetPosition + queryLength, hit.score);
        }
        for (MethylationSimilarityMatch hit : sortedHits) {
            System.out.printf("chr\tstrand\tstart-end\tquery-position\tquery\ttarget%n");
            for (int queryPosition = 0; queryPosition < queryLength; queryPosition++) {
                final float queryRate = query.query.getFloat(queryPosition);
                final float targetRate = hit.targetPosition + queryPosition < target.size() ?
                        target.getFloat(hit.targetPosition + queryPosition) :
                        Float.NaN;
                if (queryRate == queryRate && targetRate == targetRate) {
                    // print only matching positions:
                    System.out.printf("%s\t%c\t%d-%d\t%d\t%f\t%f%n", chromosome, strand,
                            query.positionStart, query.positionEnd, queryPosition,
                            queryRate,
                            targetRate);
                }

            }
        }

    }

    private HitBoundedPriorityQueue bruteForceSearch(Query query, FloatArrayList target, int maxBestHits) {
        int queryLength = query.query.size();
        int targetLength = target.size();
        HitBoundedPriorityQueue results = new HitBoundedPriorityQueue(maxBestHits);

        ProgressLogger pg = new ProgressLogger();
        //  int delta=(queryLength / 10);
        int delta = 1;
        pg.expectedUpdates = targetLength / delta;

        pg.start("brute force");
        IntArrayList queryActivePositions;
        IntArrayList targetActivePositions;
        queryActivePositions = getActivePositions(query);
        //  targetActivePositions = getActivePositions(target);

        final int earlyStop = Integer.MAX_VALUE;
        for (int targetPosition = 0; targetPosition < Math.min(earlyStop, targetLength); targetPosition += delta) {
            float score = 0;
            int lastQueryPosition = 0;

            for (int queryPosition : queryActivePositions) {


                final float queryRate = query.query.getFloat(queryPosition);
                final int targetIndex = targetPosition + queryPosition;
                if (targetIndex >= targetLength) break;
                final float targetRate = target.getFloat(targetIndex);

                if (targetRate == targetRate) {
                    final float scoreContribution = 1 - (targetRate - queryRate) * (targetRate - queryRate);
                    //    System.out.printf("matching targetIndex %d indexInQuery %d scoreContribution: %f %n", targetIndex,
                    //          queryPosition, scoreContribution);
                    score += scoreContribution;
                } else {
                    score -= 0.25 * 0.25;
                }
                // penalize positions found in the target that are not found in the query
               // score -= Math.max(0, (queryPosition - lastQueryPosition - 1)) * 0.5f * 0.5f;
                lastQueryPosition = queryPosition;

            }
            results.enqueue(targetPosition, score);
            pg.lightUpdate();

        }

        pg.stop("brute force");
        return results;

    }

    private IntArrayList getActivePositions(Query query) {
        IntArrayList queryActivePositions;
        queryActivePositions = new IntArrayList();
        for (int queryIndex = 0; queryIndex < query.query.size(); queryIndex++) {
            final float rate = query.query.getFloat(queryIndex);
            if (rate == rate) {
                queryActivePositions.add(queryIndex);
            }
        }
        return queryActivePositions;
    }


    private FloatArrayList load(String inputFilename, String chromosomeSelected, char strandSelected) throws IOException {
        TSVReader reader = new TSVReader(new FileReader(inputFilename), '\t');
        FloatArrayList methylationRates = new FloatArrayList();
        int lastPosition = 0;
        boolean loadedSome = false;
        while (reader.hasNext()) {
            if (reader.isCommentLine()) {
                reader.skip();
            } else {
                reader.next();
                String chr = reader.getString();
                if (chr.equalsIgnoreCase(chromosomeSelected)) {
                    int position = reader.getInt();
                    char strand = reader.getString().charAt(0);
                    // ignore type of site.
                    reader.getString();
                    if (strandSelected == strand) {
                        int methylatedReadCount = reader.getInt();
                        int totalReadCount = reader.getInt();
                        float percentMethylation = ((float) methylatedReadCount) / ((float) totalReadCount);
                        while (lastPosition < position) {
                            methylationRates.add(Float.NaN);
                            lastPosition++;
                        }
                        // add element at index 'position'

                        methylationRates.add(percentMethylation);
                        lastPosition++;
                        loadedSome = true;
                    }
                } else {
                    if (loadedSome == false) {
                        reader.skip();
                    } else {
                        // past the chosen chromsome already.
                        break;
                    }
                }


            }
        }
        return methylationRates;
    }


}
