/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
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

import edu.cornell.med.icb.parsers.FastaParser;
import edu.mssm.crover.cli.CLI;
import it.unimi.dsi.lang.MutableString;

import java.io.*;
import java.util.Random;
import java.util.zip.GZIPInputStream;

/**
 * A utility to plant indels in a genome.
 *
 * @author Fabien Campagne
 *         Date: 5/31/11
 *         Time: 2:23 PM
 */
public class PlantIndels {
    private String refChoice;
    private String outputFilename;
    private boolean insert = false;
    private int numIndels;

    public static void main(final String args[]) throws IOException {
        final String fastaReference = CLI.getOption(args, "-r", null);
        final String outputFilename = CLI.getOption(args, "-o", "out.fq");
        final String regionIndelTruth = CLI.getOption(args, "-t", "indel-truth.tsv");
        // reference sequence to use
        final String refChoice = CLI.getOption(args, "-c", "22");
        int from = CLI.getIntOption(args, "-s", 0);
        int to = CLI.getIntOption(args, "-e", 0);

        if (to < from) {
            System.err.println("argument to -e must be larger than argument to -s");
            System.exit(1);
        }
        final int indelLength = CLI.getIntOption(args, "-l", 3);
        PlantIndels processor = new PlantIndels();
        processor.outputFilename = outputFilename;
        processor.insert = !CLI.isKeywordGiven(args, "--read-deletion");
        processor.numIndels=CLI.getIntOption(args,"-n",5);
        processor.process(refChoice, fastaReference, from, to, indelLength, regionIndelTruth);


    }

    public void process(String refChoice, String fastaReference, int from, int to, int indelLength, String regionIndelTruth) throws IOException {
        this.refChoice = refChoice;

        FastaParser parser = new FastaParser(fastaReference.endsWith(".gz") ?

                new InputStreamReader(new GZIPInputStream(new FileInputStream(fastaReference))) :
                new FileReader(fastaReference)
        );
        PrintWriter writer = null;
        PrintWriter indelTruthWriter = null;
        try {

            writer = new PrintWriter(new FileWriter(outputFilename));
            indelTruthWriter = new PrintWriter(new FileWriter(regionIndelTruth));
            MutableString bases = new MutableString();
            MutableString description = new MutableString();
            MutableString accession = new MutableString();
            while (parser.hasNext()) {

                parser.next(description, bases);

                FastaParser.guessAccessionCode(description, accession);
                if (accession.equalsIgnoreCase(refChoice)) {

                    // from and to are zero-based
                    System.out.println("matching " + description);
                    writer.append(String.format(">%s%n", accession));
                    writer.append(bases.substring(0, from));

                    writer.append(plantIndels(bases.substring(from, to), from, indelTruthWriter, indelLength));
                    writer.append(bases.substring(to, bases.length()));
                }
            }
        } finally {
            if (writer != null) {
                writer.close();
            }
            if (indelTruthWriter != null) {
                indelTruthWriter.close();
            }
        }
    }

    /**
     * @param lo lower limit of range
     * @param hi upper limit of range
     * @return a random integer in the range <STRONG>lo</STRONG>,
     *         <STRONG>lo</STRONG>+1, ... ,<STRONG>hi</STRONG>
     */
    public int choose(Random random, final int lo, final int hi) {
        final int randomValue = (int) ((long) lo + (long) ((1L + (long) hi - (long) lo) * random.nextDouble()));
        assert randomValue <= hi;
        assert randomValue >= lo;


        return randomValue;
    }

    private MutableString plantIndels(MutableString segmentBases, int from, PrintWriter indelTruthWriter, int indelLength) {
        Random random = new Random();
        int minStart = 0;
        int positionOffset = 0;
        // we plant 5 indels in the region:
        for (int i = 0; i < numIndels; i++) {
            int indelStartPosition = chooseNext(minStart, segmentBases, indelLength, random);

            if (insert) {
                // insert bases in reference, they will end up in the read.
                String indelBases = "ACTG".substring(0, indelLength);
                segmentBases.insert(indelStartPosition, indelBases);
                indelTruthWriter.append(String.format("READ_INSERTION\t%d\t%d\t%s\t%s%n",
                        indelStartPosition + from, indelLength , indelBases,
                        context(segmentBases, indelStartPosition, indelLength)));
                positionOffset += indelLength;
            } else {
                // remove bases from reference, they will appear as read deletion
                MutableString indelBases = segmentBases.substring(indelStartPosition, indelStartPosition + indelLength);
           //     System.out.printf("before: %s%n", context(segmentBases, indelStartPosition, 0));
                CharSequence end = segmentBases.substring(indelStartPosition + indelLength).copy();
                segmentBases.setLength(indelStartPosition);
                segmentBases.append(end);

                indelTruthWriter.append(String.format("READ_DELETION\t%d\t%d\t%s\t%s%n",
                        indelStartPosition + from , indelLength, indelBases,
                        context(segmentBases, indelStartPosition, 0)));
                positionOffset -= indelLength;

               // System.out.printf("after:  %s%n", context(segmentBases, indelStartPosition, 0));
            }
            minStart = indelStartPosition + indelLength;
            if (minStart > segmentBases.length() - indelLength) break;
        }
        indelTruthWriter.flush();
        return segmentBases;
    }

    private CharSequence context(MutableString segmentBases, int indelStartPosition, int indelLength) {
        return segmentBases.subSequence(Math.max(0, indelStartPosition - 5),
                Math.min(segmentBases.length(),                  indelLength + indelStartPosition + 5));
    }

    private int chooseNext(int minStart, MutableString segmentBases, int indelLength, Random random) {
        int start;
        do {
            start = choose(random, 0, Math.max(0, segmentBases.length() - 1 - indelLength));
        } while (start < minStart);
        return start;
    }
}
