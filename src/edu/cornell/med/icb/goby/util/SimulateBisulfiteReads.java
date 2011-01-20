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

import edu.mssm.crover.cli.CLI;
import edu.cornell.med.icb.parsers.FastaParser;
import edu.rit.mp.buf.DoubleArrayBuf;

import java.io.*;
import java.util.Random;
import java.text.Format;

import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.fastutil.doubles.DoubleList;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.doubles.DoubleIterator;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;
import it.unimi.dsi.io.LineIterator;
import it.unimi.dsi.io.FastBufferedReader;
import org.apache.commons.io.IOUtils;

/**
 * @author Fabien Campagne
 *         Date: Jan 20, 2011
 *         Time: 1:22:45 PM
 */
public class SimulateBisulfiteReads {
    private int readLength;
    private Random random;
    private String refChoice;
    private String outputFilename;
    private Int2DoubleMap methylationForward;
    private Int2DoubleMap methylationReverse;
    private String regionTrueRates;

    public SimulateBisulfiteReads() {
        random = new Random();

    }

    public static void main(String[] args) throws IOException {

        String fastaReference = CLI.getOption(args, "-r", null);
        String outputFilename = CLI.getOption(args, "-o", "out.fa");
        String regionTrueRates = CLI.getOption(args, "-t", "true-methylation.tsv");
        // reference sequence to use
        String refChoice = CLI.getOption(args, "-c", "22");
        int from = CLI.getIntOption(args, "-s", 0);
        int to = CLI.getIntOption(args, "-e", 0);
        if (to < from) {
            System.err.println("argument to -e must be larger than argument to -s");
            System.exit(1);
        }
        int readLength = CLI.getIntOption(args, "-l", 40);
        String methylationRateFilename = CLI.getOption(args, "-m", "methylation-rates.tsv");
        SimulateBisulfiteReads processor = new SimulateBisulfiteReads();
        processor.readLength = readLength;
        processor.outputFilename = outputFilename;
        processor.regionTrueRates = regionTrueRates;
        processor.process(refChoice, fastaReference, from, to, methylationRateFilename);

    }

    private void process(String refChoice, String fastaReference, int from, int to, String methylationRateFilename) throws IOException {
        this.refChoice = refChoice;
        DoubleList methylationRates = load(methylationRateFilename);
        FastaParser parser = new FastaParser(new FileReader(fastaReference));
        MutableString bases = new MutableString();
        MutableString description = new MutableString();
        MutableString accession = new MutableString();
        while (parser.hasNext()) {

            parser.next(description, bases);

            FastaParser.guessAccessionCode(description, accession);
            if (accession.equalsIgnoreCase(refChoice)) {
                process(methylationRates, bases.substring(from, to));
            }
        }
    }

    private void process(DoubleList methylationRates, CharSequence segmentBases) throws IOException {
        PrintWriter writer = new PrintWriter(new FileWriter(outputFilename));
        PrintWriter trueRateWriter = new PrintWriter(new FileWriter(regionTrueRates));
        System.out.printf("segmentBases.length: %d %n", segmentBases.length());
        methylationForward = new Int2DoubleOpenHashMap();
        methylationReverse = new Int2DoubleOpenHashMap();

        final int segmentLength = segmentBases.length();
        DoubleIterator it = methylationRates.iterator();
        // prepare methylation rates for each C position of the segment, forward strand:
        for (int i = 0; i < segmentBases.length(); i++) {
            if (segmentBases.charAt(i) == 'C') {
                if (!it.hasNext()) {
                    it = methylationRates.iterator();
                }
                methylationForward.put(i, it.nextDouble());
                trueRateWriter.printf("%d\t%g\t+1%n", i, methylationForward.get(i));
            }
        }
        it = methylationRates.iterator();
        CharSequence reverseStrandSegment = reverseComplement(segmentBases);
        for (int i = reverseStrandSegment.length()-1; i >=0; i--) {
            if (reverseStrandSegment.charAt(i) == 'C') {
                if (!it.hasNext()) {
                    it = methylationRates.iterator();
                }
                final int positionInSegment =  i;
                methylationReverse.put(positionInSegment, it.nextDouble());
                trueRateWriter.printf("%d\t%g\t-1%n", positionInSegment, methylationReverse.get(positionInSegment));
            }
        }
        for (int repeatCount = 0; repeatCount < 10000; repeatCount++) {
            int startReadPosition = choose(0, segmentBases.length() - 1 - readLength);
            boolean matchedReverseStrand = random.nextBoolean();
            final CharSequence selectedReadRegion = segmentBases.subSequence(startReadPosition, startReadPosition + readLength);
            CharSequence readBases = matchedReverseStrand ? reverseComplement(selectedReadRegion) : selectedReadRegion;


            MutableString sequence = new MutableString();
            MutableString log = new MutableString();

            for (int i = 0; i < readLength; i++) {
                final int basePositionInSegment = i;
                char base = readBases.charAt(basePositionInSegment);
                if (base == 'C') {

                    boolean isBaseMethylated = random.nextDouble() < getMethylationRateAtPosition(segmentLength,
                            matchedReverseStrand,
                            i, startReadPosition);
                    if (!isBaseMethylated) {
                        // bases that are not methylated are changed to T through the bisulfite and PCR conversion steps
                        base = 'T';

                    } else {
                        // bases that are methylated are protected and stay C on the forward strand. They would also
                        // be seen as G on the opposite strand if the sequencing protocol did not respect strandness
                        log.append("met:");
                        log.append(i + startReadPosition);
                        log.append(' ');
                    }


                }
                sequence.append(base);
            }
            writer.printf(">%d reference: %s startPosition: %d strand: %s %s%n%s%n", repeatCount, refChoice,
                    startReadPosition, matchedReverseStrand ? "-1" : "+1", log, sequence);


        }
        writer.close();
        trueRateWriter.close();
    }

    private double getMethylationRateAtPosition(int segmentLength,
                                                boolean matchedReverseStrand,
                                                int positionInRead,
                                                int startReadPosition) {
        return matchedReverseStrand ? methylationReverse.get(segmentLength - (startReadPosition + readLength - positionInRead)) :
                methylationForward.get(positionInRead + startReadPosition);
    }

    private CharSequence reverseComplement(CharSequence selectedReadRegion) {
        MutableString s = new MutableString();
        s.setLength(selectedReadRegion.length());
        int j = 0;
        for (int i = s.length() - 1; i >= 0; i--) {
            char base = selectedReadRegion.charAt(i);
            switch (base) {
                case 'A':
                    base = 'T';
                    break;
                case 'C':
                    base = 'G';
                    break;
                case 'T':
                    base = 'C';
                    break;
                case 'G':
                    base = 'C';
                    break;
                default:
                    base = 'N';
                    break;
            }
            s.charAt(j++, base);

        }
        return s;
    }

    private DoubleList load(String methylationRateFilename) throws FileNotFoundException {

        DoubleList result = new DoubleArrayList();
        File toOpen = new File(methylationRateFilename);
        if (!toOpen.exists() || !toOpen.canRead()) {
            System.err.printf("methylation file cannot be read %s", methylationRateFilename);
            System.exit(10);
        }
        LineIterator it = new LineIterator(new FastBufferedReader(new FileReader(methylationRateFilename)));

        while (it.hasNext()) {
            MutableString mutableString = it.next();

            result.add(Double.valueOf(mutableString.toString()));
        }
        return result;
    }

    /**
     * @param lo lower limit of range
     * @param hi upper limit of range
     * @return a random integer in the range <STRONG>lo</STRONG>,
     *         <STRONG>lo</STRONG>+1, ... ,<STRONG>hi</STRONG>
     */
    public int choose(final int lo, final int hi) {
        final int random = (int) ((long) lo + (long) ((1L + (long) hi - (long) lo) * this.random.nextDouble()));
        assert random <= hi;
        assert random >= lo;


        return random;
    }


}
