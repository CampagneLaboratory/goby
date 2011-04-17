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

import edu.cornell.med.icb.goby.reads.QualityEncoding;
import edu.cornell.med.icb.parsers.FastaParser;
import edu.mssm.crover.cli.CLI;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.doubles.DoubleIterator;
import it.unimi.dsi.fastutil.doubles.DoubleList;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;
import it.unimi.dsi.io.FastBufferedReader;
import it.unimi.dsi.io.LineIterator;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.FilenameUtils;

import java.io.*;
import java.util.Random;

/**
 * Creates fastq files with simulated reads. Simulation can create bisulfite treated reads (or not treated)
 * with mutations on C bases at some rate.
 *
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
    private boolean bisulfiteTreatment;
    private long seed = 232424434L;
    private int numRepeats;

    public SimulateBisulfiteReads() {
        random = new Random(seed);

    }

    public static void main(String[] args) throws IOException {

        String fastaReference = CLI.getOption(args, "-r", null);
        String outputFilename = CLI.getOption(args, "-o", "out.fq");
        String regionTrueRates = CLI.getOption(args, "-t", "true-methylation.tsv");
        // reference sequence to use
        String refChoice = CLI.getOption(args, "-c", "22");
        int from = CLI.getIntOption(args, "-s", 0);
        int to = CLI.getIntOption(args, "-e", 0);
        if (to < from) {
            System.err.println("argument to -e must be larger than argument to -s");
            System.exit(1);
        }
        int readLength = CLI.getIntOption(args, "-l", 50);
        String methylationRateFilename = CLI.getOption(args, "-m", "methylation-rates.tsv");
        SimulateBisulfiteReads processor = new SimulateBisulfiteReads();
        final boolean bisulfite = CLI.isKeywordGiven(args, "--bisulfite");
        if (bisulfite) {
            System.out.println("Bisulfite treatment activated.");
        }
        processor.bisulfiteTreatment = bisulfite;
        processor.readLength = readLength;
        processor.outputFilename = outputFilename;
        processor.regionTrueRates = regionTrueRates;
        processor.numRepeats = CLI.getIntOption(args, "-n", 10000);

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
                process(methylationRates, bases.substring(from, to), from);
            }
        }
    }

    private void process(DoubleList methylationRates, CharSequence segmentBases, int from) throws IOException {
        PrintWriter writer = new PrintWriter(new FileWriter(outputFilename));
        final String trueRateFilename = FilenameUtils.removeExtension(outputFilename) + "-true-methylation.tsv";
        PrintWriter trueRateWriter = new PrintWriter(new FileWriter(trueRateFilename));
        // System.out.printf("segmentBases.length: %d %n", segmentBases.length());
        methylationForward = new Int2DoubleOpenHashMap();
        methylationReverse = new Int2DoubleOpenHashMap();
        trueRateWriter.printf("%s\t%s\ts\t%s\t%s%n", "index", "methylationRate", "chromosome", "position");
        final int segmentLength = segmentBases.length();
        DoubleIterator it = methylationRates.iterator();
        // prepare methylation rates for each C position of the segment, forward strand:
        for (int i = 0; i < segmentBases.length(); i++) {
            if (segmentBases.charAt(i) == 'C') {
                if (!it.hasNext()) {
                    it = methylationRates.iterator();
                }
                final double value = it.nextDouble();
                methylationForward.put(i, value);
                assert getMethylationRateAtPosition(segmentLength, false, i, 0) == value;
                trueRateWriter.printf("%d\t%g\t+1\t%s\t%d%n", i, methylationForward.get(i), refChoice, i + from);
            }
        }
        it = methylationRates.iterator();
        CharSequence reverseStrandSegment = reverseComplement(segmentBases);
        // same for reverse strand:
        for (int i = reverseStrandSegment.length() - 1; i >= 0; i--) {
            if (reverseStrandSegment.charAt(i) == 'C') {
                if (!it.hasNext()) {
                    it = methylationRates.iterator();
                }
                final double value = it.nextDouble();
                methylationReverse.put(i, value);
                final double getterValue = getMethylationRateAtPosition(segmentLength, true, i, 0);
                assert getterValue == value : "getter must work for reverse strand";
                trueRateWriter.printf("%d\t%g\t-1\t%s\t%d%n", i, methylationReverse.get(i),
                        refChoice, i + from);
            }
        }
        for (int repeatCount = 0; repeatCount < numRepeats; repeatCount++) {
            int startReadPosition = choose(0, segmentBases.length() - 1 - readLength);
            boolean matchedReverseStrand = random.nextBoolean();
            final CharSequence selectedReadRegion = segmentBases.subSequence(startReadPosition, startReadPosition + readLength);
            CharSequence readBases = matchedReverseStrand ? reverseComplement(selectedReadRegion) : selectedReadRegion;


            MutableString sequenceInitial = new MutableString();
            MutableString sequenceTreated = new MutableString();
            MutableString log = new MutableString();

            for (int i = 0; i < readLength; i++) {
                final int basePositionInSegment = i;
                char base = readBases.charAt(basePositionInSegment);
                sequenceInitial.append(base);
                if (base == 'C') {

                    final boolean b = random.nextDouble() < getMethylationRateAtPosition(segmentLength,
                            matchedReverseStrand,
                            i, startReadPosition);
                    boolean isBaseMethylated = b;
                    if (!isBaseMethylated) {
                        // bases that are not methylated are changed to T through the bisulfite and PCR conversion steps
                        if (bisulfiteTreatment) {
                            base = 'T';
                        }

                    } else {
                        if (!bisulfiteTreatment) {
                            // mutate base to G
                            // introduce mutation C -> G
                            base = 'G';
                        }
                        // bases that are methylated are protected and stay C on the forward strand. They would also
                        // be seen as G on the opposite strand if the sequencing protocol did not respect strandness
                        log.append(bisulfiteTreatment ? "met: " : "mut: ");
                        log.append(i + startReadPosition + from);
                        log.append(' ');

                        log.append("read-index: ");
                        log.append(matchedReverseStrand ? readLength - i : i);
                        log.append(' ');


                    }
                }
                sequenceTreated.append(base);
            }
            MutableString coveredPositions = new MutableString();
            MutableString qualityScores = new MutableString();
            for (int i = 0; i < readLength; i++) {
                final char c = QualityEncoding.ILLUMINA.phredQualityScoreToAsciiEncoding((byte) 40);
                qualityScores.append(c);

            }
            for (int i = startReadPosition + from; i < startReadPosition + from + readLength; i++) {
                coveredPositions.append(i);
                coveredPositions.append(" ");

            }
            writer.printf("@%d reference: %s startPosition: %d strand: %s %s %s%n%s%n+%n%s%n",
                    repeatCount,
                    refChoice,

                    startReadPosition, matchedReverseStrand ? "-1" : "+1", log,
                    coveredPositions, sequenceTreated,
                    qualityScores);
        }
        writer.close();
        trueRateWriter.close();
    }

    private double getMethylationRateAtPosition(int segmentLength,
                                                boolean matchedReverseStrand,
                                                int positionInRead,
                                                int startReadPosition) {
        return matchedReverseStrand ?
                methylationReverse.get(positionInRead + startReadPosition) :
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
                    base = 'A';
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
