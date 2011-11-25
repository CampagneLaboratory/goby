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

package edu.cornell.med.icb.goby.modes;

import com.google.protobuf.ByteString;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.reads.Reads;
import edu.cornell.med.icb.goby.reads.ReadsReader;
import edu.cornell.med.icb.goby.reads.ReadsWriter;
import it.unimi.dsi.bits.BitVector;
import it.unimi.dsi.bits.LongArrayBitVector;
import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.io.LineIterator;
import org.apache.log4j.Logger;

import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;

/**
 * Trims adapter sequences from reads.
 *
 * @author Fabien Campagne
 *         Date: June 6 2011
 *         Time: 9:40 PM
 */
public class TrimMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "trim";
    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Trims reads to remove adapter sequences.";

    private String inputFilename;
    private String outputFilename;
    private static final Logger LOG = Logger.getLogger(TrimMode.class);
    private byte[] buffer = new byte[10000];
    private String adapterFilename;
    private boolean complementAdapters;
    private int numTrimmedLeft;
    private int numTrimmedRight;
    private int numContained;
    private double numSequencesInInput;
    private int numTrimmed;
    private int minRightLength = 0;
    private int minLeftLength = 0;


    /**
     * {@inheritDoc}
     */
    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }


    /**
     * Configure.
     *
     * @param args command line arguments
     * @return this object for chaining
     * @throws java.io.IOException error parsing
     * @throws com.martiansoftware.jsap.JSAPException
     *                             error parsing
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args)
            throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        inputFilename = jsapResult.getString("input");
        outputFilename = jsapResult.getString("output");
        /**
         * File with one line per adapter sequence.
         */
        adapterFilename = jsapResult.getString("adapters");
        complementAdapters = jsapResult.getBoolean("complement");
        minLeftLength = jsapResult.getInt("min-left-length");
        minRightLength = jsapResult.getInt("min-right-length");
        return this;
    }

    BitVector vector = LongArrayBitVector.getInstance();
    double maxCount;

    @Override
    public void execute() throws IOException {

        ReadsReader reader = null;
        ReadsWriter writer = null;
        ProgressLogger progress = new ProgressLogger(LOG);
        try {

            reader = new ReadsReader(inputFilename);
            LineIterator lines = new LineIterator(new FileReader(adapterFilename));
            ObjectArrayList<MutableString> adapterList = new ObjectArrayList<MutableString>();
            while (lines.hasNext()) {
                String next = lines.nextLine();
                adapterList.add(new MutableString(next));
            }
            MutableString[] adapters;
            if (complementAdapters) {
                adapters = addComplementAdapters(adapterList);
            } else {
                adapters = adapterList.toArray(new MutableString[adapterList.size()]);
            }
            progress.start();
            writer = new ReadsWriter(new FileOutputStream(outputFilename));

            ByteArrayList newQualScores = new ByteArrayList();
            ByteArrayList newPairQualScores = new ByteArrayList();
            MutableString sequence = new MutableString();
            MutableString sequencePair = new MutableString();

            for (final Reads.ReadEntry entry : reader) {
                //      observe(counters, entry.getSequence(), entry.getReadIndex());
                ReadsReader.decodeSequence(entry, sequence);

                ByteString qualityScores = entry.getQualityScores();
                newQualScores.clear();
                MutableString seq1 = trim(adapters, newQualScores, sequence, qualityScores);
                MutableString pairSeq = null;

                numSequencesInInput++;
                if (entry.hasSequencePair()) {
                    newPairQualScores.clear();

                    ReadsReader.decodeSequence(entry, sequencePair, true);

                    ByteString pairQualityScores = entry.getQualityScoresPair();
                    pairSeq = trim(adapters, newPairQualScores, sequencePair, pairQualityScores);
                    numSequencesInInput++;
                }

                //    System.out.printf(">seq%n%s%n", c);
                Reads.ReadEntry.Builder builder = Reads.ReadEntry.newBuilder();
                builder = builder.mergeFrom(entry).setSequence(ReadsWriter.encodeSequence(seq1, buffer)).setReadLength(seq1.length());
                if (sequence.length() != seq1.length()) {
                    numTrimmed++;
                    final byte[] bytes1 = newQualScores.toByteArray();
                    builder = builder.setQualityScores(ByteString.copyFrom(bytes1));
                    assert builder.getQualityScores().size() == builder.getSequence().size() : "sequence length and quality scores must match.";
                }

                if (entry.hasSequencePair()) {
                    builder = builder.mergeFrom(entry)
                            .setSequencePair(ReadsWriter.encodeSequence(pairSeq, buffer))
                            .setReadLength(pairSeq.length());

                    if (sequencePair.length() != pairSeq.length()) {
                        numTrimmed++;
                        builder = builder.setQualityScoresPair(ByteString.copyFrom(newPairQualScores.toByteArray()));
                        assert builder.getQualityScoresPair().size() == builder.getSequencePair().size() : "sequence length and quality scores must match.";

                    }
                }

                writer.appendEntry(builder);
                progress.lightUpdate();
            }
            progress.stop();

            final int numSequencesTrimmed = numTrimmed;
            double percent = 100d * numSequencesTrimmed;
            percent /= numSequencesInInput;
            System.out.printf("Number of reads trimmed %d (%g %% of input sequences), including: %n" +
                    "left: %d (%g%%)%n" +
                    "right: %d (%g%%), %n" +
                    "fully contained: %d (%g%%)%n",
                    numSequencesTrimmed, percent, numTrimmedLeft,
                    percent(numTrimmedLeft, numSequencesTrimmed),
                    numTrimmedRight,
                    percent(numTrimmedRight, numSequencesTrimmed),
                    numContained,
                    percent(numContained, numSequencesTrimmed));
            System.out.flush();
        } finally {
            if (writer != null) {
                writer.close();
            }

        }

        progress.stop();
    }

    private double percent(double a, double b) {
        return a / b * 100;
    }

    protected MutableString trim(MutableString[] adapters, ByteArrayList newQualScores, MutableString sequence, ByteString qualityScores) {
        final int length = sequence.length();
        MutableString a = contains(length, sequence, qualityScores, newQualScores, adapters);
        MutableString b = trimLeft(length, a, qualityScores, newQualScores, adapters);
        return trimRight(length, b, qualityScores, newQualScores, adapters);
    }

    protected void convert(ByteString bytes, MutableString sequence) {
        int length = bytes.size();
        sequence.setLength(length);

        for (int pos = 0; pos < length; pos++) {
            sequence.charAt(pos, (char) bytes.byteAt(pos));
        }

    }

    protected MutableString[] addComplementAdapters(ObjectArrayList<MutableString> adapterList) {
        ObjectArrayList<MutableString> result = new ObjectArrayList<MutableString>();
        for (MutableString adapter : adapterList) {
            result.add(adapter);
            result.add(complement(adapter));
        }
        return result.toArray(new MutableString[result.size()]);
    }

    private MutableString complement(MutableString input) {
        MutableString result = new MutableString();
        result.setLength(input.length());
        // return the complement of input:
        for (int i = 0; i < input.length(); i++) {
            char base = input.charAt(i);
            switch (base) {
                case 'A':
                    base = 'T';
                    break;
                case 'C':
                    base = 'G';
                    break;
                case 'G':
                    base = 'C';
                    break;
                case 'T':
                    base = 'A';
                    break;
            }
            result.charAt(i, base);
        }
        return result;
    }


    protected MutableString trimRight(int length, MutableString sequence,
                                      ByteString qualityScores,
                                      ByteArrayList newQualScores,
                                      MutableString[] adapters) {

        int currentLength = sequence.length();
        for (MutableString adapter : adapters) {
            final int adaptLength = adapter.length();

            for (int j = minRightLength; j < adaptLength; j++) {
                if (sequence.endsWith(adapter.subSequence(j, adaptLength))) {
                    final int trimedLength = adaptLength - j;
                    if (trimedLength > 10) {
                        System.out.printf("%d bases matching right %s %s %n", trimedLength, sequence, adapter);
                    }
                    if (currentLength == length) {
                        copy(qualityScores, newQualScores);
                    }
                    if (!qualityScores.isEmpty()) {
                        newQualScores.removeElements(currentLength - trimedLength, Math.min(currentLength + 1, newQualScores.size()));
                    }
                    numTrimmedRight++;
                    return sequence.substring(0, currentLength - trimedLength);


                }

            }
        }
        return sequence;
    }

    protected MutableString trimLeft(int length, MutableString sequence, ByteString qualityScores, ByteArrayList newQualScores, MutableString[] adapters) {
        int currentLength = sequence.length();

        for (MutableString adapter : adapters) {
            final int adaptLength = adapter.length();
            for (int j = adaptLength; j >= minLeftLength; --j) {
                if (sequence.startsWith(adapter.subSequence(0, j))) {
                    final int trimedLength = j;
                    if (trimedLength > 10) {
                        System.out.printf("%d bases matching left %s %s %n", trimedLength, sequence, adapter);
                    }
                    if (currentLength == length) { // previously unchanged, we need to copy quality score to the list representation for editing.
                        copy(qualityScores, newQualScores);
                    }
                    newQualScores.removeElements(0, trimedLength);
                    numTrimmedLeft++;
                    return sequence.substring(trimedLength, currentLength);


                }
            }
        }
        return sequence;
    }

    protected MutableString contains(int length, MutableString sequence, ByteString qualityScores, ByteArrayList newQualScores, MutableString[] adapters) {

        for (MutableString adapter : adapters) {
            final int index = sequence.indexOf(adapter);
            if (index >= 0) {
                System.out.printf("adapter %s contained entirely in sequence %s%n", adapter, sequence);
                copy(qualityScores, newQualScores);

                final int end = adapter.length() + index;
                newQualScores.removeElements(index, end);
                numContained++;
                return sequence.delete(index, end);
            }

        }
        return sequence;
    }

    private void copy(ByteString qualityScores, ByteArrayList newQualScores) {
        for (byte qual : qualityScores.toByteArray()) {
            newQualScores.add(qual);
        }
    }


    public static void main
            (
                    final String[] args) throws IOException, JSAPException {
        new TrimMode().configure(args).execute();
    }

    public void setMinLengthLeft(int i) {
        this.minLeftLength = i;
    }

    public void setMinLengthRight(int i) {
        this.minRightLength = i;
    }
}
