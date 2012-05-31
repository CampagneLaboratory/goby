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

import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.algorithmic.data.Read;
import edu.cornell.med.icb.goby.reads.Reads;
import edu.cornell.med.icb.goby.reads.ReadsReader;
import edu.cornell.med.icb.goby.util.IOUtil;
import it.unimi.dsi.bits.BitVector;
import it.unimi.dsi.bits.LongArrayBitVector;
import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectIterator;
import it.unimi.dsi.io.FastBufferedReader;
import it.unimi.dsi.io.LineIterator;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.util.IntHyperLogLogCounterArray;
import org.apache.commons.io.IOUtils;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

/**
 * Detect branch points in set of reads that differ in cardinality between two groups of samples.
 *
 * @author campagne
 *         Date: 9/30/11
 *         Time: 1:39 PM
 */
public class BranchingDetectorMode extends AbstractGobyMode {

    /**
     * The mode name.
     */
    private static final String MODE_NAME = "branch-point-detector";
    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Detect branch points in set of reads that differ in cardinality between two groups of samples.";


    private String outputFilename;

    private String[] groupAFilenames;
    private String[] groupBFilenames;
    private int kmerLength;
    private Int2IntMap codedGramToGramIndex=new Int2IntOpenHashMap();

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


        groupAFilenames = flatten(jsapResult.getString("group1"));
        groupBFilenames = flatten(jsapResult.getString("group2"));
        outputFilename = jsapResult.getString("output");
        kmerLength = jsapResult.getInt("kmer-length");

        return this;
    }

    private String[] flatten(String group1) {
        ObjectArrayList<String> groups = new ObjectArrayList<String>();
        try {

            LineIterator it = new LineIterator(new FastBufferedReader(new FileReader(group1)));
            while (it.hasNext()) {
                MutableString next = it.next();
                groups.add(next.toString());
            }
        } catch (FileNotFoundException e) {
            return new String[0];
        }
        return groups.toArray(new String[groups.size()]);
    }

    @Override
    public void execute() throws IOException {
        System.out.printf("Running with parameter k=%d %n",kmerLength);
        tallyKmers(groupAFilenames);
        System.out.printf("Found %d kmers ", frequencies.size());
        tallyKmers(groupBFilenames);
        prune(frequencies);
        // allocate counters:
        counters = new IntHyperLogLogCounterArray[2];
        for (int groupIndex : new int[]{0, 1}) {
            final long numDistinctElements = frequencies.size();
            final double estimateStandardDeviationSought = 0.1;
            final int numberOfCounters = frequencies.size();
            System.out.println("numberOfCounters " + numberOfCounters);
            counters[groupIndex] = new IntHyperLogLogCounterArray(numberOfCounters,
                    IntHyperLogLogCounterArray.registerSize(numDistinctElements),
                    IntHyperLogLogCounterArray.log2NumberOfRegisters(estimateStandardDeviationSought));

        }
       codedGramToGramIndex=new Int2IntOpenHashMap(frequencies.size());
        ObjectIterator<Int2IntMap.Entry> it = frequencies.int2IntEntrySet().fastIterator();
        int gramIndex = 0;

        while (it.hasNext()) {
            Int2IntMap.Entry kGramEntry = it.next();
            codedGramToGramIndex.put(kGramEntry.getIntKey(),  gramIndex++);
        }

        tallyBranchs(groupAFilenames, 0);
        tallyBranchs(groupBFilenames, 1);
        compareGroups(0, 1);
        System.out.println(counters);
    }

    private void compareGroups(int groupIndexA, int groupIndexB) {
        ObjectIterator<Int2IntMap.Entry> it = frequencies.int2IntEntrySet().fastIterator();
        int prunedCount = 0;
        while (it.hasNext()) {
            Int2IntMap.Entry kGramEntry = it.next();
            final int codedGram = kGramEntry.getIntKey();
            final int gramIndex = codedGramToGramIndex.get(codedGram);
            final double countGroup1 = counters[groupIndexA].count(gramIndex);
            final double countGroup2 = counters[groupIndexB].count(gramIndex);
            final double averageCount=(countGroup1+countGroup2)/2.0;
            if (averageCount>5 && Math.abs(countGroup1 - countGroup2)>(averageCount*0.5)) {
                System.out.printf("kgram %d %s has %g successors in group 1 and %g successors in group 2%n",
                        codedGram, decode(codedGram),
                        countGroup1, countGroup2);
            }
        }
    }


    private void tallyBranchs(String[] filenames, int groupIndex) {

        // initialize counters for each group:

        for (final String inputFilename : filenames) {
            ReadsReader reader = null;
            try {

                reader = new ReadsReader(inputFilename);
                final IntArrayFIFOQueue queue = new IntArrayFIFOQueue(kmerLength);
                for (final Reads.ReadEntry read : reader) {


                    queue.clear();
                    final byte[] sequence = read.getSequence().toByteArray();
                    for (int i = 0; i < sequence.length - kmerLength; i++) {
                        final int codedGram = codeGram(sequence, i, kmerLength);

                        if (queue.size() == kmerLength) { // queue at capacity, the first to dequeue
                            // is kmer bases in the past.
                            final int previousGram = queue.dequeue();
                            if (previousGram != -1 && codedGram != -1) {
                                int gramIndex = codedGramToGramIndex.get(previousGram);
                                counters[groupIndex].add(gramIndex, codedGram);
                                if (gramIndex == 0) {
                                    //System.out.printf("previous=%d next=%d groupIndex=%d %n",previousGram, codedGram, groupIndex);
                                }
                            }
                        }
                        queue.enqueue(codedGram);
                    }
                }
            } catch (IOException e) {

                IOUtil.closeQuietly(reader);

            } finally {
                IOUtil.closeQuietly(reader);
            }
        }
    }

    private void prune
            (Int2IntOpenHashMap
                     frequencies) {
        int pruningThreshold = 3;
        ObjectIterator<Int2IntMap.Entry> it = frequencies.int2IntEntrySet().fastIterator();
        int prunedCount = 0;
        while (it.hasNext()) {
            Int2IntMap.Entry kGramEntry = it.next();
            if (kGramEntry.getIntValue() <= pruningThreshold) {
                frequencies.remove(kGramEntry.getIntKey());
                //  System.out.printf("Pruned %d frequency %d%n", kGramEntry.getIntKey(), kGramEntry.getIntValue());
                ++prunedCount;
            }
        }
        System.out.printf("Pruned %d kGrams.%n", prunedCount);
    }

    private void tallyKmers
            (String[] filenames) {
        for (final String inputFilename : filenames) {
            ReadsReader reader = null;
            try {
                reader = new ReadsReader(inputFilename);
                for (Reads.ReadEntry read : reader) {
                    final byte[] sequence = read.getSequence().toByteArray();
                    for (int i = 0; i < sequence.length - kmerLength; i++) {
                        final int codedGram = codeGram(sequence, i, kmerLength);
                        incrementFrequency(codedGram);

                    }
                }
            } catch (IOException e) {

                IOUtil.closeQuietly(reader);

            }
        }
    }

    Int2IntOpenHashMap frequencies = new Int2IntOpenHashMap();
    it.unimi.dsi.util.IntHyperLogLogCounterArray[] counters;

    private void incrementFrequency
            (
                    int codedGram) {
        if (codedGram == -1) {
            return;
        }

        int freq = frequencies.get(codedGram);
        frequencies.put(codedGram, freq + 1);
    }

    private BitVector vector = LongArrayBitVector.getInstance();

    private int codeGram
            (
                    final byte[] sequence,
                    final int i,
                    final int kmerLength) {
        assert kmerLength <= 16 : "The kmer must fit in a 32bit integer.";
        vector.clear();
        for (int j = 0; j < kmerLength; ++j) {
            switch (sequence[i + j]) {
                case 'A':
                    vector.add(0);
                    vector.add(0);
                    //     System.out.println('A');
                    break;
                case 'C':
                    vector.add(0);
                    vector.add(1);
                    //     System.out.println('C');
                    break;

                case 'T':
                    vector.add(1);
                    vector.add(0);
                    //   System.out.println('T');
                    break;
                case 'G':
                    vector.add(1);
                    vector.add(1);
                    //  System.out.println('G');
                    break;
                case 'N':
                    // coded gram is -1 if the sequence of the gram contains an N at any position
                    return -1;
            }
        }
        int aLong = (int) vector.getLong(0, kmerLength * 2 - 1);
        //  System.out.println(aLong);
        return aLong;
    }

    MutableString decodedGram = new MutableString();
    long[] value = new long[1];

    private String decode(int codedGram) {
        decodedGram.setLength(0);
        value[0] = codedGram;
        LongArrayBitVector v = LongArrayBitVector.wrap(value);
        int j = kmerLength * 2 + 1;
        while (--j > 0) {
            boolean first = v.getBoolean(j);
            boolean second = v.getBoolean(j);
            if (first) {
                if (second) {
                    decodedGram.append('G');
                } else {
                    decodedGram.append('T');
                }
            } else {
                if (second) {
                    decodedGram.append('C');
                } else {
                    decodedGram.append('A');
                }
            }
        }
        return decodedGram.toString();
    }

    public static void main
            (
                    final String[] args) throws IOException, JSAPException {
        new BranchingDetectorMode().configure(args).execute();
    }
}