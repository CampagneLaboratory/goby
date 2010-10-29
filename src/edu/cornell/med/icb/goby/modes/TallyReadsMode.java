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

package edu.cornell.med.icb.goby.modes;

import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.reads.CompressedRead;
import edu.cornell.med.icb.goby.reads.ReadSet;
import edu.cornell.med.icb.goby.reads.Reads;
import edu.cornell.med.icb.goby.reads.ReadsReader;
import edu.cornell.med.icb.goby.reads.SequenceDigests;
import it.unimi.dsi.bits.BitVector;
import it.unimi.dsi.bits.LongArrayBitVector;
import it.unimi.dsi.fastutil.ints.Int2ByteMap;
import it.unimi.dsi.fastutil.ints.Int2ByteOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import it.unimi.dsi.io.OutputBitStream;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.logging.ProgressLogger;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;

/**
 * Tally the number of times sequences appear in a set of read files. Exact sequence comparison
 * is performed.
 *
 * @author Fabien Campagne
 *         Date: May 4 2009
 *         Time: 12:28 PM
 */
public class TallyReadsMode extends AbstractGobyMode {
    private String inputFilename;
    private String outputBasename;

    /**
     * The mode name.
     */
    private static final String MODE_NAME = "tally-reads";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Tally the number of times sequences appear "
            + "in a set of read files. Exact sequence comparison is performed.";

    private boolean colorSpace;
    private final int MAX_PROCESS_READS = Integer.MAX_VALUE;

    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }

    /**
     * Configure.
     *
     * @param args command line arguments
     * @return this object for chaining
     * @throws IOException error parsing
     * @throws JSAPException error parsing
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args)
            throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        inputFilename = jsapResult.getString("input");
        outputBasename = jsapResult.getString("output");
        colorSpace = jsapResult.getBoolean("color-space");
        return this;
    }


    @Override
    public void execute() throws IOException {
        final MutableString sequence = new MutableString();

        final ProgressLogger progress = new ProgressLogger();
        progress.start("first pass: starting to collect inspect hash");
        progress.displayFreeMemory = true;



        int readLength = 0;
        Int2ByteMap readRedundant = new Int2ByteOpenHashMap();
        int numReads = 0;
        SequenceDigests sd = null;
        {


            final ReadsReader readsReader = new ReadsReader(new FileInputStream(inputFilename));
            for (final Reads.ReadEntry readEntry : readsReader) {
                readLength = readEntry.getReadLength();
                if (sd == null) {
                    sd = new SequenceDigests(readLength, true);
                }

                final int hashCode = sd.digest(readEntry.getSequence(), 0, readLength);
                byte occ = readRedundant.get(hashCode);
                if (occ < Byte.MAX_VALUE - 1) {
                    occ += 1;
                }
                readRedundant.put(hashCode, occ);
                numReads++;
                progress.lightUpdate();

                if (numReads > MAX_PROCESS_READS) {
                    break;
                }
            }
            readsReader.close();
            System.gc();
        }
        progress.stop("first pass finished.");
        assert sd != null : "at least one read must have been processed.";

        int redundantHash = 0;
        final IntSet inspectHashcodes = new IntOpenHashSet();
        for (final int key : readRedundant.keySet()) {
            if (readRedundant.get(key) >= 2) {
                inspectHashcodes.add(key);
                redundantHash++;
            }

        }

        System.out.printf("Found %d redundant hashCodes. %n", redundantHash);
        readRedundant = null;
        progress.expectedUpdates = numReads;
        progress.start("second pass: actual read redundancy evaluation.");
        final int numberOfReads = numReads;
        numReads = 0;
        // readRedundant contains a lower bound on the number of times the read has been seen. Since we keyed
        // on hashcode, hashcode collisions make the set of redundant reads somewhat unreliable. This is not too bad because
        // we will now calculate redundancy exactly for those reads that passed the initial screen.

        {
            byte[] byteBuffer = new byte[readLength];
            // set initial capacity according to first pass:
            final Object2IntMap<CompressedRead> tallies = new Object2IntOpenHashMap<CompressedRead>(redundantHash + 1);
            tallies.defaultReturnValue(0);
            //   IntSet otherReadIndices = new IntOpenHashSet(numReads+1);
            final BitVector otherReadIndices = LongArrayBitVector.ofLength(numberOfReads);

            final ReadsReader readsReader = new ReadsReader(new FileInputStream(inputFilename));
            for (final Reads.ReadEntry readEntry : readsReader) {

                byteBuffer = toByteBuffer(sequence, byteBuffer, readEntry);

                final CompressedRead read = new CompressedRead(byteBuffer.clone());


                read.readIndex = readEntry.getReadIndex();
                // we keep track of all the read indices. The ones that go into otherReadIndices
                // will be given a multiplicity of 1.
                otherReadIndices.set(read.readIndex, true);

                if (inspectHashcodes.contains(sd.digest(readEntry.getSequence(), 0, readLength))) {
                    final int count = tallies.getInt(read) + 1;
                    tallies.put(read, count);
                    if (count > 1) {
                        // We have seen this sequence already.
                        // Remove the read from the filter. The first read with the sequence will be searched and
                        // results adjusted according to multiplicity (count)
                        otherReadIndices.set(read.readIndex, false);
                    }
                }

                progress.lightUpdate();

                numReads++;


                //   if (numReads % reportEveryNReads == 1 && numReads != 1) {
                //    printStats(tallies, numReads, false, otherReadIndices);
                //    }

                if (numReads > MAX_PROCESS_READS) {
                    break;
                }
            }
            readsReader.close();
            System.gc();
            printStats(tallies, numReads, true, otherReadIndices);

        }

        progress.stop("second pass");
        System.exit(0);
    }

    private byte[] toByteBuffer(final MutableString sequence, byte[] byteBuffer, final Reads.ReadEntry readEntry) throws IOException {
        ReadsReader.decodeSequence(readEntry, sequence);
        final int i = sequence.length();
        if ((int) (i / 4) + 1 != byteBuffer.length) {
            byteBuffer = new byte[((int) (i / 4)) + 1]; //2 bits per base require four times less space.
        } else {
            Arrays.fill(byteBuffer, (byte) 0);
        }
        final OutputBitStream compressed = new OutputBitStream(byteBuffer);
        int index = 0;
        for (final char c : sequence.array()) {
            if (colorSpace) {

                if (index != 0) {
                    // color space encoded as color space:
                    encodeColorSpace(c, compressed);
                } else {
                    // first base encoded as ATCG:
                    encode(c, compressed);
                }

            } else {
                encode(c, compressed);
            }
            ++index;
            if (index >= i) {
                break;
            }
        }
        compressed.flush();
        return byteBuffer;
    }

    public static void toByteBuffer(final CharSequence sequence,
                                          final byte[] byteBuffer,
                                          final boolean colorSpace,
                                          final int maxReadLength) throws IOException {


        final OutputBitStream compressed = new OutputBitStream(byteBuffer);
        int index = 0;
        for (int i = 0; i < maxReadLength; i++) {
            final char c = sequence.charAt(i);
            if (colorSpace) {

                if (index != 0) {
                    // color space encoded as color space:
                    encodeColorSpace(c, compressed);
                } else {
                    // first base encoded as ATCG:
                    encode(c, compressed);
                }

            } else {
                encode(c, compressed);
            }
            ++index;
            if (index >= maxReadLength) {
                break;
            }
        }
        compressed.flush();

    }

    private void printStats(final Object2IntMap<CompressedRead> tallies, final int numReads, final boolean writeTallies, final BitVector otherReadIndices) {
        final IntList counts = new IntArrayList();
        int sum = 0;
        int num = 0;
        for (final int count : tallies.values()) {
            if (count > 1) {
                counts.add(count);
                sum += count;
                num++;
            }
        }

        Collections.sort(counts);

        if (writeTallies) {
            final ReadSet set = new ReadSet();
            set.smallestStoredMultiplicity(1);
            for (int readIndex = 0; readIndex < otherReadIndices.size(); ++readIndex) {
               if (otherReadIndices.getBoolean(readIndex)) {
                   set.add(readIndex, 1);
               }

            }
            for (final CompressedRead read : tallies.keySet()) {

                final int count = tallies.getInt(read);
                if (count > 1 && otherReadIndices.get(read.readIndex)) {
                    set.add(read.readIndex, count);

                }
            }


            try {
                set.save(outputBasename, "keep");
                System.out.printf("Saved filter with %d elements %n", set.size());
            } catch (IOException e) {
                System.out.println("Error saving read set: " + e);
                System.exit(1);
            }
        }
        System.out.println("Number of reads: " + numReads);
        System.out.printf("Number of redundant reads: %d %n", num);
        System.out.printf("Redunduncy sum: %d %n", sum);
        // we still need to map the first redundant read:
        final int avoidableMappings = sum - counts.size();
        System.out.printf("Could avoid: %d alignments %n", avoidableMappings);
        System.out.printf("Fraction of redundant reads: %3.3g %% %n", (100d * ((double) num) / (double) numReads));
        System.out.printf("Fraction of time saved: %3.3g %% %n", (100d * ((double) avoidableMappings) / (double) numReads));
    }


    public static void main(final String[] args) throws IOException, JSAPException {
        new TallyReadsMode().configure(args).execute();
    }

    public static void encode(final int c,
                              final OutputBitStream compressed) throws IOException {
        switch (c) {
            case 'A':
                compressed.writeBit(1);
                compressed.writeBit(1);

                break;
            case 'C':
                compressed.writeBit(0);
                compressed.writeBit(1);

                break;
            case 'T':
                compressed.writeBit(1);
                compressed.writeBit(0);

                break;
            case 'G':
                compressed.writeBit(0);
                compressed.writeBit(0);

                break;
            default:
                if (c != 'N') {
                    System.out.printf("Unknown base detected: %c %n", (char) c);
                }
                compressed.writeBit(0);
                compressed.writeBit(0);

                break;
        }
    }

    public static void encodeColorSpace(final int c,
                                        final OutputBitStream compressed) throws IOException {
        switch (c) {
            case '0':
                compressed.writeBit(1);
                compressed.writeBit(1);

                break;
            case '1':
                compressed.writeBit(0);
                compressed.writeBit(1);

                break;
            case '2':
                compressed.writeBit(1);
                compressed.writeBit(0);

                break;
            case '3':
                compressed.writeBit(0);
                compressed.writeBit(0);

                break;
            default:
                System.out.printf("Unknown base detected: %c %n", (char) c);
                compressed.writeBit(0);
                compressed.writeBit(0);

                break;
        }
    }
}
