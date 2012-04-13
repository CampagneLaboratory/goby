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

import com.google.protobuf.ByteString;
import com.google.protobuf.GeneratedMessage;
import com.google.protobuf.Message;
import edu.cornell.med.icb.goby.algorithmic.compression.FastArithmeticCoder;
import edu.cornell.med.icb.goby.algorithmic.compression.FastArithmeticDecoder;
import edu.cornell.med.icb.goby.compression.ProtobuffCollectionHandler;
import edu.cornell.med.icb.goby.util.dynoptions.DynamicOptionClient;
import edu.cornell.med.icb.goby.util.dynoptions.RegisterThis;
import it.unimi.dsi.bits.Fast;
import it.unimi.dsi.compression.CodeWordCoder;
import it.unimi.dsi.compression.Decoder;
import it.unimi.dsi.compression.HuffmanCodec;
import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.fastutil.io.FastByteArrayInputStream;
import it.unimi.dsi.fastutil.objects.*;
import it.unimi.dsi.io.InputBitStream;
import it.unimi.dsi.io.OutputBitStream;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.*;
import java.util.Arrays;
import java.util.List;

/**
 * A handler for collections that contain alignment entries.
 *
 * @author Fabien Campagne
 *         Date: 3/3/12
 *         Time: 11:45 AM
 */
public class AlignmentCollectionHandler implements ProtobuffCollectionHandler {
    /**
     * Used to log informational and debug messages.
     */
    private static final Log LOG = LogFactory.getLog(AlignmentCollectionHandler.class);
    @RegisterThis
    public static DynamicOptionClient doc = new DynamicOptionClient(AlignmentCollectionHandler.class,
            "stats-filename:string, the file where to append statistics to:compress-stats.tsv",
            "debug-level:integer, a number between zero and 2. Numbers larger than zero activate debugging. 1 writes stats to stats-filename.:0",
            "basename:string, a basename for the file being converted.:",
            "ignore-read-origin:boolean, When this flag is true do not compress read origin/read groups.:false"

    );
    private String statsFilename;
    private String basename;
    private PrintWriter statsWriter;
    private static final IntArrayList EMPTY_LIST = new IntArrayList();
    private boolean storeReadOrigins = true;
    private final boolean useArithmeticCoding = true;
    private final boolean useHuffmanCoding = false;
    private int streamVersion;


    public static DynamicOptionClient doc() {
        return doc;
    }

    private int previousPosition;
    private int previousTargetIndex;
    private int deltaPosIndex = 0;
    private int qualScoreIndex = 0;
    private int debug = 0;
    /**
     * This variable keeps track of the number of chunks compressed or decompressed.
     */
    private int chunkIndex = 0;
    //Two types of encoding currently supported for query indices:
    private static final int DELTA_ENCODING_SCHEME = 0;
    private static final int MINIMAL_BINARY_ENCODING_SCHEME = 1;
    private static final int MISSING_VALUE = -1;
    private boolean multiplicityFieldsAllMissing = true;
    private long writtenBits;
    private long writtenBases;
    private static final int NO_VALUE = MISSING_VALUE;
    private int varToQualLengthIndex = 0;
    private boolean useTemplateBasedCompression = true;
    private static final int LOG2_8 = Fast.mostSignificantBit(8);

    public AlignmentCollectionHandler() {
        for (int length = 0; length < qualArrays.length; length++) {
            qualArrays[length] = new byte[length];
        }
        debug = doc().getInteger("debug-level");
        storeReadOrigins = !doc().getBoolean("ignore-read-origin");
        statsFilename = doc().getString("stats-filename");
        basename = doc().getString("basename");

        if (debug(1)) {
            try {
                final boolean appending = new File(statsFilename).exists();
                final FileWriter fileWriter = appending ? new FileWriter(statsFilename, true) : new FileWriter(statsFilename);
                statsWriter = new PrintWriter(fileWriter);
                if (!appending) {
                    statsWriter.print("basename\tchunkIndex\tlabel\tnumElements\ttotalBitsWritten\tBitsPerElement\n");
                }

            } catch (FileNotFoundException e) {
                LOG.error("Cannot open stats file", e);
            } catch (IOException e) {
                LOG.error("Cannot open stats file", e);
            }
        }
    }

    @Override
    public int getType() {
        return TYPE_ALIGNMENTS;
    }

    @Override
    public GeneratedMessage parse(final InputStream uncompressedStream) throws IOException {
        return Alignments.AlignmentCollection.parseFrom(uncompressedStream);
    }

    int numChunksProcessed = 0;
    /**
     * The version of the stream that this class reads and writes.
     */
    public static final int VERSION = 6;

    @Override
    public Message compressCollection(final Message collection, final ByteArrayOutputStream compressedBits) throws IOException {
        reset();
        final Alignments.AlignmentCollection alignmentCollection = (Alignments.AlignmentCollection) collection;
        final Alignments.AlignmentCollection.Builder remainingCollection = Alignments.AlignmentCollection.newBuilder();
        final int size = alignmentCollection.getAlignmentEntriesCount();
        int indexInReducedCollection = 0;
        collectStrings(alignmentCollection.getAlignmentEntriesCount(), alignmentCollection);
        for (int index = 0; index < size; index++) {
            final Alignments.AlignmentEntry entry = alignmentCollection.getAlignmentEntries(index);

            final Alignments.AlignmentEntry transformed = transform(index, indexInReducedCollection, entry);
            if (transformed != null) {
                remainingCollection.addAlignmentEntries(transformed);
                indexInReducedCollection++;
                //          System.out.println("not a duplicate");
            } else {

            }
        }
        final OutputBitStream outputBitStream = new OutputBitStream(compressedBits);
        //    System.out.println("queryIndex="+((Alignments.AlignmentCollection) collection).getAlignmentEntries(0).getQueryIndex());
        writeCompressed(outputBitStream);
        outputBitStream.flush();
        writtenBits += outputBitStream.writtenBits();
        if (numChunksProcessed++ % 200 == 0) {
            displayStats();
        }
        ++chunkIndex;
        return remainingCollection.build();
    }

    private void collectStrings(int size, Alignments.AlignmentCollection alignmentCollection) {

        int indexInStrings = 0;
        for (int index = 0; index < size; index++) {
            final Alignments.AlignmentEntry entry = alignmentCollection.getAlignmentEntries(index);
            numSoftClipRightBases.add(entry.hasSoftClippedBasesRight() ? entry.getSoftClippedBasesRight().length() : MISSING_VALUE);
            numSoftClipLeftBases.add(entry.hasSoftClippedBasesLeft() ? entry.getSoftClippedBasesLeft().length() : MISSING_VALUE);
        }

        boolean finished1 = false;
        boolean finished2 = false;
        while (!(finished1 && finished2)) {
            finished1 = true;
            finished2 = true;
            for (int index = 0; index < size; index++) {
                final Alignments.AlignmentEntry entry = alignmentCollection.getAlignmentEntries(index);

                final int numSoftClipLeftBasesInt = numSoftClipLeftBases.getInt(index);
                if (numSoftClipLeftBasesInt != MISSING_VALUE) {
                    if (indexInStrings < numSoftClipLeftBasesInt) {
                        final String softClippedBasesLeft = entry.getSoftClippedBasesLeft();
                        softClipLeftBases.add(softClippedBasesLeft.charAt(indexInStrings));
                        finished1 = false;
                    }
                }
                final int numSoftClipRightBasesInt = numSoftClipRightBases.getInt(index);
                if (numSoftClipRightBasesInt != MISSING_VALUE) {
                    if (indexInStrings < numSoftClipRightBasesInt) {
                        final String softClippedBasesRight = entry.getSoftClippedBasesRight();

                        softClipRightBases.add(softClippedBasesRight.charAt(indexInStrings));
                        finished2 = false;
                    }
                }
            }
            indexInStrings++;
        }

    }

    private void restoreStrings(Alignments.AlignmentCollection.Builder alignmentCollection) {
        final int size = alignmentCollection.getAlignmentEntriesCount();
        int indexInStrings = 0;
        boolean finished1 = false;
        boolean finished2 = false;
        ObjectArrayList<MutableString> softClipsLeft = new ObjectArrayList<MutableString>();
        ObjectArrayList<MutableString> softClipsRight = new ObjectArrayList<MutableString>();
        softClipsLeft.size(size);
        softClipsRight.size(size);
        for (int index = 0; index < size; index++) {

            if (index < numSoftClipLeftBases.size()) {
                final int stringLength = numSoftClipLeftBases.getInt(index);
                if (stringLength != MISSING_VALUE) {

                    final MutableString mutableString = new MutableString();
                    mutableString.setLength(stringLength);
                    softClipsLeft.set(index, mutableString);
                }
            }

            if (index < numSoftClipRightBases.size()) {
                final int stringLength = numSoftClipRightBases.getInt(index);
                if (stringLength != MISSING_VALUE) {

                    final MutableString mutableString = new MutableString();
                    mutableString.setLength(stringLength);
                    softClipsRight.set(index, mutableString);
                }
            }

        }
        int iLeft = 0;
        int iRight = 0;
        while (!(finished1 && finished2)) {
            finished1 = true;
            finished2 = true;
            for (int index = 0; index < size; index++) {
                final int numSoftClipLeftBasesInt = numSoftClipLeftBases.getInt(index);
                if (numSoftClipLeftBasesInt != MISSING_VALUE && indexInStrings < numSoftClipLeftBasesInt) {
                    final MutableString mutableString = softClipsLeft.get(index);
                    if (mutableString != null) {
                        final int anInt = softClipLeftBases.getInt(iLeft++);

                        mutableString.setCharAt(indexInStrings, (char) anInt);
                        finished1 = false;
                    }
                }
                final int numSoftClipRightBasesInt = numSoftClipRightBases.getInt(index);

                if (numSoftClipRightBasesInt != MISSING_VALUE && indexInStrings < numSoftClipRightBasesInt) {
                    final MutableString mutableString = softClipsRight.get(index);
                    if (mutableString != null) {
                        final int anInt = softClipRightBases.getInt(iRight++);

                        mutableString.setCharAt(indexInStrings, (char) anInt);
                        finished2 = false;
                    }
                }
            }
            indexInStrings++;
        }

        for (int index = 0; index < size; index++) {
            Alignments.AlignmentEntry.Builder builder = alignmentCollection.getAlignmentEntriesBuilder(index);
            {
                final MutableString mutableString = softClipsLeft.get(index);
                if (mutableString != null) {
                    builder.setSoftClippedBasesLeft(mutableString.toString());
                }
            }
            {
                final MutableString mutableString = softClipsRight.get(index);
                if (mutableString != null) {
                    builder.setSoftClippedBasesRight(mutableString.toString());
                }
            }

        }
    }


    @Override
    public Message decompressCollection(Message reducedCollection, byte[] compressedBytes) throws IOException {
        reset();
        //TODO optimize away the copy:
        byte[] moreRoom = new byte[compressedBytes.length + 100];
        System.arraycopy(compressedBytes, 0, moreRoom, 0, compressedBytes.length);

        final Alignments.AlignmentCollection alignmentCollection = (Alignments.AlignmentCollection) reducedCollection;
        final Alignments.AlignmentCollection.Builder result = Alignments.AlignmentCollection.newBuilder();
        final InputBitStream bitInput = new InputBitStream(new FastByteArrayInputStream(moreRoom));
        final int numEntriesInChunk = alignmentCollection.getAlignmentEntriesCount();

        final int streamVersion = decompressBits(bitInput, numEntriesInChunk);
        int originalIndex = 0;
        for (int templateIndex = 0; templateIndex < numEntriesInChunk; templateIndex++) {
            final int templatePositionIndex = varPositionIndex;
            final int templateVarFromToIndex = varFromToIndex;
            final int templateVarHasToQualsIndex = varToQualLengthIndex;
            while (multiplicities.get(templateIndex) >= 1) {
                result.addAlignmentEntries(
                        andBack(templateIndex, originalIndex, alignmentCollection.getAlignmentEntries(templateIndex), streamVersion));
                if (multiplicities.get(templateIndex) >= 1) {
                    // go back to the indices for the template:
                    varPositionIndex = templatePositionIndex;
                    varFromToIndex = templateVarFromToIndex;
                }
                originalIndex++;
            }
        }
        restoreStrings(result);
        ++chunkIndex;
        return result.build();
    }


    @Override
    public void setUseTemplateCompression(final boolean useTemplateCompression) {
        useTemplateBasedCompression = useTemplateCompression;
    }


    public void displayStats() {
        if (debug(1)) {
            for (String label : typeToNumEntries.keySet()) {
                int n = typeToNumEntries.getInt(label);
                long written = typeToWrittenBits.getLong(label);
                double average = (double) written / (double) n;
                LOG.info
                        (String.format("encoded %d %s in %d bits, average %g bits /element. ", n, label,
                                written, average));
                statsWriter.write(String.format("%s\t%d\t%s\t%d\t%d\t%g%n", basename, chunkIndex, label, n, written, divide(written, n)));
                statsWriter.flush();
            }
            LOG.info(String.format("entries aggregated with multiplicity= %d", countAggregatedWithMultiplicity));
            LOG.info(String.format("Overall: bits per aligned bases= %g", divide(writtenBits, writtenBases)));
        }
    }

    private double divide(long a, long b) {
        return ((double) a / (double) b);
    }

    protected final boolean debug(int level) {
        return debug >= level;
    }


    private void writeInts(final String label, final IntList list, final OutputBitStream out) throws IOException {
        final long writtenStart = out.writtenBits();

        for (final int value : list) {
            out.writeInt(value, 32);
        }
        final long writtenStop = out.writtenBits();
        final long written = writtenStop - writtenStart;
        recordStats(label, list, written);
    }

    private void writeQueryIndices(final String label, final IntList list, final OutputBitStream out) throws IOException {
        boolean success = tryWriteDeltas(label, list, out);
        if (success) {
            return;
        } else {
            out.writeDelta(MINIMAL_BINARY_ENCODING_SCHEME);
        }

        int min = Integer.MAX_VALUE;
        int max = Integer.MIN_VALUE;
        final int size = list.size();
        //   System.out.printf("encoding, chunk=%d delta-positions.size=%d%n", chunkIndex, deltaPositions.size());
        for (final int value : list) {

            min = Math.min(value, min);
            max = Math.max(value, max);

        }

        out.writeNibble(size);
        if (size == 0) {
            return;
        }
        out.writeNibble(min);
        out.writeNibble(max);
        // add one to each value, since we cannot write zeroes in minimal binary.

        final long writtenStart = out.writtenBits();
        final int b = max - min + 1;
        final int log2b = Fast.mostSignificantBit(b);
        for (final int value : list) {

            out.writeMinimalBinary(value - min, b, log2b);
            // out.writeLongMinimalBinary(value-min, max-min+1);
        }
        //    out.flush();
        if (debug(1)) {
            //   out.flush();
            final long writtenStop = out.writtenBits();
            final long written = writtenStop - writtenStart;
            recordStats(label, list, written);
        }
    }

    /**
     * Write a list with Rice/Golomb coding
     *
     * @param label
     * @param list
     * @param out
     * @throws IOException
     */
    public void writeRiceCoding(String label, IntList list, OutputBitStream out) throws IOException {
        final long writtenStart = out.writtenBits();
        out.writeNibble(list.size());
        for (final int value : list) {

            out.writeGolomb(value, 8, LOG2_8);
        }
        if (debug(1)) {
            System.err.flush();
            final long writtenStop = out.writtenBits();
            final long written = writtenStop - writtenStart;
            recordStats(label, list, written);
        }
    }

    /**
     * Try to write query indices as delta. If the number of unique deltas is larger than 10% of list size, do not
     * write anything and return false. Otherwise, write as delta and return true.
     *
     * @param label
     * @param list
     * @param out
     * @return
     * @throws java.io.IOException
     */
    private boolean tryWriteDeltas(String label, IntList list, OutputBitStream out) throws IOException {
        if (list.size() == 0) {
            return false;
        }
        final IntArrayList deltas = new IntArrayList();
        int first = list.getInt(0);
        // write the first value as is:
        int previous = first;
        int index = 0;
        for (int value : list) {
            if (index > 0) {
                deltas.add(Fast.int2nat(value - previous));
                previous = value;
            }
            ++index;
        }

        final IntSet tokens = getTokens(deltas);
        //   System.out.printf("tokenSize=%d listSize=%d%n", tokens.size(), list.size());
        if (divide(tokens.size(), list.size()) > 0.2f) {
            return false;
        } else {
            //     System.out.println("Using delta encoding scheme");
            out.writeDelta(DELTA_ENCODING_SCHEME);
            out.writeNibble(first);
            writeArithmetic(label, deltas, out);
            return true;
        }

    }

    private float divide(int a, int b) {
        return ((float) a) / ((float) b);
    }

    private void decodeQueryIndices(final String label, final int numEntriesInChunk, final InputBitStream bitInput, final IntList list) throws IOException {
        switch (bitInput.readDelta()) {
            case MINIMAL_BINARY_ENCODING_SCHEME:
                readMinimalUnary(label, numEntriesInChunk, bitInput, list);
                break;
            case DELTA_ENCODING_SCHEME:
                readAsDeltas(label, numEntriesInChunk, bitInput, list);
                break;
        }
    }

    private void readAsDeltas(String label, int numEntriesInChunk, InputBitStream bitInput, IntList list) throws IOException {
        IntArrayList deltas = new IntArrayList();
        int previous = bitInput.readNibble();
        decodeArithmetic(label, numEntriesInChunk - 1, bitInput, deltas);
        list.add(previous);
        for (int delta : deltas) {
            final int newValue = Fast.nat2int(delta) + previous;
            list.add(newValue);
            previous = newValue;
        }

    }

    private void readMinimalUnary(final String label, final int numEntriesInChunk, final InputBitStream bitInput, final IntList list) throws IOException {
        final int size = bitInput.readNibble();
        if (size == 0) {
            return;
        }
        final int min = bitInput.readNibble();
        final int max = bitInput.readNibble();
        final int b = max - min + 1;
        final int log2b = Fast.mostSignificantBit(b);
        for (int i = 0; i < size; i++) {
            final int reducedReadIndex = bitInput.readMinimalBinary(max - min + 1, log2b);
            list.add(reducedReadIndex + min);
        }
        //  bitInput.flush();
    }

    private void writeNibble(String label, IntList list, OutputBitStream out) throws IOException {
        long writtenStart = out.writtenBits();
        for (int value : list) {
            out.writeNibble(value);
        }
        long writtenStop = out.writtenBits();
        long written = writtenStop - writtenStart;
        recordStats(label, list, written);
    }

    Object2IntMap<String> typeToNumEntries = new Object2IntAVLTreeMap<String>();
    Object2LongMap<String> typeToWrittenBits = new Object2LongAVLTreeMap<String>();


    protected final void decodeArithmetic(final String label, final int numEntriesInChunk, final InputBitStream bitInput, final IntList list) throws IOException {

        if (debug(2)) {
            System.err.flush();
            System.err.println("\nreading " + label + " with available=" + bitInput.available());
            System.err.flush();
        }
        boolean useRunLength = false;
        if (streamVersion >= 5) {
            // version 5 and up choose to use run length encoding for each field:
            useRunLength = bitInput.readBit() == 1;
        }
        if (useRunLength) {
            final IntArrayList encodedLengths = new IntArrayList();
            final IntArrayList encodedValues = new IntArrayList();
            decodeArithmeticInternal(bitInput, encodedLengths);
            decodeArithmeticInternal(bitInput, encodedValues);

            decodeRunLengths(encodedLengths, encodedValues, list);
        } else {
            decodeArithmeticInternal(bitInput, list);
        }


    }

    private void decodeArithmeticInternal(InputBitStream bitInput, IntList list) throws IOException {
        final int size = bitInput.readNibble();
        if (size == 0) {
            return;
        }
        final boolean hasNegatives = bitInput.readBit() == 1;
        final int numTokens = bitInput.readNibble();

        final int[] distinctvalue = new int[numTokens];
        for (int i = 0; i < numTokens; i++) {
            // -1 makes 0 symbol -1 (missing value) again
            final int token = bitInput.readNibble();
            final int anInt = hasNegatives ? Fast.nat2int(token) : token - 1;
            distinctvalue[i] = anInt;
        }
        if (hasNegatives) {
            // we must sort the symbol values again since the bijection has permuted them
            Arrays.sort(distinctvalue);
        }
        decode(bitInput, list, size, numTokens, distinctvalue);
    }

    protected final void writeArithmetic(final String label, final IntList list, OutputBitStream out) throws IOException {
        if (debug(2)) {
            System.err.flush();
            System.err.println("\nwriting " + label);
            System.err.flush();
        }
        final long writtenStart = out.writtenBits();
        IntArrayList encodedLengths = new IntArrayList();
        IntArrayList encodedValues = new IntArrayList();
        encodeRunLengths(list, encodedLengths, encodedValues);
        if (runLengthEncoding(label, list, encodedLengths, encodedValues)) {
            out.writeBit(1);
            encodeArithmeticInternal(label + "Lengths", encodedLengths, out);
            encodeArithmeticInternal(label + "Values", encodedValues, out);
        } else {
            out.writeBit(0);
            encodeArithmeticInternal(label, list, out);
        }
        if (debug(1)) {
            System.err.flush();
            final long writtenStop = out.writtenBits();
            final long written = writtenStop - writtenStart;
            recordStats(label, list, written);
        }
    }

    private boolean runLengthEncoding(String label, IntList list, IntArrayList encodedLengths, IntArrayList encodedValues) {

        final boolean result = encodedLengths.size() > 10 && (encodedLengths.size() + encodedValues.size()) < list.size() ;
        if (result) {
      //          System.out.println("Using run-length encoding for "+label);
        }
        return result;
    }

    private boolean encodeArithmeticInternal(String label, IntList list, OutputBitStream out) throws IOException {
        out.writeNibble(list.size());        // LIST.size
        if (list.isEmpty()) {
            // no list to write.
            return true;
        }
        final IntSet distinctSymbols = getTokens(list);

        final int[] symbolValues = distinctSymbols.toIntArray();
        final boolean hasNegativeValues = hasNegatives(symbolValues);
        out.writeBit(hasNegativeValues);
        out.writeNibble(distinctSymbols.size());    //LIST.distinct-values.size()
        for (final int token : distinctSymbols) {
            // +1 makes -1 (missing value) symbol 0 so it can be written Nibble:

            final int anInt = token;
            final int nat = hasNegativeValues ? Fast.int2nat(anInt) : anInt + 1;

            out.writeNibble(nat);
        }

        encode(label, list, out, distinctSymbols, symbolValues);
        return false;
    }

    private void encode(String label, final IntList list, final OutputBitStream out, final IntSet distinctSymbols, final int[] symbolValues) throws IOException {
        if (useArithmeticCoding) {
            final FastArithmeticCoder coder = new FastArithmeticCoder(distinctSymbols.size());
            for (final int dp : list) {
                final int symbolCode = Arrays.binarySearch(symbolValues, dp);
                assert symbolCode >= 0 : "symbol code must exist.";
                coder.encode(symbolCode, out);
            }
            coder.flush(out);
        } else if (useHuffmanCoding) {
            final int[] frequencies = frequencies(list, symbolValues);
            final HuffmanCodec codec = new HuffmanCodec(frequencies);
            final CodeWordCoder coder = codec.coder();
            for (int freq : frequencies) {
                out.writeNibble(freq);
            }
            for (final int dp : list) {
                final int symbolCode = Arrays.binarySearch(symbolValues, dp);
                assert symbolCode >= 0 : "symbol code must exist.";
                coder.encode(symbolCode, out);
            }
            coder.flush(out);
        }
    }

    private void decode(final InputBitStream bitInput, final IntList list, final int size, final int numTokens, int[] distinctvalue) throws IOException {
        if (useArithmeticCoding) {
            final FastArithmeticDecoder decoder = new FastArithmeticDecoder(numTokens);
            for (int i = 0; i < size; i++) {
                final int tokenValue = distinctvalue[decoder.decode(bitInput)];
                list.add(tokenValue);
            }
            decoder.reposition(bitInput);
        } else if (useHuffmanCoding) {
            final int[] frequencies = new int[distinctvalue.length];
            for (int i = 0; i < frequencies.length; i++) {
                frequencies[i] = bitInput.readNibble();
            }
            final HuffmanCodec codec = new HuffmanCodec(frequencies);
            final Decoder decoder = codec.decoder();
            for (int i = 0; i < size; i++) {
                final int tokenValue = distinctvalue[decoder.decode(bitInput)];
                list.add(tokenValue);
            }
            // decoder.reposition(bitInput);
        }
    }

    // return the frequencies of symbols in the list
    private int[] frequencies(IntList list, int[] symbolValues) {
        final int[] freqs = new int[symbolValues.length];
        for (final int value : list) {
            final int index = Arrays.binarySearch(symbolValues, value);
            freqs[index] += 1;
        }
        return freqs;
    }

    private boolean hasNegatives(final int[] symbolValues) {
        for (final int val : symbolValues) {
            if (val < -1) {
                return true;
            }
        }
        return false;
    }

    private void recordStats(String label, IntList list, long written) {
        if (debug(1)) {
            double average = ((double) written) / list.size();
            typeToNumEntries.put(label, list.size() + typeToNumEntries.getInt(label));
            typeToWrittenBits.put(label, written + typeToWrittenBits.getLong(label));
        }
    }

    final IntSortedSet tokenSet = new IntAVLTreeSet();

    private final IntSortedSet getTokens(final IntList list) {
        tokenSet.clear();
        for (final int value : list) {
            tokenSet.add(value);
        }
        return tokenSet;
    }


    private IntList deltaPositions = new IntArrayList();
    private IntList deltaTargetIndices = new IntArrayList();
    private IntList queryLengths = new IntArrayList();
    private IntList mappingQualities = new IntArrayList();
    private IntList matchingReverseStrand = new IntArrayList();
    private IntList multiplicity = new IntArrayList();
    private IntList numberOfIndels = new IntArrayList();
    private IntList numberOfMismatches = new IntArrayList();
    private IntList queryAlignedLengths = new IntArrayList();
    private IntList targetAlignedLengths = new IntArrayList();
    private IntList queryIndices = new IntArrayList();
    private IntList queryPositions = new IntArrayList();
    private IntList fragmentIndices = new IntArrayList();
    private IntList variationCount = new IntArrayList();
    private IntList insertSizes = new IntArrayList();
    private IntList fromLengths = new IntArrayList();
    private IntList toLengths = new IntArrayList();
    private IntList varPositions = new IntArrayList();
    private IntList varReadIndex = new IntArrayList();
    private IntList varFromTo = new IntArrayList();
    private IntList varQuals = new IntArrayList();
    private IntList varToQualLength = new IntArrayList();
    private IntList allReadQualityScores = new IntArrayList();
    private IntList sampleIndices = new IntArrayList();
    private IntList readOriginIndices = new IntArrayList();
    private IntList pairFlags = new IntArrayList();
    private IntList scores = new IntArrayList();
    private IntArrayList numReadQualityScores = new IntArrayList();

    private IntArrayList multiplicities = new IntArrayList();

    private IntList numSoftClipLeftBases = new IntArrayList();
    private IntList numSoftClipRightBases = new IntArrayList();
    private IntList softClipLeftBases = new IntArrayList();
    private IntList softClipRightBases = new IntArrayList();


    private int decompressBits(InputBitStream bitInput, final int numEntriesInChunk) throws IOException {
        streamVersion = bitInput.readDelta();

        assert streamVersion <= VERSION : "FATAL: The stream version cannot have been written with a more recent version of Goby (The hybrid chunk codec cannot not support forward compatibility of the compressed stream).";
        multiplicityFieldsAllMissing = bitInput.readBit() == 1;

        decodeArithmetic("deltaPositions", numEntriesInChunk, bitInput, deltaPositions);
        decodeArithmetic("deltaTargetIndices", numEntriesInChunk, bitInput, deltaTargetIndices);
        decodeArithmetic("queryLengths", numEntriesInChunk, bitInput, queryLengths);
        decodeArithmetic("mappingQualities", numEntriesInChunk, bitInput, mappingQualities);
        decodeArithmetic("matchingReverseStrand", numEntriesInChunk, bitInput, matchingReverseStrand);
        decodeArithmetic("numberOfIndels", numEntriesInChunk, bitInput, numberOfIndels);
        decodeArithmetic("numberOfMismatches", numEntriesInChunk, bitInput, numberOfMismatches);
        decodeArithmetic("insertSize", numEntriesInChunk, bitInput, insertSizes);
        decodeArithmetic("queryAlignedLength", numEntriesInChunk, bitInput, queryAlignedLengths);
        decodeArithmetic("targetAlignedLength", numEntriesInChunk, bitInput, targetAlignedLengths);
        decodeArithmetic("queryPositions", numEntriesInChunk, bitInput, queryPositions);
        decodeArithmetic("fragmentIndex", numEntriesInChunk, bitInput, fragmentIndices);
        decodeArithmetic("variationCount", numEntriesInChunk, bitInput, variationCount);
        decodeArithmetic("varPositions", numEntriesInChunk, bitInput, varPositions);
        decodeArithmetic("fromLengths", numEntriesInChunk, bitInput, fromLengths);
        decodeArithmetic("toLengths", numEntriesInChunk, bitInput, toLengths);
        decodeArithmetic("varReadIndex", numEntriesInChunk, bitInput, varReadIndex);
        decodeArithmetic("varFromTo", numEntriesInChunk, bitInput, varFromTo);
        decodeArithmetic("varQuals", numEntriesInChunk, bitInput, varQuals);
        decodeArithmetic("varHasToQuals", numEntriesInChunk, bitInput, varToQualLength);
        decodeArithmetic("multiplicities", numEntriesInChunk, bitInput, multiplicities);
        pairLinks.read(numEntriesInChunk, bitInput);
        forwardSpliceLinks.read(numEntriesInChunk, bitInput);
        backwardSpliceLinks.read(numEntriesInChunk, bitInput);

        decodeQueryIndices("queryIndices", numEntriesInChunk, bitInput, queryIndices);

        if (streamVersion >= 2) {

            decodeArithmetic("numReadQualityScores", numEntriesInChunk, bitInput, numReadQualityScores);
            decodeArithmetic("allReadQualityScores", numEntriesInChunk, bitInput, allReadQualityScores);
        }
        if (streamVersion >= 3) {

            decodeArithmetic("sampleIndices", numEntriesInChunk, bitInput, sampleIndices);
            decodeArithmetic("readOriginIndices", numEntriesInChunk, bitInput, readOriginIndices);
        }
        if (streamVersion >= 4) {

            decodeArithmetic("pairFlags", numEntriesInChunk, bitInput, pairFlags);
            decodeArithmetic("scores", numEntriesInChunk, bitInput, scores);
        }
        if (streamVersion >= 6) {

            decodeArithmetic("softClipLeftBasesNum", numEntriesInChunk, bitInput, numSoftClipLeftBases);
            decodeArithmetic("softClipRightBasesNum", numEntriesInChunk, bitInput, numSoftClipRightBases);
            decodeArithmetic("softClipLeftBases", numEntriesInChunk, bitInput, softClipLeftBases);
            decodeArithmetic("softClipRightBases", numEntriesInChunk, bitInput, softClipRightBases);
        }
        return streamVersion;
    }


    private void decodeRunLengths(final IntArrayList encodedLengths, final IntArrayList encodedValues, final IntList list) {

        final int size = encodedLengths.size();
        int index = 0;
        int valueIndex = 0;
        for (index = 0; index < size; index++) {
            final int runLength = encodedLengths.get(index);
            for (int j = 0; j < runLength; j++) {

                list.add(encodedValues.get(valueIndex));
            }
            valueIndex += 1;
        }
    }

    private void encodeRunLengths(IntList list, IntArrayList encodedLengths, IntArrayList encodedValues) {
        if (list.size() == 0) {
            return;
        }
        int previous = list.get(0);
        int runLength = 1;

        final int size = list.size();
        for (int index = 1; index < size; index++) {
            int value = list.get(index);
            if (value == previous) {
                runLength += 1;
            } else {
                encodedLengths.add(runLength);
                encodedValues.add(previous);
                runLength = 1;
            }
            previous = value;
        }

        encodedLengths.add(runLength);
        encodedValues.add(previous);

    }

    private void writeCompressed(final OutputBitStream out) throws IOException {
        //   out.writeNibble(0);
        out.writeDelta(VERSION);
        out.writeBit(multiplicityFieldsAllMissing);

        writeArithmetic("positions", deltaPositions, out);
        writeArithmetic("targets", deltaTargetIndices, out);
        writeArithmetic("queryLengths", queryLengths, out);
        writeArithmetic("mappingQualities", mappingQualities, out);
        writeArithmetic("matchingReverseStrand", matchingReverseStrand, out);
        writeArithmetic("numberOfIndels", numberOfIndels, out);
        writeArithmetic("numberOfMismatches", numberOfMismatches, out);
        writeArithmetic("insertSize", insertSizes, out);
        writeArithmetic("queryAlignedLength", queryAlignedLengths, out);
        writeArithmetic("targetAlignedLength", targetAlignedLengths, out);
        writeArithmetic("queryPositions", queryPositions, out);
        writeArithmetic("fragmentIndex", fragmentIndices, out);
        writeArithmetic("variationCount", variationCount, out);
        writeArithmetic("varPositions", varPositions, out);
        writeArithmetic("fromLengths", fromLengths, out);
        writeArithmetic("toLengths", toLengths, out);
        writeArithmetic("varReadIndex", varReadIndex, out);
        writeArithmetic("varFromTo", varFromTo, out);
        writeArithmetic("varQuals", varQuals, out);
        writeArithmetic("varHasToQuals", varToQualLength, out);
        writeArithmetic("multiplicities", multiplicities, out);
        pairLinks.write(out);
        forwardSpliceLinks.write(out);
        backwardSpliceLinks.write(out);

        writeQueryIndices("queryIndices", queryIndices, out);

        writeArithmetic("numReadQualityScores", numReadQualityScores, out);
        writeArithmetic("allReadQualityScoresLengths", allReadQualityScores, out);

        writeArithmetic("sampleIndices", sampleIndices, out);
        writeArithmetic("readOriginIndices", readOriginIndices, out);
        writeArithmetic("pairFlags", pairFlags, out);
        writeArithmetic("scores", scores, out);

        writeArithmetic("softClipLeftBasesNum", numSoftClipLeftBases, out);
        writeArithmetic("softClipRightBasesNum", numSoftClipRightBases, out);
        writeArithmetic("softClipLeftBases", softClipLeftBases, out);
        writeArithmetic("softClipRightBases", softClipRightBases, out);
    }

    private void reset() {
        multiplicityFieldsAllMissing = true;

        previousPosition = -1;
        previousTargetIndex = -1;
        deltaPositions.clear();
        deltaTargetIndices.clear();
        queryLengths.clear();
        mappingQualities.clear();
        matchingReverseStrand.clear();
        multiplicity.clear();
        numberOfIndels.clear();
        queryAlignedLengths.clear();
        targetAlignedLengths.clear();
        numberOfMismatches.clear();
        insertSizes.clear();
        queryIndices.clear();
        queryPositions.clear();
        fragmentIndices.clear();
        queryIndices.clear();
        variationCount.clear();
        varPositions.clear();
        fromLengths.clear();
        toLengths.clear();
        varReadIndex.clear();
        varFromTo.clear();
        varQuals.clear();
        varQualIndex = 0;
        varPositionIndex = 0;
        varFromToIndex = 0;
        varToQualLength.clear();
        varToQualLengthIndex = 0;
        multiplicities.clear();
        countAggregatedWithMultiplicity = 0;
        previousPartial = null;
        deltaPosIndex = 0;
        pairLinks.reset();
        forwardSpliceLinks.reset();
        backwardSpliceLinks.reset();
        qualScoreIndex = 0;
        numReadQualityScores.clear();
        allReadQualityScores.clear();
        sampleIndices.clear();
        readOriginIndices.clear();
        pairFlags.clear();
        scores.clear();

        numSoftClipLeftBases.clear();
        numSoftClipRightBases.clear();
        softClipLeftBases.clear();
        softClipRightBases.clear();
    }

    private final LinkInfo pairLinks = new LinkInfo(this, "pairs");
    private final LinkInfo forwardSpliceLinks = new LinkInfo(this, "forward-splice");
    private final LinkInfo backwardSpliceLinks = new LinkInfo(this, "backward-splice");

    /**
     * An empty sequence variation.
     */
    private final Alignments.SequenceVariation EMPTY_SEQ_VAR = Alignments.SequenceVariation.newBuilder().build();
    private Alignments.AlignmentEntry previousPartial;
    private int countAggregatedWithMultiplicity;

    private Alignments.AlignmentEntry transform(final int index, int indexInReducedCollection, final Alignments.AlignmentEntry source) {
        final Alignments.AlignmentEntry.Builder result = Alignments.AlignmentEntry.newBuilder(source);
        final int position = source.getPosition();
        final int targetIndex = source.getTargetIndex();

        // clear the strings we collected earlier:
        result.clearSoftClippedBasesLeft();
        result.clearSoftClippedBasesRight();

        if (index > 0 && targetIndex == previousTargetIndex) {
            result.clearPosition();
            result.clearTargetIndex();

            deltaPositions.add(position - previousPosition);
            deltaTargetIndices.add(targetIndex - previousTargetIndex);
        }

        final int queryIndex = source.getQueryIndex();

        queryIndices.add(queryIndex);
        previousPosition = position;
        previousTargetIndex = targetIndex;

        if (debug(1) && source.hasQueryLength()) {
            writtenBases += source.getQueryAlignedLength();
        }
        result.clearQueryIndex();

        recordVariationQualitiesAndClear(source, result, result.getSequenceVariationsList());

        final boolean entryMatchingReverseStrand = source.getMatchingReverseStrand();
        Alignments.RelatedAlignmentEntry link = pairLinks.code(source.hasPairAlignmentLink(), entryMatchingReverseStrand,
                source.getPairAlignmentLink());
        if (link == null) {
            result.clearPairAlignmentLink();
        } else {
            result.setPairAlignmentLink(link);
        }

        link = forwardSpliceLinks.code(source.hasSplicedForwardAlignmentLink(), entryMatchingReverseStrand, source.getSplicedForwardAlignmentLink());
        if (link == null) {
            result.clearSplicedForwardAlignmentLink();
        } else {
            result.setSplicedForwardAlignmentLink(link);
        }

        link = backwardSpliceLinks.code(source.hasSplicedBackwardAlignmentLink(), entryMatchingReverseStrand, source.getSplicedBackwardAlignmentLink());
        if (link == null) {
            result.clearSplicedBackwardAlignmentLink();
        } else {
            result.setSplicedBackwardAlignmentLink(link);
        }
        final Alignments.AlignmentEntry partial = result.clone().build();
        if (source.hasReadQualityScores()) {

            final ByteString quals = source.getReadQualityScores();
            final int size = quals.size();
            numReadQualityScores.add(size);
            for (int i = 0; i < size; i++) {
                allReadQualityScores.add(quals.byteAt(i));
            }
        } else {
            numReadQualityScores.add(0);
        }
        if (previousPartial != null && indexInReducedCollection >= 1 && fastEquals(previousPartial, partial)) {
            //   System.out.println("same");
            //  print(partial);
            int m = multiplicities.get(indexInReducedCollection - 1);
            multiplicities.set(indexInReducedCollection - 1, m + 1);
            // do not add this one, we just increased the multiplicity of the previous one.
            countAggregatedWithMultiplicity++;
            //      System.out.printf("Returning for template match to previous, current queryIndex=%d%n",queryIndex);
            return null;
        } else {
            previousPartial = partial;
            multiplicityFieldsAllMissing &= !source.hasMultiplicity();
            multiplicities.add(Math.max(1, source.getMultiplicity()));
        }
        //System.out.printf("encoding query-index=%d varPositionIndex=%d %n",queryIndex, varPositionIndex);

        queryLengths.add(source.hasQueryLength() ? source.getQueryLength() : MISSING_VALUE);
        mappingQualities.add(source.hasMappingQuality() ? source.getMappingQuality() : MISSING_VALUE);
        matchingReverseStrand.add(source.hasMatchingReverseStrand() ? source.getMatchingReverseStrand() ? 1 : 0 : MISSING_VALUE);
        numberOfIndels.add(source.hasNumberOfIndels() ? source.getNumberOfIndels() : MISSING_VALUE);
        numberOfMismatches.add(source.hasNumberOfMismatches() ? source.getNumberOfMismatches() : MISSING_VALUE);
        if (source.hasInsertSize()) {
            final int readPos = source.getPosition();
            final int matePos = source.getPairAlignmentLink().getPosition();
            final int length = source.getTargetAlignedLength();
            final int pos1 = source.getMatchingReverseStrand() ? length + readPos : readPos + 1;
            final int pos2 = EntryFlagHelper.isMateReverseStrand(source) ? length + matePos : matePos + 1;
            final int insertSize = source.getInsertSize();
            int insertSizeDiff = pos2 - pos1 - insertSize;
            // reverse:  insertSize= (pos2-pos1) - isd
            if (insertSize != 0) {
                /* System.out.printf("insertSize %d length= %d %c %c readPos %d matePos %d  pos1 %d pos2 %d  pos2-pos1 %d matePos-readPos  %d  insertSizeDiff %d %n",
                       source.getInsertSize(), length,
                       source.getMatchingReverseStrand() ? '+' : '-',
                       EntryFlagHelper.isMateReverseStrand(source) ? '+' : '-',
                       readPos, matePos, pos1, pos2, pos2 - pos1, matePos - readPos, insertSizeDiff);
                */
            }

            if (insertSize == 0) {
                insertSizeDiff = MISSING_VALUE;

            }
            insertSizes.add(source.hasInsertSize() ? insertSizeDiff : MISSING_VALUE);
        } else {
            insertSizes.add(MISSING_VALUE);
        }

        queryAlignedLengths.add(source.hasQueryAlignedLength() ? source.getQueryAlignedLength() : MISSING_VALUE);
        targetAlignedLengths.add(source.hasTargetAlignedLength() ? source.getTargetAlignedLength() : MISSING_VALUE);
        fragmentIndices.add(source.hasFragmentIndex() ? source.getFragmentIndex() : MISSING_VALUE);
        variationCount.add(source.getSequenceVariationsCount());
        queryPositions.add(source.hasQueryPosition() ? source.getQueryPosition() : MISSING_VALUE);
        sampleIndices.add(source.hasSampleIndex() ? source.getSampleIndex() : MISSING_VALUE);
        readOriginIndices.add(source.hasReadOriginIndex() && storeReadOrigins ? source.getReadOriginIndex() : MISSING_VALUE);
        pairFlags.add(source.hasPairFlags() ? source.getPairFlags() : MISSING_VALUE);
        scores.add(source.hasScore() ? Float.floatToIntBits(source.getScore()) : MISSING_VALUE);


        result.clearQueryLength();
        result.clearMappingQuality();
        result.clearMatchingReverseStrand();
        result.clearMultiplicity();
        result.clearNumberOfIndels();
        result.clearNumberOfMismatches();
        result.clearInsertSize();
        result.clearQueryAlignedLength();
        result.clearTargetAlignedLength();

        result.clearQueryPosition();
        result.clearFragmentIndex();
        result.clearReadQualityScores();
        result.clearSampleIndex();
        result.clearReadOriginIndex();
        result.clearPairFlags();
        result.clearScore();
        boolean canFullyRemoveThisOne = true;
        boolean canFullyRemoveCollection = true;
        int seqVarIndex = 0;

        for (final Alignments.SequenceVariation seqVar : result.getSequenceVariationsList()) {

            encodeVar(source.getMatchingReverseStrand(), source.getQueryLength(), seqVar);
            Alignments.SequenceVariation.Builder varBuilder = Alignments.SequenceVariation.newBuilder(seqVar);
            varBuilder.clearPosition();
            varBuilder.clearFrom();
            varBuilder.clearTo();
            varBuilder.clearToQuality();
            varBuilder.clearReadIndex();
            if (!fastEquals(EMPTY_SEQ_VAR, varBuilder.build())) {
                canFullyRemoveThisOne = false;
                canFullyRemoveCollection = false;
            }
            if (canFullyRemoveThisOne) {
                result.removeSequenceVariations(seqVarIndex);
                seqVarIndex--;
            }
            seqVarIndex++;
        }
        if (canFullyRemoveCollection) {
            result.clearSequenceVariations();
        }


        final Alignments.AlignmentEntry alignmentEntry = result.build();
        //   System.out.println(alignmentEntry);
        return alignmentEntry;
    }

    private boolean fastEquals(final Object o1, final Object o2) {
        return useTemplateBasedCompression && o1.equals(o2);

    }


    private void recordVariationQualitiesAndClear(Alignments.AlignmentEntry source, Alignments.AlignmentEntry.Builder result, List<Alignments.SequenceVariation> sequenceVariationsList) {
        if (source.hasReadQualityScores()) {
            for (final Alignments.SequenceVariation seqVar : sequenceVariationsList) {
                varToQualLength.add(0);
            }
        } else

        {
            int index = 0;
            for (final Alignments.SequenceVariation seqVar : sequenceVariationsList) {
                final String from = seqVar.getFrom();

                final ByteString toQualities = seqVar.getToQuality();
                final boolean hasToQuals = seqVar.hasToQuality();
                final int toQualSize = hasToQuals ? toQualities.size() : 0;
                varToQualLength.add(toQualSize);

                for (int i = 0; i < toQualSize; i++) {
                    varQuals.add(toQualities.byteAt(i));
                }
                Alignments.SequenceVariation.Builder varBuilder = Alignments.SequenceVariation.newBuilder(seqVar);
                varBuilder.clearToQuality();
                result.setSequenceVariations(index, varBuilder.buildPartial());
                index++;
            }
        }

    }

    private void print(Alignments.AlignmentEntry result) {

        System.out.println(result);
    }

    private void encodeVar(boolean entryOnReverseStrand, int queryLenth, final Alignments.SequenceVariation seqVar) {

        final String from = seqVar.getFrom();
        final String to = seqVar.getTo();
        final int fromLength = from.length();
        final int toLength = to.length();
        final int position = seqVar.getPosition();
        varPositions.add(position);
        final int readIndex = seqVar.getReadIndex();
        final int recodedReadIndex = entryOnReverseStrand ? readIndex - (queryLenth - position) + 5 : 5 + position - readIndex;

        // System.out.printf("%c CODING readIndex=%d position=%d queryLength=%d recodedReadIndex=%d %n",
        // entryOnReverseStrand ? '+' : '-', readIndex, position, queryLenth, recodedReadIndex);

        varReadIndex.add(recodedReadIndex);
        fromLengths.add(fromLength);
        toLengths.add(toLength);
        final int maxLength = Math.max(fromLength, toLength);
        for (int i = 0; i < maxLength; i++) {

            final char baseFrom = i < fromLength ? from.charAt(i) : '\0';
            final char baseTo = i < toLength ? to.charAt(i) : '\0';
            final byte byteFrom = (byte) baseFrom;
            final byte byteTo = (byte) baseTo;
            varFromTo.add(byteFrom << 8 | byteTo);

        }
    }

    MutableString from = new MutableString();
    MutableString to = new MutableString();


    private Alignments.AlignmentEntry andBack(final int index, int originalIndex, final Alignments.AlignmentEntry reduced, int streamVersion) {
        final Alignments.AlignmentEntry.Builder result = Alignments.AlignmentEntry.newBuilder(reduced);

        final int multiplicity = multiplicities.get(index);
        final int k = multiplicity - 1;

        multiplicities.set(index, k);
        //if (k > 1) {
        if (!multiplicityFieldsAllMissing) {
            result.setMultiplicity(1);
        }
        final int queryIndex = queryIndices.getInt(originalIndex);
        result.setQueryIndex(queryIndex);
        // System.out.printf("decoding query-index=%d (originalIndex=%d) varPositionIndex=%d %n",queryIndex,originalIndex, varPositionIndex);

        if (originalIndex == 0 || reduced.hasPosition() || reduced.hasTargetIndex()) {
            previousPosition = reduced.getPosition();
            previousTargetIndex = reduced.getTargetIndex();
        } else {


            final int deltaPos = deltaPositions.getInt(deltaPosIndex);
            final int deltaTarget = deltaTargetIndices.getInt(deltaPosIndex);
            final int position = previousPosition + deltaPos;
            final int targetIndex = previousTargetIndex + deltaTarget;
            result.setPosition(position);
            result.setTargetIndex(targetIndex);
            previousPosition += deltaPos;
            previousTargetIndex += deltaTarget;
            deltaPosIndex++;
        }
        if (streamVersion >= 2) {
            final int numReadQualScores = numReadQualityScores.get(index);
            if (numReadQualScores > 0) {

                final byte[] scores = new byte[numReadQualScores];
                for (int i = 0; i < numReadQualScores; i++) {
                    scores[i] = (byte) allReadQualityScores.getInt(qualScoreIndex++);
                }
                result.setReadQualityScores(ByteString.copyFrom(scores));

            }
        }
        int anInt = mappingQualities.getInt(index);
        if (anInt != MISSING_VALUE) {
            result.setMappingQuality(anInt);
        }
        anInt = fragmentIndices.getInt(index);
        if (anInt != MISSING_VALUE) {
            result.setFragmentIndex(anInt);
        }
        anInt = matchingReverseStrand.getInt(index);
        if (anInt != MISSING_VALUE) {
            result.setMatchingReverseStrand(anInt == 1);
        }
        anInt = numberOfMismatches.getInt(index);
        if (anInt != MISSING_VALUE) {
            result.setNumberOfMismatches(anInt);
        }

        anInt = numberOfIndels.getInt(index);
        if (anInt != MISSING_VALUE) {
            result.setNumberOfIndels(anInt);

        }
        final int queryLength = queryLengths.getInt(index);
        if (queryLength != MISSING_VALUE) {
            result.setQueryLength(queryLength);
        }
        anInt = queryAlignedLengths.getInt(index);
        if (anInt != MISSING_VALUE) {
            result.setQueryAlignedLength(anInt);
        }
        anInt = queryPositions.getInt(index);
        if (anInt != MISSING_VALUE) {
            result.setQueryPosition(anInt);
        }
        anInt = targetAlignedLengths.getInt(index);
        if (anInt != MISSING_VALUE) {
            result.setTargetAlignedLength(anInt);
        }
        anInt = sampleIndices.getInt(index);
        if (anInt != MISSING_VALUE) {
            result.setSampleIndex(anInt);
        }
        anInt = readOriginIndices.getInt(index);
        if (anInt != MISSING_VALUE) {
            result.setReadOriginIndex(anInt);
        }
        anInt = pairFlags.getInt(index);
        if (anInt != MISSING_VALUE) {
            result.setPairFlags(anInt);
        }
        anInt = scores.getInt(index);
        if (anInt != MISSING_VALUE) {
            result.setScore(Float.intBitsToFloat(anInt));
        }
        final boolean entryMatchingReverseStrand = result.hasMatchingReverseStrand() ? result.getMatchingReverseStrand() : false;
        Alignments.RelatedAlignmentEntry link = pairLinks.decode(originalIndex, entryMatchingReverseStrand, reduced.getPairAlignmentLink());
        if (link != null) {
            result.setPairAlignmentLink(link);
        }
        link = forwardSpliceLinks.decode(originalIndex, entryMatchingReverseStrand, reduced.getSplicedForwardAlignmentLink());
        if (link != null) {
            result.setSplicedForwardAlignmentLink(link);
        }
        link = backwardSpliceLinks.decode(originalIndex, entryMatchingReverseStrand, reduced.getSplicedBackwardAlignmentLink());
        if (link != null) {
            result.setSplicedBackwardAlignmentLink(link);
        }

        decodeInsertSize(result, index);

        final boolean templateHasSequenceVariations = reduced.getSequenceVariationsCount() > 0;
        final int numVariations = variationCount.getInt(index);

        for (int varIndex = 0; varIndex < numVariations; varIndex++) {
            final Alignments.SequenceVariation template = templateHasSequenceVariations ? reduced.getSequenceVariations(varIndex) : null;
            final Alignments.SequenceVariation.Builder varBuilder = templateHasSequenceVariations ?
                    Alignments.SequenceVariation.newBuilder(template) : Alignments.SequenceVariation.newBuilder();

            from.setLength(0);
            to.setLength(0);

            final int fromLength = fromLengths.getInt(varPositionIndex);
            final int toLength = toLengths.getInt(varPositionIndex);
            final int position = varPositions.getInt(varPositionIndex);
            varBuilder.setPosition(position);


            final int recodedReadIndex = varReadIndex.getInt(varPositionIndex);
            final int readIndex = entryMatchingReverseStrand ? recodedReadIndex + (queryLength - position) - 5 : -recodedReadIndex + position + 5;
            varBuilder.setReadIndex(readIndex);
            //  System.out.printf("%c DECODING position=%d queryLength=%d recodedReadIndex=%d readIndex=%d  %n",
            //         entryMatchingReverseStrand ? '+' : '-', position, queryLength, recodedReadIndex, readIndex);


            final int toQualLength = varToQualLength.getInt(varToQualLengthIndex);


            varToQualLengthIndex++;
            final byte[] quals = getQualArray(toQualLength);
            ++varPositionIndex;
            final int maxLength = Math.max(fromLength, toLength);
            for (int i = 0; i < maxLength; i++) {

                final int fromTo = varFromTo.getInt(varFromToIndex++);
                if (i < fromLength) {
                    from.append((char) (fromTo >> 8));
                }
                if (i < toLength) {
                    to.append((char) (fromTo & 0xFF));
                }
                if (i < toQualLength) {
                    if (varQualIndex < varQuals.size()) {

                        quals[i] = (byte) varQuals.getInt(varQualIndex);
                        ++varQualIndex;
                    }

                }
            }
            varBuilder.setFrom(from.toString());
            varBuilder.setTo(to.toString());
            if (toQualLength > 0) {
                varBuilder.setToQuality(ByteString.copyFrom(quals));
            }

            if (templateHasSequenceVariations) {
                result.setSequenceVariations(varIndex, varBuilder);
            } else {
                result.addSequenceVariations(varBuilder);
            }

        }
        if (result.hasReadQualityScores()) {

            final ByteString readQualScores = result.getReadQualityScores();
            // put toQual back on entries:
            for (int varIndex = 0; varIndex < numVariations; varIndex++) {

                final Alignments.SequenceVariation.Builder seqVarBuilder = result.getSequenceVariationsBuilder(varIndex);
                final String toBases = seqVarBuilder.getTo();

                final byte[] toQuals = new byte[toBases.length()];
                int indelOffset = 0;
                for (int l = 0; l < toBases.length(); ++l) {
                    final int i = l + seqVarBuilder.getReadIndex() - 1 - indelOffset;
                    final byte b = i >= readQualScores.size() ? 0 : readQualScores.byteAt(i);
                    final boolean ignoreBase = toBases.charAt(l) == '-';
                    toQuals[l] = ignoreBase ? 0 : b;

                    if (ignoreBase) {
                        indelOffset++;
                    }
                }
                seqVarBuilder.setToQuality(ByteString.copyFrom(toQuals));

                result.setSequenceVariations(varIndex, seqVarBuilder);
            }

        }
        return result.build();
    }

    /**
     * Decode the insert size given the already stored positions, an arithmetic expression linking position and insert size,
     * and any observed difference.
     *
     * @param result
     * @param index
     */
    private void decodeInsertSize(final Alignments.AlignmentEntry.Builder result, final int index) {
        final int anInt = insertSizes.getInt(index);
        if (anInt != MISSING_VALUE) {

            final int readPos = result.getPosition();
            final int matePos = result.getPairAlignmentLink().getPosition();
            final int length = result.getTargetAlignedLength();
            final int pos1 = result.getMatchingReverseStrand() ? length + readPos : readPos + 1;
            final int pos2 = EntryFlagHelper.isMateReverseStrand(result.getPairFlags()) ? length + matePos : matePos + 1;

            final int insertSize = pos2 - pos1 - anInt;
            // reverse:  insertSize= (pos2-pos1) - isd
            result.setInsertSize(insertSize);

        }
    }

    // pre-allocated arrays, size 1 to 100.
    byte[][] qualArrays = new byte[100][];

    private byte[] getQualArray(final int toQualLength) {
        return qualArrays[toQualLength];
    }

    private int varQualIndex = 0;
    private int varPositionIndex = 0;
    private int varFromToIndex = 0;


}
