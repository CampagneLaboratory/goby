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
import it.unimi.dsi.bits.Fast;
import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.fastutil.io.FastByteArrayInputStream;
import it.unimi.dsi.fastutil.objects.Object2IntAVLTreeMap;
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2LongAVLTreeMap;
import it.unimi.dsi.fastutil.objects.Object2LongMap;
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
    private static DynamicOptionClient doc = new DynamicOptionClient(AlignmentCollectionHandler.class,
            "stats-filename:string, the file where to append statistics to:",
            "debug-level:integer, a number between zero and 2. Numbers larger than zero activate debugging. 1 writes stats to stats-filename.:0",
            "basename:string, a basename for the file being converted.:"
    );
    private String statsFilename;
    private String basename;
    private PrintWriter statsWriter;
    private static final IntArrayList EMPTY_LIST = new IntArrayList();


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
    private int varToQualsLength = 0;
    private boolean useTemplateBasedCompression=true;
    private static final int LOG2_8 = Fast.mostSignificantBit(8);

    public AlignmentCollectionHandler() {
        for (int length = 0; length < qualArrays.length; length++) {
            qualArrays[length] = new byte[length];
        }
        debug = doc().getInteger("debug-level");
        statsFilename = doc().getString("stats-filename");
        basename = doc().getString("basename");
        if (debug(1)) {
            try {
                statsWriter = new PrintWriter(new FileWriter(statsFilename, true));
                statsWriter.print("basename\tchunkIndex\tlabel\tnumElements\ttotalBitsWritten\tBitsPerElement\n");

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
    public static final int VERSION = 3;

    @Override
    public Message compressCollection(final Message collection, final ByteArrayOutputStream compressedBits) throws IOException {
        reset();
        final Alignments.AlignmentCollection alignmentCollection = (Alignments.AlignmentCollection) collection;
        final Alignments.AlignmentCollection.Builder remainingCollection = Alignments.AlignmentCollection.newBuilder();
        final int size = alignmentCollection.getAlignmentEntriesCount();
        int indexInReducedCollection = 0;
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
            final int templateVarHasToQualsIndex = varToQualsLength;
            while (multiplicities.get(templateIndex) >= 1) {
                result.addAlignmentEntries(
                        andBack(templateIndex, originalIndex, alignmentCollection.getAlignmentEntries(templateIndex), streamVersion));
                if (multiplicities.get(templateIndex) >= 1) {
                    // go back to the indices for the template:
                    varPositionIndex = templatePositionIndex;
                    varFromToIndex = templateVarFromToIndex;
                    varToQualsLength = templateVarHasToQualsIndex;
                }
                originalIndex++;
            }
        }
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
        if (debug(1)) {
            //   out.flush();
            final long writtenStop = out.writtenBits();
            final long written = writtenStop - writtenStart;
            recordStats(label, list, written);
        }
    }

    /**
     * Write a list with Rice/Golomb coding
     * @param label
     * @param list
     * @param out
     * @throws IOException
     */
    public void writeRiceCoding(String label, IntList list, OutputBitStream out) throws IOException {
        final long writtenStart=out.writtenBits();
        out.writeNibble(list.size());
        for (final int value : list) {

            out.writeGolomb(value,8,LOG2_8);
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
        list.add(previous);
        decodeArithmetic(label, numEntriesInChunk - 1, bitInput, deltas);

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

        for (int i = 0; i < size; i++) {
            final int reducedReadIndex = bitInput.readMinimalBinary(max - min + 1);
            list.add(reducedReadIndex + min);
        }

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
        if (numEntriesInChunk == 0) {
            return;
        }
        // TODO see if we can avoid reading the number of elements in some cases.
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
        final FastArithmeticDecoder decoder = new FastArithmeticDecoder(numTokens);
        for (int i = 0; i < size; i++) {
            final int tokenValue = distinctvalue[decoder.decode(bitInput)];
            list.add(tokenValue);
        }
        decoder.reposition(bitInput);

    }

    protected final void writeArithmetic(final String label, final IntList list, OutputBitStream out) throws IOException {
        if (debug(2)) {
            System.err.flush();
            System.err.println("\nwriting " + label);
            System.err.flush();
        }
        final long writtenStart = out.writtenBits();
        out.writeNibble(list.size());        // LIST.size
        if (list.isEmpty()) {
            // no list to write.
            return;
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

        final FastArithmeticCoder coder = new FastArithmeticCoder(distinctSymbols.size());
        for (final int dp : list) {
            final int symbolCode = Arrays.binarySearch(symbolValues, dp);
            assert symbolCode >= 0 : "symbol code must exist.";
            coder.encode(symbolCode, out);
        }
        coder.flush(out);
        if (debug(1)) {
            System.err.flush();
            final long writtenStop = out.writtenBits();
            final long written = writtenStop - writtenStart;
            recordStats(label, list, written);
        }
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

    IntSortedSet tokenSet = new IntAVLTreeSet();

    private IntSortedSet getTokens(IntList list) {
        tokenSet.clear();
        for (int value : list) {
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
    private IntArrayList numReadQualityScores = new IntArrayList();

    private IntArrayList multiplicities = new IntArrayList();


    private int decompressBits(InputBitStream bitInput, final int numEntriesInChunk) throws IOException {
        final int streamVersion = bitInput.readDelta();
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
        return streamVersion;
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
        writeArithmetic("allReadQualityScores", allReadQualityScores, out);
        writeArithmetic("sampleIndices", sampleIndices, out);
        writeArithmetic("readOriginIndices", readOriginIndices, out);
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
        varToQualsLength = 0;
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

        recordVariationQualitiesAndClear(result, result.getSequenceVariationsList());

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
        insertSizes.add(source.hasInsertSize() ? source.getInsertSize() : MISSING_VALUE);
        queryAlignedLengths.add(source.hasQueryAlignedLength() ? source.getQueryAlignedLength() : MISSING_VALUE);
        targetAlignedLengths.add(source.hasTargetAlignedLength() ? source.getTargetAlignedLength() : MISSING_VALUE);
        fragmentIndices.add(source.hasFragmentIndex() ? source.getFragmentIndex() : MISSING_VALUE);
        variationCount.add(source.getSequenceVariationsCount());
        queryPositions.add(source.hasQueryPosition() ? source.getQueryPosition() : MISSING_VALUE);
        sampleIndices.add(source.hasSampleIndex() ? source.getSampleIndex() : MISSING_VALUE);
        readOriginIndices.add(source.hasReadOriginIndex() ? source.getReadOriginIndex() : MISSING_VALUE);

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
        //     print(result);


        final Alignments.AlignmentEntry alignmentEntry = result.build();
  //    System.out.println(alignmentEntry);
        return alignmentEntry;
    }

    private boolean fastEquals(final Object o1, final Object o2) {
        return useTemplateBasedCompression && o1.equals(o2);

    }


    private void recordVariationQualitiesAndClear(Alignments.AlignmentEntry.Builder result, List<Alignments.SequenceVariation> sequenceVariationsList) {

        int index = 0;
        for (final Alignments.SequenceVariation seqVar : sequenceVariationsList) {
            final String from = seqVar.getFrom();

            final ByteString toQualities = seqVar.getToQuality();
            final boolean hasToQuals = seqVar.hasToQuality();
            final int toQualSize = toQualities.size();
            varToQualLength.add(toQualSize);
            final int toQualsLength = hasToQuals ? seqVar.getToQuality().size() : 0;
            for (int i = 0; i < toQualsLength; i++) {
                if (hasToQuals && i < toQualSize) {
                    varQuals.add(toQualities.byteAt(i));
                }
            }
            Alignments.SequenceVariation.Builder varBuilder = Alignments.SequenceVariation.newBuilder(seqVar);
            varBuilder.clearToQuality();
            result.setSequenceVariations(index, varBuilder.buildPartial());
            index++;
        }
    }

    private void print(Alignments.AlignmentEntry result) {

        System.out.println(result);
    }

    private void encodeVar(boolean entryOnReverseStrand, int queryLenth, final Alignments.SequenceVariation seqVar) {

        final String from = seqVar.getFrom();
        final String to = seqVar.getTo();
        final ByteString toQualities = seqVar.getToQuality();
        final int fromLength = from.length();
        final int toLength = to.length();
        final boolean hasToQuals = seqVar.hasToQuality();
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
        anInt = insertSizes.getInt(index);
        if (anInt != MISSING_VALUE) {
            result.setInsertSize(anInt);
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

            final int toQualLength = varToQualLength.getInt(varToQualsLength);
            varToQualsLength++;

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
                        //  if (quals[i] != NO_VALUE) {
                        ++varQualIndex;
                    }
                    // }
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
        return result.build();
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
