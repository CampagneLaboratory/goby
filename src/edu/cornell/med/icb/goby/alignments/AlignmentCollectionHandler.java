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
import com.google.protobuf.CodedInputStream;
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
import java.util.Date;
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
            "ignore-read-origin:boolean, When this flag is true do not compress read origin/read groups.:false",
            "symbol-dependency-order:integer, The number of past symbols to consider when estimating compression models (0 or 1):0",
            "enable-domain-optimizations:boolean, When this flag is true we use compression methods that are domain specific, and can increase further compression. For instance, setting this flag to true will compress related-alignment-links very efficiently if they link entries in the same chunk.:true"

    );
    private String statsFilename;
    private String basename;
    private PrintWriter statsWriter;
    private static final IntArrayList EMPTY_LIST = new IntArrayList();
    private boolean storeReadOrigins = true;
    private final boolean useArithmeticCoding = true;
    private final boolean useHuffmanCoding = false;

    private int numReadQualScoresIndex;
    public boolean enableDomainOptimizations = false;
    private boolean useTemplateBasedCompression = true;
    private int linkOffsetOptimizationIndex;
    /* Special symbol used to indicate we need to reset previousVar position to 1.*/

    private static final int RESET_VAR_POS = 0;
    private int insertSizeIndex = 0;
    /**
     * Time stamp determined at the time the stats file is opened. Identifies all the stats written in the same
     * execution run.
     */
    private long timeStamp;

    /**
     * Specify the order (0 or 1) for modeling symbol dependencies.
     * @param order 1 (probability of a symbol depends on the previous symbol and the current symbol) 0 (only current symbol considered)
     */
    public void setSymbolDependencyOrder(int order) {
        this.encodeSymbolDependencyOrder = order;
    }

    /**
     * The number of previous symbols to consider when estimating a model probability.
     */
    private int encodeSymbolDependencyOrder = 0;
    private int decodeSymbolDependencyOrder;

    public boolean isEnableDomainOptimizations() {
        return enableDomainOptimizations;
    }

    public void setEnableDomainOptimizations(boolean enableDomainOptimizations) {
        this.enableDomainOptimizations = enableDomainOptimizations;
    }

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
    static final int MISSING_VALUE = -1;
    private boolean multiplicityFieldsAllMissing = true;
    private long writtenBits;
    private long writtenBases;
    private static final int NO_VALUE = MISSING_VALUE;
    private int varToQualLengthIndex = 0;

    private static final int LOG2_8 = Fast.mostSignificantBit(8);

    public AlignmentCollectionHandler() {
        for (int length = 0; length < qualArrays.length; length++) {
            qualArrays[length] = new byte[length];
        }
        debug = doc().getInteger("debug-level");
        storeReadOrigins = !doc().getBoolean("ignore-read-origin");
        enableDomainOptimizations = doc().getBoolean("enable-domain-optimizations");
        statsFilename = doc().getString("stats-filename");
        basename = doc().getString("basename");
        encodeSymbolDependencyOrder = doc().getInteger("symbol-dependency-order");

        if (debug(1)) {
            try {
                final boolean appending = new File(statsFilename).exists();
                final FileWriter fileWriter = appending ? new FileWriter(statsFilename, true) : new FileWriter(statsFilename);
                timeStamp = new Date().getTime();
                statsWriter = new PrintWriter(fileWriter);
                if (!appending) {
                    statsWriter.print("timeStamp\tbasename\tchunkIndex\tlabel\tnumElements\ttotalBitsWritten\tBitsPerElement\tPercentageOfCompressed\n");
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
        final CodedInputStream codedInput = CodedInputStream.newInstance(uncompressedStream);
        codedInput.setSizeLimit(Integer.MAX_VALUE);

        return Alignments.AlignmentCollection.parseFrom(codedInput);
    }

    int numChunksProcessed = 0;
    /**
     * The version of the stream that this class reads and writes.
     */

    public static final int VERSION = 13;
    private int streamVersion;

    @Override
    public Message compressCollection(final Message collection, final ByteArrayOutputStream compressedBits) throws IOException {
        reset();
        final Alignments.AlignmentCollection alignmentCollection = (Alignments.AlignmentCollection) collection;
        final Alignments.AlignmentCollection.Builder remainingCollection = Alignments.AlignmentCollection.newBuilder();
        final int size = alignmentCollection.getAlignmentEntriesCount();
        int indexInReducedCollection = 0;
        collectLinkLists(alignmentCollection);
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

    Int2ObjectMap<IntArrayList> queryIndexToPositionList = new Int2ObjectOpenHashMap<IntArrayList>();
    Int2ObjectMap<IntArrayList> queryIndex2EntryIndices = new Int2ObjectOpenHashMap<IntArrayList>();
    Int2ObjectMap<IntArrayList> queryIndex2FragmentIndices = new Int2ObjectOpenHashMap<IntArrayList>();

    private void collectLinkLists(Alignments.AlignmentCollectionOrBuilder alignmentCollection) {
        int entryIndex = 0;
        // collect positions for each query index:
        for (final Alignments.AlignmentEntry entry : alignmentCollection.getAlignmentEntriesList()) {

            final int queryIndex = entry.getQueryIndex();
            IntArrayList list = queryIndexToPositionList.get(queryIndex);
            IntArrayList listOfIndices = queryIndex2EntryIndices.get(queryIndex);
            IntArrayList listOfFragmentIndices = queryIndex2FragmentIndices.get(queryIndex);
            if (list == null) {
                list = new IntArrayList();
                listOfIndices = new IntArrayList();
                listOfFragmentIndices = new IntArrayList();
                queryIndexToPositionList.put(queryIndex, list);
                queryIndex2EntryIndices.put(queryIndex, listOfIndices);
                queryIndex2FragmentIndices.put(queryIndex, listOfFragmentIndices);
            }


            list.add(entry.getPosition());
            listOfIndices.add(entryIndex);
            listOfFragmentIndices.add(entry.getFragmentIndex());
            entryIndex++;


        }

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
                        if (entry.hasSoftClippedQualityLeft()) {
                            softClipLeftQualityScores.add(entry.getSoftClippedQualityLeft().byteAt(indexInStrings));
                        }
                        finished1 = false;
                    }
                }
                final int numSoftClipRightBasesInt = numSoftClipRightBases.getInt(index);
                if (numSoftClipRightBasesInt != MISSING_VALUE) {
                    if (indexInStrings < numSoftClipRightBasesInt) {
                        final String softClippedBasesRight = entry.getSoftClippedBasesRight();

                        softClipRightBases.add(softClippedBasesRight.charAt(indexInStrings));
                        if (entry.hasSoftClippedQualityRight()) {
                            softClipRightQualityScores.add(entry.getSoftClippedQualityRight().byteAt(indexInStrings));
                        }
                        finished2 = false;
                    }
                }
            }
            indexInStrings++;
        }
        // must clear the attributes we have collected, but cannot do this on the read-only PB entries.
    }

    private void restoreStrings(Alignments.AlignmentCollection.Builder alignmentCollection) {
        final int size = alignmentCollection.getAlignmentEntriesCount();
        int indexInStrings = 0;
        boolean finished1 = false;
        boolean finished2 = false;
        final ObjectArrayList<MutableString> softClipsLeft = new ObjectArrayList<MutableString>();
        final ObjectArrayList<MutableString> softClipsRight = new ObjectArrayList<MutableString>();
        final ObjectArrayList<byte[]> softClipsLeftQual = new ObjectArrayList<byte[]>();
        final ObjectArrayList<byte[]> softClipsRightQual = new ObjectArrayList<byte[]>();
        softClipsLeft.size(size);
        softClipsLeftQual.size(size);
        softClipsRight.size(size);
        softClipsRightQual.size(size);
        for (int index = 0; index < size; index++) {

            if (index < numSoftClipLeftBases.size()) {
                final int stringLength = numSoftClipLeftBases.getInt(index);
                if (stringLength != MISSING_VALUE) {

                    final MutableString mutableString = new MutableString();
                    mutableString.setLength(stringLength);

                    softClipsLeft.set(index, mutableString);
                    final byte[] leftQuals = new byte[stringLength];
                    softClipsLeftQual.set(index, leftQuals);
                }
            }

            if (index < numSoftClipRightBases.size()) {
                final int stringLength = numSoftClipRightBases.getInt(index);
                if (stringLength != MISSING_VALUE) {

                    final MutableString mutableString = new MutableString();
                    mutableString.setLength(stringLength);
                    softClipsRight.set(index, mutableString);
                    final byte[] rightQuals = new byte[stringLength];
                    softClipsRightQual.set(index, rightQuals);
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
                        final int anInt = softClipLeftBases.getInt(iLeft);
                        mutableString.setCharAt(indexInStrings, (char) anInt);
                        if (streamVersion >= 10) {
                            final byte qual = (byte) softClipLeftQualityScores.getInt(iLeft);
                            softClipsLeftQual.get(index)[indexInStrings] = qual;
                        }
                        iLeft += 1;
                        finished1 = false;
                    }
                }
                final int numSoftClipRightBasesInt = numSoftClipRightBases.getInt(index);

                if (numSoftClipRightBasesInt != MISSING_VALUE && indexInStrings < numSoftClipRightBasesInt) {
                    final MutableString mutableString = softClipsRight.get(index);
                    if (mutableString != null) {
                        final int anInt = softClipRightBases.getInt(iRight);
                        mutableString.setCharAt(indexInStrings, (char) anInt);
                        if (streamVersion >= 10) {
                            final byte qual = (byte) softClipRightQualityScores.getInt(iRight);
                            softClipsRightQual.get(index)[indexInStrings] = qual;
                        }
                        iRight += 1;
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
                    if (streamVersion >= 10) {
                        builder.setSoftClippedQualityLeft(ByteString.copyFrom(softClipsLeftQual.get(index)));
                    }
                }
            }
            {
                final MutableString mutableString = softClipsRight.get(index);
                if (mutableString != null) {
                    builder.setSoftClippedBasesRight(mutableString.toString());
                    if (streamVersion >= 10) {
                        builder.setSoftClippedQualityRight(ByteString.copyFrom(softClipsRightQual.get(index)));
                    }
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
        restoreLinks(result);
        ++chunkIndex;
        return result.build();
    }

    private void restoreLinks(Alignments.AlignmentCollection.Builder alignmentCollection) {
        queryIndexToPositionList.clear();
        collectLinkLists(alignmentCollection);

        final int size = alignmentCollection.getAlignmentEntriesCount();
        insertSizeIndex = 0;
        for (int index = 0; index < size; index++) {
            final Alignments.AlignmentEntry.Builder entry = alignmentCollection.getAlignmentEntriesBuilder(index);
            final IntArrayList positionList = queryIndexToPositionList.get(entry.getQueryIndex());
            final IntArrayList fragmentList = queryIndex2FragmentIndices.get(entry.getQueryIndex());
            final int queryIndex = entry.getQueryIndex();
            if (entry.hasPairAlignmentLink() && entry.getPairAlignmentLink().hasOptimizedIndex()) {

                final Alignments.RelatedAlignmentEntry.Builder linkBuilder = entry.getPairAlignmentLinkBuilder();
                recoverLink(alignmentCollection, entry, positionList, fragmentList, queryIndex, linkBuilder);
            }

            if (entry.hasSplicedForwardAlignmentLink() && entry.getSplicedForwardAlignmentLink().hasOptimizedIndex()) {

                final Alignments.RelatedAlignmentEntry.Builder linkBuilder = entry.getSplicedForwardAlignmentLinkBuilder();
                recoverLink(alignmentCollection, entry, positionList, fragmentList, queryIndex, linkBuilder);
            }
            if (entry.hasSplicedBackwardAlignmentLink() && entry.getSplicedBackwardAlignmentLink().hasOptimizedIndex()) {

                final Alignments.RelatedAlignmentEntry.Builder linkBuilder = entry.getSplicedBackwardAlignmentLinkBuilder();
                recoverLink(alignmentCollection, entry, positionList, fragmentList, queryIndex, linkBuilder);
            }
            // we need to update insert size because the optimization messed up the mate position used by domain optimization:

            recalculateInsertSize(entry, insertSizeIndex++);
        }
    }

    private void recoverLink(final Alignments.AlignmentCollection.Builder alignmentCollection, final Alignments.AlignmentEntry.Builder entry,
                             final IntArrayList positionList, final IntArrayList fragmentList, final int queryIndex,
                             final Alignments.RelatedAlignmentEntry.Builder linkBuilder) {
        final int indexOffset = linkBuilder.getOptimizedIndex();
        final int position = entry.getPosition();
        int index = 0;
        int thisEntryIndex = -1;
        for (final int value : positionList) {
            if (value == position && fragmentList.get(index) == entry.getFragmentIndex()) {
                thisEntryIndex = index;
            }
            index++;
        }
        final int optimizedIndex = indexOffset + thisEntryIndex;

        //   final int linkIndex = Fast.nat2int(linkOffset) + thisEntryIndex;
        final IntArrayList listOfIndices = queryIndex2EntryIndices.get(queryIndex);
        final int otherEntryIndex = listOfIndices.get(optimizedIndex);
        final Alignments.AlignmentEntry otherEntry = alignmentCollection.getAlignmentEntries(otherEntryIndex);
        linkBuilder.setPosition(otherEntry.getPosition());
        linkBuilder.setFragmentIndex(otherEntry.getFragmentIndex());
        linkBuilder.clearOptimizedIndex();
    }

    int findIndex(IntSortedSet list, int position) {
        int index = 0;
        for (int value : list.toIntArray()) {
            if (value == position) return index;
            index += 1;
        }
        return -1;

    }


    @Override
    public void setUseTemplateCompression(final boolean useTemplateCompression) {
        useTemplateBasedCompression = useTemplateCompression;
    }


    public void displayStats() {
        if (debug(1) && statsWriter != null) {

            final double overallBpb = divide(writtenBits, writtenBases);
            double totalBpb = overallBpb;
            long sumN = 0;
            long sumWritten = 0;
            for (String label : typeToNumEntries.keySet()) {
                int n = typeToNumEntries.getInt(label);
                long written = typeToWrittenBits.getLong(label);
                sumWritten += written;
                sumN += n;
            }
            for (String label : typeToNumEntries.keySet()) {
                int n = typeToNumEntries.getInt(label);
                long written = typeToWrittenBits.getLong(label);
                double average = (double) written / (double) n;
                double averageBitPerBase = (double) written / (double) writtenBases;
                // LOG.info
                //       (String.format("encoded %d %s in %d bits, average %g bits /element. ", n, label,
                //             written, average));
                final double percentageOfTotal = divide(written, sumWritten) * 100.0;
                LOG.info
                        (String.format("encoded %d %s in %d bits, average %g bpb. %2.2f%% ", n, label,
                                written, averageBitPerBase, percentageOfTotal));
                statsWriter.write(String.format("%d\t%s\t%d\t%s\t%d\t%d\t%g\t%.2g%n", timeStamp, basename, chunkIndex, label, n, written,
                        divide(written, n), percentageOfTotal));
                totalBpb += averageBitPerBase;

            }
            LOG.info(String.format("entries aggregated with multiplicity= %d", countAggregatedWithMultiplicity));
            LOG.info(String.format("Overall: bits per aligned bases= %g", overallBpb));
            statsWriter.write(String.format("%d\t%s\t%d\t%s\t%d\t%d\t%g\t%d%n", timeStamp, basename, chunkIndex, "overall-bits-per-base", sumN, sumWritten,
                    overallBpb, 100));
            statsWriter.flush();

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
        final boolean success = tryWriteDeltas(label, list, out);
        final long writtenStart = out.writtenBits();

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
        final long writtenStart = out.writtenBits();

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
            writeArithmetic(null, deltas, out);
            final long writtenStop = out.writtenBits();
            final long written = writtenStop - writtenStart;
            recordStats(label + "-as-deltas", deltas, written);
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
            encodedLengths.clear();
            encodedValues.clear();
            decodeArithmeticInternal(bitInput, encodedLengths);
            decodeArithmeticInternal(bitInput, encodedValues);

            decodeRunLengths(encodedLengths, encodedValues, list);
        } else {
            decodeArithmeticInternal(bitInput, list);
        }


    }

    final IntArrayList encodedLengths = new IntArrayList();
    final IntArrayList encodedValues = new IntArrayList();

    private void decodeArithmeticInternal(InputBitStream bitInput, IntList list) throws IOException {
        final int size = bitInput.readNibble();
        if (size == 0) {
            return;
        }
        final boolean directEncoding = streamVersion >= 12 ? bitInput.readBit() == 1 : false;
        final boolean hasNegatives = bitInput.readBit() == 1;

        final int numTokens = bitInput.readNibble();
        final int[] distinctvalue = new int[numTokens];
        if (!directEncoding) {

            for (int i = 0; i < numTokens; i++) {
                // -1 makes 0 symbol -1 (missing value) again
                final int token = bitInput.readNibble();
                final int anInt = hasNegatives ? Fast.nat2int(token) : token - 1;
                distinctvalue[i] = anInt;
            }
        } else {
            // read the first token and reconstruct the others in sequence:
            final int token = bitInput.readNibble();
            final int anInt = hasNegatives ? Fast.nat2int(token) : token - 1;
            distinctvalue[0] = anInt;
            for (int i = 1; i < numTokens; i++) {
                distinctvalue[i] = distinctvalue[i - 1] + 1;
            }
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
        encodedLengths.clear();
        encodedValues.clear();
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
        float ratioListSizes = ((float) encodedLengths.size() + encodedValues.size()) / (float) list.size();
        final boolean result = encodedLengths.size() > 10 && ratioListSizes < 1;
        /*
        if (result) {
                  System.out.println("Using run-length encoding for "+label + " encodedLengths.size() "+encodedLengths.size() );
        } */
        return result;
    }

    private boolean encodeArithmeticInternal(String label, IntList list, OutputBitStream out) throws IOException {
        out.writeNibble(list.size());        // LIST.size
        if (list.isEmpty()) {
            // no list to write.
            return true;
        }
        final IntSortedSet distinctSymbols = getTokens(list);
        int minSymbol = distinctSymbols.firstInt();
        int maxSymbol = distinctSymbols.lastInt();
        int symbolRange = maxSymbol - minSymbol;
        int numDistinctSymbols = distinctSymbols.size();
        final boolean hasNegativeValues = minSymbol < 0;
        if (symbolRange + 1 > numDistinctSymbols) {
            final int[] symbolValues = distinctSymbols.toIntArray();
            // new in version 12:
            out.writeBit(false);
            out.writeBit(hasNegativeValues);
            out.writeNibble(distinctSymbols.size());    //LIST.distinct-values.size()
            for (final int token : distinctSymbols) {
                // +1 makes -1 (missing value) symbol 0 so it can be written Nibble:

                final int anInt = token;
                final int nat = hasNegativeValues ? Fast.int2nat(anInt) : anInt + 1;

                out.writeNibble(nat);
            }

            encode(label, list, out, distinctSymbols, symbolValues,true,-1);
            return false;
        } else {
            //  System.out.printf("%s symbolRange=%d symbolSize=%d%n", label, symbolRange, numDistinctSymbols);

            // in this branch, we know the symbols are consecutive in values, so we write the min symbol, and the number of symbols.
            out.writeBit(true);
            out.writeBit(hasNegativeValues);
            out.writeNibble(distinctSymbols.size());    //LIST.distinct-values.size()

            // +1 makes -1 (missing value) symbol 0 so it can be written Nibble:

            final int anInt = minSymbol;
            final int nat = hasNegativeValues ? Fast.int2nat(anInt) : anInt + 1;

            out.writeNibble(nat);
            final int[] symbolValues = distinctSymbols.toIntArray();
            encode(label, list, out, distinctSymbols, symbolValues, false, minSymbol);
            return false;
        }
    }

    private void encode(String label, final IntList list, final OutputBitStream out, final IntSet distinctSymbols, final int[] symbolValues,final boolean useBinarySearch, final int minSymbol) throws IOException {
        if (useArithmeticCoding) {
            final int numSymbols = distinctSymbols.size();
            int thisSymbolDependencyOrder= encodeSymbolDependencyOrder;
            switch (thisSymbolDependencyOrder) {
                case 0:
                    final FastArithmeticCoder coder = new FastArithmeticCoder(numSymbols);
                    for (final int dp : list) {
                        final int symbolCode =useBinarySearch? Arrays.binarySearch(symbolValues, dp) : dp-minSymbol;
                        assert symbolCode >= 0 : "symbol code must exist.";
                        coder.encode(symbolCode, out);
                    }
                    coder.flush(out);
                    break;
                case 1:
                    final FastArithmeticCoder[] coders = new FastArithmeticCoder[numSymbols];
                    for (int i = 0; i < numSymbols; i++) {
                        coders[i] = new FastArithmeticCoder(numSymbols);
                    }
                    int last = 0;
                    for (final int dp : list) {

                        final int symbolCode = useBinarySearch? Arrays.binarySearch(symbolValues, dp) : dp-minSymbol;
                        assert symbolCode >= 0 : "symbol code must exist.";

                        coders[last].encode(symbolCode, out);
                        last = symbolCode;
                    }
                    for (int i = 0; i < numSymbols; i++) {
                        coders[i].flush(out);
                    }
            }
        } else if (useHuffmanCoding) {
            final int[] frequencies = frequencies(list, symbolValues);
            final HuffmanCodec codec = new HuffmanCodec(frequencies);
            final CodeWordCoder coder = codec.coder();
            for (int freq : frequencies) {
                out.writeNibble(freq);
            }
            for (final int dp : list) {
                final int symbolCode = useBinarySearch? Arrays.binarySearch(symbolValues, dp) : dp-minSymbol;
                assert symbolCode >= 0 : "symbol code must exist.";
                coder.encode(symbolCode, out);
            }
            coder.flush(out);
        }
    }

    /**
     * Here, we know that symbol values are consecutive, so we don't need to use binarySearch to code list values into symbols.
     * We use direct access instead.
     *
     * @param label
     * @param list
     * @param out
     * @throws IOException
     */
    private void encodeDirect(String label, final IntList list, final OutputBitStream out, int minSymbol, int numSymbols) throws IOException {
     //TODO add symbolDependencyOrder
        if (useArithmeticCoding) {
            final FastArithmeticCoder coder = new FastArithmeticCoder(numSymbols);
            for (final int dp : list) {
                final int symbolCode = dp - minSymbol;
                assert symbolCode >= 0 : "symbol code must exist.";
                coder.encode(symbolCode, out);
            }
            coder.flush(out);
        } else if (useHuffmanCoding) {
            final int[] frequencies = frequenciesDirect(list, minSymbol, numSymbols);
            final HuffmanCodec codec = new HuffmanCodec(frequencies);
            final CodeWordCoder coder = codec.coder();
            for (int freq : frequencies) {
                out.writeNibble(freq);
            }
            for (final int dp : list) {
                final int symbolCode = dp - minSymbol;
                assert symbolCode >= 0 : "symbol code must exist.";
                coder.encode(symbolCode, out);
            }
            coder.flush(out);
        }
    }

    private void decode(final InputBitStream bitInput, final IntList list, final int size, final int numTokens, int[] distinctvalue) throws IOException {
        if (useArithmeticCoding) {
            switch (decodeSymbolDependencyOrder) {
                case 0:
                    final FastArithmeticDecoder decoder = new FastArithmeticDecoder(numTokens);
                    for (int i = 0; i < size; i++) {
                        final int tokenValue = distinctvalue[decoder.decode(bitInput)];
                        list.add(tokenValue);
                    }
                    decoder.reposition(bitInput);
                    break;
                case 1:
                    final FastArithmeticDecoder decoders[] = new FastArithmeticDecoder[numTokens];
                    for (int i = 0; i < numTokens; i++) {
                        decoders[i] = new FastArithmeticDecoder(numTokens);
                    }
                    int last = 0;
                    for (int i = 0; i < size; i++) {
                        final int tokenValue = distinctvalue[decoders[last].decode(bitInput)];
                        list.add(tokenValue);
                        last = tokenValue;
                    }
                    for (int i = 0; i < numTokens; i++) {
                        decoders[i].reposition(bitInput);
                    }

            }
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

    // return the frequencies of symbols in the list
    private int[] frequenciesDirect(IntList list, int minSymbol, int numSymbols) {
        final int[] freqs = new int[numSymbols];
        for (final int value : list) {
            final int index = value - minSymbol;
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
        if (debug(1) && label != null) {
            double average = ((double) written) / list.size();
            typeToNumEntries.put(label, list.size() + typeToNumEntries.getInt(label));
            typeToWrittenBits.put(label, written + typeToWrittenBits.getLong(label));
            /* if (label.contains("delta"))
               System.out.printf("written: %d %n", written);
            */
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
    private IntList softClipLeftQualityScores = new IntArrayList();
    private IntList softClipRightQualityScores = new IntArrayList();

    protected IntList linkOffsetOptimization = new IntArrayList();


    private int decompressBits(InputBitStream bitInput, final int numEntriesInChunk) throws IOException {
        streamVersion = bitInput.readDelta();
        if (streamVersion<=12) {
            decodeSymbolDependencyOrder=0;
        } else {
            decodeSymbolDependencyOrder= bitInput.readDelta();
        }
        assert streamVersion <= VERSION : "FATAL: The stream version cannot have been written with a more recent version of Goby (The hybrid chunk codec cannot not support forward compatibility of the compressed stream).";
        if (streamVersion >= 8) {
            enableDomainOptimizations = bitInput.readBit() == 1;
        }
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
        if (streamVersion >= 10) { // in version 10 we try to store positions as deltas:
            decodeVarPositions("varPositions-as-deltamod", numEntriesInChunk, bitInput, varPositions);
        } else {
            decodeArithmetic("varPositions", numEntriesInChunk, bitInput, varPositions);
        }
        decodeArithmetic("fromLengths", numEntriesInChunk, bitInput, fromLengths);
        decodeArithmetic("toLengths", numEntriesInChunk, bitInput, toLengths);
        if (enableDomainOptimizations && streamVersion >= 10) {
            decodeToLength(fromLengths, toLengths);
        }
        decodeArithmetic("varReadIndex", numEntriesInChunk, bitInput, varReadIndex);
        decodeArithmetic("varFromTo", numEntriesInChunk, bitInput, varFromTo);
        decodeArithmetic("varQuals", numEntriesInChunk, bitInput, varQuals);
        decodeArithmetic("varToQualLength", numEntriesInChunk, bitInput, varToQualLength);
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
        if (streamVersion >= 7) {
            decodeArithmetic("linkOffsetOptimization", numEntriesInChunk, bitInput, linkOffsetOptimization);
        }
        if (streamVersion >= 9) {
            decodeArithmetic("softClipLeftQualityScores", numEntriesInChunk, bitInput, softClipLeftQualityScores);
            decodeArithmetic("softClipRightQualityScores", numEntriesInChunk, bitInput, softClipRightQualityScores);
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
        out.writeDelta(encodeSymbolDependencyOrder);
        out.writeBit(enableDomainOptimizations);
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
        writeVarPositions("varPositions", varPositions, out);
        writeArithmetic("fromLengths", fromLengths, out);
        if (enableDomainOptimizations) {
            encodeToLength(fromLengths, toLengths);
        }
        writeArithmetic("toLengths", toLengths, out);
        writeArithmetic("varReadIndex", varReadIndex, out);
        writeArithmetic("varFromTo", varFromTo, out);
        writeArithmetic("varQuals", varQuals, out);
        writeArithmetic("varToQualLength", varToQualLength, out);
        writeArithmetic("multiplicities", multiplicities, out);
        pairLinks.write(out);
        forwardSpliceLinks.write(out);
        backwardSpliceLinks.write(out);

        writeQueryIndices("queryIndices", queryIndices, out);

        writeArithmetic("numReadQualityScores", numReadQualityScores, out);
        writeArithmetic("allReadQualityScores", allReadQualityScores, out);

        writeArithmetic("sampleIndices", sampleIndices, out);
        writeArithmetic("readOriginIndices", readOriginIndices, out);
        writeArithmetic("pairFlags", pairFlags, out);
        writeArithmetic("scores", scores, out);

        writeArithmetic("softClipLeftBasesNum", numSoftClipLeftBases, out);
        writeArithmetic("softClipRightBasesNum", numSoftClipRightBases, out);
        writeArithmetic("softClipLeftBases", softClipLeftBases, out);
        writeArithmetic("softClipRightBases", softClipRightBases, out);

        writeArithmetic("linkOffsetOptimization", linkOffsetOptimization, out);
        writeArithmetic("softClipLeftQualityScores", softClipLeftQualityScores, out);
        writeArithmetic("softClipRightQualityScores", softClipRightQualityScores, out);
    }

    /**
     * Model toLength as a function of fromLengths. Values are recoded in place.
     *
     * @param fromLengths
     * @param toLengths
     */
    public void encodeToLength(final IntList fromLengths, final IntList toLengths) {
        final int min = Math.min(fromLengths.size(), toLengths.size());
        for (int i = 0; i < min; i++) {
            // diff=fromLength -toLength
            final int toLength = toLengths.getInt(i);
            final int fromLength = fromLengths.getInt(i);
            final int diff = Fast.int2nat(fromLength - toLength);
            toLengths.set(i, diff);
        }
    }

    /**
     * Model toLength as a function of fromLengths. Values are recoded in place.
     *
     * @param fromLengths
     * @param toLengths
     */
    public void decodeToLength(final IntList fromLengths, final IntList toLengths) {
        for (int i = 0; i < fromLengths.size(); i++) {

            final int diff = Fast.nat2int(toLengths.getInt(i));
            final int fromLength = fromLengths.getInt(i);
            final int toLength = -(diff - fromLength);
            toLengths.set(i, toLength);
        }
    }

    final protected IntArrayList varPositionDeltaMods = new IntArrayList();

    private void writeVarPositions(final String label, final IntList varPositionsList, final OutputBitStream out) throws IOException {
        if (!enableDomainOptimizations) {
            writeArithmetic(label, varPositionsList, out);
        } else {
            deltaModTransform(varPositionsList);
            writeArithmetic(label + "-as-deltamod", varPositionDeltaMods, out);
        }
    }

    public IntArrayList deltaModTransform(IntList varPositionsList) {
        varPositionDeltaMods.clear();
        int previous = 0;

        for (final int value : varPositionsList) {

            final int pos_plus_one = value + 1;
            final int diff = pos_plus_one - previous;
            if (diff > 0) {
                varPositionDeltaMods.add(diff);
            } else {
                varPositionDeltaMods.add(RESET_VAR_POS);
                varPositionDeltaMods.add(pos_plus_one);
                previous = 0;

            }
            previous = value;

        }
        return varPositionDeltaMods;
    }

    protected final void decodeVarPositions(final String label, final int numEntriesInChunk, final InputBitStream bitInput, final IntList varPositionsList) throws IOException {
        if (!enableDomainOptimizations) {
            decodeArithmetic(label, numEntriesInChunk, bitInput, varPositionsList);
            return;
        }
        if (debug(2)) {
            System.err.flush();
            System.err.println("\nreading " + label + " with available=" + bitInput.available());
            System.err.flush();
        }
        varPositionDeltaMods.clear();
        decodeArithmetic(label, numEntriesInChunk, bitInput, varPositionDeltaMods);
        decodeDeltaModTransform(varPositionDeltaMods, varPositionsList);

    }

    public void decodeDeltaModTransform(final IntList input, final IntList varPositionsList) {
        int previous = 0;
        for (final int diff : input) {
            final int pos_plus_one = diff + previous;
            final int value = pos_plus_one - 1;
            if (diff == RESET_VAR_POS) {
                previous = 0;
                continue;
            } else {

                previous = value;
            }
            varPositionsList.add(value);
        }
        // varPositionsList.add(previous);
    }

    private void reset() {

        insertSizeIndex = 0;
        multiplicityFieldsAllMissing = true;
        queryIndexToPositionList.clear();
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
        numReadQualScoresIndex = 0;
        allReadQualityScores.clear();
        sampleIndices.clear();
        readOriginIndices.clear();
        pairFlags.clear();
        scores.clear();

        numSoftClipLeftBases.clear();
        numSoftClipRightBases.clear();
        softClipLeftBases.clear();
        softClipRightBases.clear();
        softClipLeftQualityScores.clear();
        softClipRightQualityScores.clear();

        linkOffsetOptimization.clear();
        linkOffsetOptimizationIndex = 0;
        queryIndex2EntryIndices.clear();
        queryIndex2FragmentIndices.clear();
        queryIndexToPositionList.clear();

        varPositionDeltaMods.clear();
    }

    private final LinkInfo pairLinks = new LinkInfo(this, "pairs");
    private final LinkInfo forwardSpliceLinks = new LinkInfo(this, "forward-splice");
    private final LinkInfo backwardSpliceLinks = new LinkInfo(this, "backward-splice");

    /**
     * An empty sequence variation.
     */
    private static final Alignments.SequenceVariation EMPTY_SEQ_VAR = Alignments.SequenceVariation.newBuilder().build();
    private Alignments.AlignmentEntry previousPartial;
    private int countAggregatedWithMultiplicity;

    private Alignments.AlignmentEntry transform(final int index, int indexInReducedCollection, final Alignments.AlignmentEntry source) {
        final Alignments.AlignmentEntry.Builder result = Alignments.AlignmentEntry.newBuilder(source);
        final int position = source.getPosition();
        final int targetIndex = source.getTargetIndex();

        // clear the strings we collected earlier:
        result.clearSoftClippedBasesLeft();
        result.clearSoftClippedBasesRight();
        result.clearSoftClippedQualityLeft();
        result.clearSoftClippedQualityRight();

        if (index > 0 && targetIndex == previousTargetIndex) {
            result.clearPosition();
            result.clearTargetIndex();

            deltaPositions.add(position - previousPosition);
            deltaTargetIndices.add(targetIndex - previousTargetIndex);
            //       System.out.printf("entry delta position: %d delta-target=%d %n",position - previousPosition, targetIndex - previousTargetIndex);
        }

        final int queryIndex = source.getQueryIndex();

        queryIndices.add(queryIndex);
        previousPosition = position;
        previousTargetIndex = targetIndex;

        if (debug(1) && source.hasQueryLength()) {
            writtenBases += source.hasQueryLength() ? source.getQueryLength() : source.getQueryAlignedLength();
        }
        result.clearQueryIndex();

        recordVariationQualitiesAndClear(source, result, result.getSequenceVariationsList());

        final boolean entryMatchingReverseStrand = source.getMatchingReverseStrand();
        Alignments.RelatedAlignmentEntry link = pairLinks.code(source.hasPairAlignmentLink(), source,
                source.getPairAlignmentLink());
        if (link == null) {
            result.clearPairAlignmentLink();
        } else {
            result.setPairAlignmentLink(link);
        }

        link = forwardSpliceLinks.code(source.hasSplicedForwardAlignmentLink(), source, source.getSplicedForwardAlignmentLink());
        if (link == null) {
            result.clearSplicedForwardAlignmentLink();
        } else {
            result.setSplicedForwardAlignmentLink(link);
        }

        link = backwardSpliceLinks.code(source.hasSplicedBackwardAlignmentLink(), source, source.getSplicedBackwardAlignmentLink());
        if (link == null) {
            result.clearSplicedBackwardAlignmentLink();
        } else {
            result.setSplicedBackwardAlignmentLink(link);
        }
        if (source.hasReadQualityScores()) {

            final ByteString quals = source.getReadQualityScores();
            final int size = quals.size();
            numReadQualityScores.add(size);
            for (int i = 0; i < size; i++) {
                allReadQualityScores.add(quals.byteAt(i));
                qualScoreIndex++;
            }
            result.clearReadQualityScores();
        } else {
            numReadQualityScores.add(0);
        }

        if (source.hasInsertSize()) {
            final int readPos = source.getPosition();

            final int matePos = source.getPairAlignmentLink().getPosition();
            final int length = source.getTargetAlignedLength();
            final int pos1 = source.getMatchingReverseStrand() ? length + readPos : readPos + 1;
            final int pos2 = EntryFlagHelper.isMateReverseStrand(source) ? length + matePos : matePos + 1;
            final int insertSize = source.getInsertSize();
            int insertSizeDiff = Fast.int2nat(pos2 - pos1 - insertSize);
            // reverse:  insertSize= (pos2-pos1) - isd
            /*    if (insertSize != 0) {
       System.out.printf("insertSize %d length= %d %c %c readPos %d matePos %d  pos1 %d pos2 %d  pos2-pos1 %d matePos-readPos  %d  insertSizeDiff %d %n",
                 source.getInsertSize(), length,
                 source.getMatchingReverseStrand() ? '+' : '-',
                 EntryFlagHelper.isMateReverseStrand(source) ? '+' : '-',
                 readPos, matePos, pos1, pos2, pos2 - pos1, matePos - readPos, insertSizeDiff);

      }      */

            insertSizes.add(insertSizeDiff);

        } else {
            insertSizes.add(MISSING_VALUE);
        }
        result.clearInsertSize();
        // Fields above this line were removed before comparing source to the previous template.
        final Alignments.AlignmentEntry partial = result.clone().build();

        if (previousPartial != null && indexInReducedCollection >= 1 && fastEqualsEntry(previousPartial, partial)) {
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
        queryAlignedLengths.add(source.hasQueryAlignedLength() ?
                modelQueryAlignedLength(source.getQueryAlignedLength(), source.getTargetAlignedLength()) :
                MISSING_VALUE);
        targetAlignedLengths.add(source.hasTargetAlignedLength() ? source.getTargetAlignedLength() : MISSING_VALUE);
        fragmentIndices.add(source.hasFragmentIndex() ? source.getFragmentIndex() : MISSING_VALUE);
        variationCount.add(source.getSequenceVariationsCount());
        queryPositions.add(source.hasQueryPosition() ? source.getQueryPosition() : MISSING_VALUE);
        sampleIndices.add(source.hasSampleIndex() ? source.getSampleIndex() : MISSING_VALUE);
        readOriginIndices.add(source.hasReadOriginIndex() && storeReadOrigins ? source.getReadOriginIndex() : MISSING_VALUE);
        pairFlags.add(source.hasPairFlags() ? reduceSamFlags(source) : MISSING_VALUE);
        scores.add(source.hasScore() ? Float.floatToIntBits(source.getScore()) : MISSING_VALUE);


        result.clearQueryLength();
        result.clearMappingQuality();
        result.clearMatchingReverseStrand();
        result.clearMultiplicity();
        result.clearNumberOfIndels();
        result.clearNumberOfMismatches();
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
            assert seqVar.getPosition() >= 0 : String.format("The following entry had a sequence variation with a negative position. This is not allowed since seqVar.positions must be >=0. %s ", source.toString());

            encodeVar(source.getMatchingReverseStrand(), source.getQueryLength(), seqVar);
            Alignments.SequenceVariation.Builder varBuilder = Alignments.SequenceVariation.newBuilder(seqVar);
            varBuilder.clearPosition();
            varBuilder.clearFrom();
            varBuilder.clearTo();
            varBuilder.clearToQuality();
            varBuilder.clearReadIndex();
            if (!isEmpty(varBuilder.build())) {
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
        //       System.out.println(alignmentEntry);
        return alignmentEntry;
    }

    // mask the strand bit from pair flag (we store it separately anyway)
    private int reduceSamFlags(Alignments.AlignmentEntry source) {
        return source.getPairFlags() & (~16);
    }

    // reconstitute the sam Flag with strand information:
    private int restoreSamFlags(int samFlag, boolean matchesReverseStrand) {
        return samFlag | (matchesReverseStrand ? 16 : 0);
    }

    private boolean isEmpty(final Alignments.SequenceVariation varBuilder) {
        return useTemplateBasedCompression && fastEqualsInternal(varBuilder);
    }

    final ByteArrayOutputStream byteBufferSVO1 = new ByteArrayOutputStream();
    final ByteArrayOutputStream byteBufferSVO2 = new ByteArrayOutputStream();
    private byte[] EMPTY_SEQ_VAR_SERIALIZED;

    private boolean fastEqualsInternal(Alignments.SequenceVariation o2) {
        // The protobuf message.equals method is a performance  bottleneck when performing template compression. See
        //http://www.mail-archive.com/protobuf@googlegroups.com/msg02534.html
        // This method  first serializes the two messages to bytes, then compare the byte arrays. This
        // happens to be much faster than relying on the protobuff equals method (which relies on reflection).

        byteBufferSVO2.reset();
        try {
            if (EMPTY_SEQ_VAR_SERIALIZED == null) {
                byteBufferSVO1.reset();
                EMPTY_SEQ_VAR.writeTo(byteBufferSVO1);
                EMPTY_SEQ_VAR_SERIALIZED = byteBufferSVO1.toByteArray();
            }
            o2.writeTo(byteBufferSVO2);
        } catch (IOException e) {
            LOG.error("Error serializing in fastEqualsInternal", e);
            return false;
        }
        return Arrays.equals(EMPTY_SEQ_VAR_SERIALIZED, byteBufferSVO2.toByteArray());
    }

    public static int modelQueryAlignedLength(final int queryAlignedLength, final int targetAlignedLength) {
        return Fast.int2nat(queryAlignedLength - targetAlignedLength);
        // codedValue=   queryAlignedLength-targetAlignedLength
        //  queryAlignedLength=  codedValue+targetAlignedLength
    }

    public static int decodeQueryAlignedLength(int codedValue, int targetAlignedLength) {
        return Fast.nat2int(codedValue) + targetAlignedLength;
    }

    private boolean fastEqualsEntry(final Alignments.AlignmentEntry o1, final Alignments.AlignmentEntry o2) {
        return useTemplateBasedCompression && fastEqualsInternal(o1, o2);

    }

    final ByteArrayOutputStream byteBufferO1 = new ByteArrayOutputStream();
    final ByteArrayOutputStream byteBufferO2 = new ByteArrayOutputStream();

    private boolean fastEqualsInternal(final Alignments.AlignmentEntry o1, final Alignments.AlignmentEntry o2) {
        // The protobuf message.equals method is a performance  bottleneck when performing template compression. See
        //http://www.mail-archive.com/protobuf@googlegroups.com/msg02534.html
        // This method  first serializes the two messages to bytes, then compare the byte arrays. This
        // happens to be much faster than relying on the protobuff equals method (which relies on reflection).
        byteBufferO1.reset();
        byteBufferO2.reset();
        try {
            o1.writeTo(byteBufferO1);
            o2.writeTo(byteBufferO2);
        } catch (IOException e) {
            LOG.error("Error serializing in fastEqualsInternal", e);
            return false;
        }
        final boolean equals = Arrays.equals(byteBufferO1.toByteArray(), byteBufferO2.toByteArray());
        return equals;
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
            final int numReadQualScores = numReadQualityScores.get(numReadQualScoresIndex++);
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

        anInt = queryPositions.getInt(index);
        if (anInt != MISSING_VALUE) {
            result.setQueryPosition(anInt);
        }
        final int targetAlignedLength = targetAlignedLengths.getInt(index);
        if (targetAlignedLength != MISSING_VALUE) {
            result.setTargetAlignedLength(targetAlignedLength);
        }
        anInt = queryAlignedLengths.getInt(index);
        if (anInt != MISSING_VALUE) {
            result.setQueryAlignedLength(decodeQueryAlignedLength(anInt, targetAlignedLength));
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
            result.setPairFlags(restoreSamFlags(anInt, result.getMatchingReverseStrand()));
        }
        anInt = scores.getInt(index);
        if (anInt != MISSING_VALUE) {
            result.setScore(Float.intBitsToFloat(anInt));
        }
        Alignments.RelatedAlignmentEntry link = pairLinks.decode(originalIndex, result, reduced.getPairAlignmentLink());
        if (link != null) {
            result.setPairAlignmentLink(link);
        }
        link = forwardSpliceLinks.decode(originalIndex, result, reduced.getSplicedForwardAlignmentLink());
        if (link != null) {
            result.setSplicedForwardAlignmentLink(link);
        }
        link = backwardSpliceLinks.decode(originalIndex, result, reduced.getSplicedBackwardAlignmentLink());
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
            final boolean entryMatchingReverseStrand = result.hasMatchingReverseStrand() ? result.getMatchingReverseStrand() : false;
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
        if (streamVersion >= 10 && result.hasPairAlignmentLink() && result.getPairAlignmentLink().hasOptimizedIndex()) {
            // we need to reconstruct pair position before we can recalculate insert size.
            return;
        }

        final int anInt = insertSizes.getInt(index);
        if (anInt != MISSING_VALUE) {

            final int readPos = result.getPosition();

            final Alignments.RelatedAlignmentEntry pairAlignmentLink = result.getPairAlignmentLink();
            final int matePos = pairAlignmentLink.getPosition();
            final int length = result.getTargetAlignedLength();
            final int pos1 = result.getMatchingReverseStrand() ? length + readPos : readPos + 1;
            final int pos2 = EntryFlagHelper.isMateReverseStrand(result.getPairFlags()) ? length + matePos : matePos + 1;

            final int insertSize = pos2 - pos1 - Fast.nat2int(anInt);
            //  System.out.println("insertSize="+insertSize +" anInt: "+anInt);
            // reverse:  insertSize= (pos2-pos1) - isd
            result.setInsertSize(insertSize);


        }
    }

    private void recalculateInsertSize(Alignments.AlignmentEntry.Builder entry, int index) {
        if (streamVersion < 10) {
            return;
        }
        decodeInsertSize(entry, index);
    }

    // pre-allocated arrays, size 1 to 100.
    byte[][] qualArrays = new byte[100][];

    private byte[] getQualArray(final int toQualLength) {
        return qualArrays[toQualLength];
    }

    private int varQualIndex = 0;
    private int varPositionIndex = 0;
    private int varFromToIndex = 0;


    public int getNextLinkOptimizationOffset() {
        if (!enableDomainOptimizations) {
            return MISSING_VALUE;
        }
        return linkOffsetOptimization.get(linkOffsetOptimizationIndex++);
    }

    public void setDebugLevel(int level) {
        this.debug = level;
    }
}
