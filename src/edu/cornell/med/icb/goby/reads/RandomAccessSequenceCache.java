/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 *  This file is part of the Goby IO API.
 *
 *     The Goby IO API is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     The Goby IO API is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with the Goby IO API.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.cornell.med.icb.goby.reads;

import com.google.protobuf.ByteString;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.parsers.ReaderFastaParser;
import it.unimi.dsi.bits.BitVector;
import it.unimi.dsi.bits.LongArrayBitVector;
import it.unimi.dsi.fastutil.ints.Int2ObjectArrayMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.io.BinIO;
import it.unimi.dsi.fastutil.objects.Object2IntArrayMap;
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.io.OutputBitStream;
import it.unimi.dsi.lang.MutableString;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.Random;
import java.util.zip.GZIPInputStream;

/**
 * Load a genome into memory and provide random access to individual bases. Supports DNA (ACTG) and
 * other bases (encoded as N).
 *
 * @author Fabien Campagne
 *         Date: May 19, 2009
 *         Time: 3:08:21 PM
 */
public class RandomAccessSequenceCache implements RandomAccessSequenceInterface {
    private ObjectArrayList<LongArrayBitVector> referenceIgnoreLists;
    private Object2IntMap<String> referenceNameMap;
    private Int2ObjectMap<String> indexToNameMap;
    private ObjectArrayList<byte[]> compressedData;
    private IntList sizes;
    private static final Logger LOG = Logger.getLogger(RandomAccessSequenceCache.class);
    private String basename;
    /**
     * All queries must be contained within   minRefIndex and maxRefIndex (inclusive).
     */
    private int maxRefIndex;
    private int minRefIndex;

    public RandomAccessSequenceCache() {
        super();
        compressedData = new ObjectArrayList<byte[]>();
        referenceIgnoreLists = new ObjectArrayList<LongArrayBitVector>();
        referenceNameMap = new Object2IntOpenHashMap<String>();
        referenceNameMap.defaultReturnValue(-1);
        indexToNameMap = new Int2ObjectArrayMap<String>();
        sizes = new IntArrayList();
    }

    /**
     * Load a fasta file into this cache.
     *
     * @param reader
     * @throws IOException
     */
    public void loadFasta(final Reader reader) throws IOException {
        final ReaderFastaParser parser = new ReaderFastaParser(reader);
        final MutableString description = new MutableString();
        int refIndex = 0;
        minRefIndex=Integer.MAX_VALUE;
        maxRefIndex=Integer.MIN_VALUE;
        while (parser.hasNextSequence()) {
            int position = 0;
            parser.nextSequence(description);
            final String referenceName = description.toString().split(" ")[0];
            final int initialCapacity = 40000000;
            final ByteArrayOutputStream byteArrayOutputStream = new ByteArrayOutputStream(initialCapacity);
            final OutputBitStream compressed = new OutputBitStream(byteArrayOutputStream);

            final LongArrayBitVector referenceIgnoreList = LongArrayBitVector.ofLength(0);

            //  referenceIgnoreList.defaultReturnValue(false);
            final Reader baseReader = parser.getBaseReader();
            int c;
            while ((c = baseReader.read()) != -1) {
                encode(c, compressed, referenceIgnoreList);
                position++;
            }
            baseReader.close();
            compressed.flush();
            compressed.close();
            referenceIgnoreLists.add(referenceIgnoreList);
            referenceNameMap.put(referenceName, refIndex);
            indexToNameMap.put(refIndex, referenceName);
            updateSliceIndices(refIndex);

            refIndex++;

            final byte[] bytes = byteArrayOutputStream.toByteArray();
            LOG.debug("size of last sequence " + description + ", in bytes: " + bytes.length);

            compressedData.add(bytes);
            sizes.add(position + 1);
        }

    }

    private void updateSliceIndices(final int refIndex) {

         this.minRefIndex=Math.min(refIndex,minRefIndex);
         this.maxRefIndex=Math.max(refIndex, maxRefIndex);

    }

    /**
     * Load a fasta file into this cache.
     *
     * @param compactInput inputstream for a compact reads file.
     * @throws IOException
     */
    public void loadCompact(final InputStream compactInput) throws IOException {
        final ReadsReader parser = new ReadsReader(compactInput);
        final MutableString description = new MutableString();
        int refIndex = 0;
        Reads.ReadEntry entry;
         minRefIndex=Integer.MAX_VALUE;
        maxRefIndex=Integer.MIN_VALUE;
        while (parser.hasNext()) {
            entry = parser.next();

            final String referenceName = entry.getReadIdentifier();
            final int initialCapacity = 40000000;
            final ByteArrayOutputStream byteArrayOutputStream = new ByteArrayOutputStream(initialCapacity);
            final OutputBitStream compressed = new OutputBitStream(byteArrayOutputStream);

            final LongArrayBitVector referenceIgnoreList = LongArrayBitVector.ofLength(0);
            final ByteString seq = entry.getSequence();
            for (int position = 0; position < seq.size(); ++position) {
                final char c = (char) seq.byteAt(position);
                encode(c, compressed, referenceIgnoreList);

            }


            compressed.flush();
            compressed.close();
            referenceIgnoreLists.add(referenceIgnoreList);
            referenceNameMap.put(referenceName, refIndex);
            indexToNameMap.put(refIndex, referenceName);
            updateSliceIndices(refIndex);

            refIndex++;

            final byte[] bytes = byteArrayOutputStream.toByteArray();
            LOG.debug("size of last sequence " + description + ", in bytes: " + bytes.length);

            compressedData.add(bytes);
            sizes.add(seq.size() + 1);
        }
    }

    public void save(final String basename) throws IOException {
        BinIO.storeObject(sizes, basename + ".sizes");
        BinIO.storeObject(compressedData, basename + ".bases");
        BinIO.storeObject(referenceIgnoreLists, basename + ".ignore");
        BinIO.storeObject(referenceNameMap, basename + ".names");
    }

    @SuppressWarnings("unchecked")
    public void load(final String basename) throws IOException, ClassNotFoundException {
        sizes = (IntList) BinIO.loadObject(basename + ".sizes");
        compressedData = (ObjectArrayList<byte[]>) BinIO.loadObject(basename + ".bases");
        referenceIgnoreLists = (ObjectArrayList<LongArrayBitVector>) BinIO.loadObject(basename + ".ignore");
        referenceNameMap = (Object2IntMap<String>) BinIO.loadObject(basename + ".names");

        for (final String name : referenceNameMap.keySet()) {
            indexToNameMap.put(referenceNameMap.get(name), name);
        }
    }

    /**
     * Load a slice of this genome, contained between minReferenceId position 0 and max referenceId position 0.
     * The special values "min" and "max" can be used for  minRefId and maxRefId respectively to retrieve the entire
     * genome.
     *
     * @param basename
     * @param minRefId
     * @param maxRefId
     * @throws IOException
     * @throws ClassNotFoundException
     */
    public void load(final String basename, final String minRefId, final String maxRefId) throws IOException, ClassNotFoundException {

        load(basename);
        minRefIndex = "min".equals(minRefId) ? -1 : referenceNameMap.getInt(minRefId);
        maxRefIndex = "max".equals(maxRefId) ? referenceNameMap.size() : referenceNameMap.getInt(maxRefId);

        // remove sequences we won't need
        for (int i = 0; i < minRefIndex; i++) {
            compressedData.rem(i);
            referenceIgnoreLists.rem(i);
        }
        for (int i = maxRefIndex + 1; i < referenceNameMap.size(); i++) {
            compressedData.rem(i);
            referenceIgnoreLists.rem(i);
        }

    }

    public boolean canLoad(final String basename) {
        final String[] extensions = {
                ".sizes", ".bases", ".ignore", ".names"
        };
        for (final String extension : extensions) {
            final File part = new File(basename + extension);
            if (!part.exists()) {
                return false;
            }
        }
        return true;

    }

    public final char get(final String referenceName, final int position) {
        final int referenceIndex = getReferenceIndex(referenceName);
        return get(referenceIndex, position);
    }

    public final int getRange(final String referenceName, final int position, final int length) {
        final int referenceIndex = getReferenceIndex(referenceName);
        return getRange(referenceIndex, position, length);
    }

    @Override
    public void getRange(final int referenceIndex, final int position, final int length, final MutableString bases) {
        bases.setLength(0);
        for (int i = position; i < position + length; i++) {
            bases.append(get(referenceIndex, i));
        }
    }

    final LongArrayBitVector bits = LongArrayBitVector.getInstance();

    public int getRange(final int referenceIndex, final int position, final int length) {
        assert referenceIndex >= minRefIndex && referenceIndex <= maxRefIndex :
                String.format("referenceindex %d out of genome slice [%d-%d].", referenceIndex,
                        minRefIndex, maxRefIndex);
        final int maxSize = sizes.getInt(referenceIndex);
        assert position + length < maxSize : "position must be less than size of the reference sequence (" + maxSize + ")";


        assert length < 15 : "length must be less than 15";
        bits.clear();
        final byte[] bytes = compressedData.get(referenceIndex);
        final LongArrayBitVector ignoreList = referenceIgnoreLists.get(referenceIndex);

        for (int i = 0; i < length; i++) {
            final int offset = (position + i) * 2;
            final byte b = bytes[offset / 8];
            if (ignoreList.get(position)) {
                // a range that contain 'N' at any position is represented by -1.
                return -1;
            }
            final int c = b >> (6 - (offset % 8)) & 0x3; // two right-most bits are left, which encode
            switch (c) {
                case 2 * 1 + 1 * 1:
                    bits.add(0);
                    bits.add(0);
                    break;
                case 2 * 0 + 1 * 1:
                    bits.add(1);
                    bits.add(0);
                    break;
                case 2 * 1 + 1 * 0:
                    bits.add(0);
                    bits.add(1);
                    break;
                case 2 * 0 + 1 * 0:
                    bits.add(1);
                    bits.add(1);
                    break;
                default:
                    throw new InternalError("Should never happen");


            }
        }
        return (int) bits.bits()[0];
    }

    /**
     * Return the index of the reference sequence identified by name, or -1 if the sequence name
     * is not in the cache.
     *
     * @param referenceName The name of the sequence to get the index for
     * @return The index for the specified reference
     */

    public final int getReferenceIndex(final String referenceName) {
        return referenceNameMap.getInt(referenceName);
    }

    /**
     * Return the reference name corresponding to this index.
     *
     * @param index for the specified reference
     * @return referenceName The name of the sequence to get the index for
     */

    public final String getReferenceName(final int index) {
        return indexToNameMap.get(index);
    }

    @Override
    public int size() {
        return referenceNameMap.size();
    }

    public final char get(final int referenceIndex, final int position) {
        assert referenceIndex >= minRefIndex && referenceIndex <= maxRefIndex :
                String.format("referenceindex %d out of genome slice [%d-%d].", referenceIndex,
                        minRefIndex, maxRefIndex);

        final int maxSize = sizes.getInt(referenceIndex);
        final LongArrayBitVector ignoreList = referenceIgnoreLists.get(referenceIndex);

        assert position < maxSize : "position must be less than size of the reference sequence (" + maxSize + ")";
        assert position < ignoreList.length() : " position must be smaller than ignore list size.";

        if (!ignoreList.get(position)) {
            return decode(compressedData.get(referenceIndex), position,
                    maxSize);
        } else {
            return 'N';
        }
    }


    public int getLength(int targetIndex) {
        return sizes.getInt(targetIndex);
    }

    private char decode(final byte[] bytes, final int position, final int maxSize) {
        assert position < maxSize : "position must be less than size of the reference sequence (" + maxSize + ")";

        final int offset = position * 2;
        final byte b = bytes[offset / 8];
        final int c = b >> (6 - (offset % 8)) & 0x3; // two right-most bits are left, which encode

        switch (c) {
            case 2 * 1 + 1 * 1:
                return 'A';
            case 2 * 0 + 1 * 1:
                return 'C';
            case 2 * 1 + 1 * 0:
                return 'T';
            case 2 * 0 + 1 * 0:
                return 'G';
            default:
                throw new InternalError("This should never happen");

        }
    }

    public static void main(final String[] args) throws IOException, ClassNotFoundException {
        final RandomAccessSequenceCache cache = new RandomAccessSequenceCache();
        final String basename = "compressed-genome-cache";
        if (cache.canLoad(basename)) {
            cache.load(basename);
        } else {
            System.out.println("Loading");
            cache.loadFasta(new InputStreamReader(new GZIPInputStream(new FileInputStream("/Users/fac2003/IdeaProjects/data/Homo_sapiens.NCBI36.54.dna.toplevel.fa.gz"))));

        }
        System.out.println("Done loading.");
        if (!cache.canLoad(basename)) {
            cache.save(basename);
            System.out.println("Sequence cache written to disk with basename " + basename);
        }
        System.out.println("Searching..");
        final Random random = new Random();
        for (int i = 0; i < 1000; i++) {
            for (int j = 0; j < 10000; j++) {
                final int referenceIndex = random.nextInt(cache.numberOfSequences() - 1);
                final int position = random.nextInt(cache.size(referenceIndex) - 1);
                cache.get(referenceIndex, position);
            }
        }
        System.out.println("Done searching");
    }

    private int size(final int referenceIndex) {

        final LongArrayBitVector ignoreList = referenceIgnoreLists.get(referenceIndex);
        return Math.min(ignoreList.size(), sizes.get(referenceIndex));
    }

    public int numberOfSequences() {
        return sizes.size();
    }

    private void encode(final int c,
                        final OutputBitStream compressed,
                        final BitVector ignoreList) throws IOException {
        switch (c) {
            case 'A':
                compressed.writeBit(1);
                compressed.writeBit(1);
                ignoreList.add(false);
                break;
            case 'C':
                compressed.writeBit(0);
                compressed.writeBit(1);
                ignoreList.add(false);
                break;
            case 'T':
                compressed.writeBit(1);
                compressed.writeBit(0);
                ignoreList.add(false);
                break;
            case 'G':
                compressed.writeBit(0);
                compressed.writeBit(0);
                ignoreList.add(false);
                break;
            default:
                compressed.writeBit(0);
                compressed.writeBit(0);
                ignoreList.add(true);
                break;
        }
    }


    public int getSequenceSize(final int referenceIndex) {
        return size(referenceIndex);
    }

    public String getBasename() {
        return basename;
    }

    public void setBasename(String basename) {
        this.basename = basename;
    }
}
