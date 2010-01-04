/*
 * Copyright (C) 2009 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 * WEILL MEDICAL COLLEGE OF CORNELL UNIVERSITY MAKES NO REPRESENTATIONS
 * ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY PURPOSE. IT IS PROVIDED
 * "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY. THE WEILL MEDICAL COLLEGE
 * OF CORNELL UNIVERSITY SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY
 * THE USERS OF THIS SOFTWARE.
 */

package edu.cornell.med.icb.goby.reads;

import edu.cornell.med.icb.parsers.ReaderFastaParser;
import it.unimi.dsi.bits.BitVector;
import it.unimi.dsi.bits.LongArrayBitVector;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.io.BinIO;
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.io.OutputBitStream;
import it.unimi.dsi.lang.MutableString;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
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
public class RandomAccessSequenceCache {
    private ObjectArrayList<LongArrayBitVector> referenceIgnoreLists;
    private Object2IntMap<String> referenceNameMap;
    private ObjectArrayList<byte[]> compressedData;
    private IntList sizes;

    public RandomAccessSequenceCache() {
        super();
        compressedData = new ObjectArrayList<byte[]>();
        referenceIgnoreLists = new ObjectArrayList<LongArrayBitVector>();
        referenceNameMap = new Object2IntOpenHashMap<String>();
        referenceNameMap.defaultReturnValue(-1);
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

            refIndex++;

            final byte[] bytes = byteArrayOutputStream.toByteArray();
            System.out.println("size of last sequence " + description + ", in bytes: " + bytes.length);
            compressedData.add(bytes);
            sizes.add(position + 1);
        }
    }

    public void save(final String basename) throws IOException {
        BinIO.storeObject(sizes, basename + ".sizes");
        BinIO.storeObject(compressedData, basename + ".bases");
        BinIO.storeObject(referenceIgnoreLists, basename + ".ignore");
        BinIO.storeObject(referenceNameMap, basename + ".names");
    }

    public void load(final String basename) throws IOException, ClassNotFoundException {

        sizes = (IntList) BinIO.loadObject(basename + ".sizes");
        compressedData = (ObjectArrayList<byte[]>) BinIO.loadObject(basename + ".bases");
        referenceIgnoreLists = (ObjectArrayList<LongArrayBitVector>) BinIO.loadObject(basename + ".ignore");
        referenceNameMap = (Object2IntMap<String>) BinIO.loadObject(basename + ".names");
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

    /**
     * Return the index of the reference sequence identified by name, or -1 if the sequence name
     * is not in the cache.
     *
     * @param referenceName
     * @return
     */
    public final int getReferenceIndex(final String referenceName) {
        return referenceNameMap.getInt(referenceName);
    }

    public final char get(final int referenceIndex, final int position) {
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
        return sizes.get(referenceIndex);
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


    public int getSequenceSize(int referenceIndex) {
        return size(referenceIndex);
    }
}
