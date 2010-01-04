package edu.cornell.med.icb.goby.reads;

import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.io.InputBitStream;
import it.unimi.dsi.io.OutputBitStream;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Collections;

/**
 * A set  of read/query indices and associated multiplicity. Multiplicity is an int associated with each read. It can
 * represent anything, but is often used to stored the number of times a read is repeated in a sample file or dataset.
 * The read indices are stored with delta coding, and multiplicity values with gamma coding, keeping the file representation
 * of the set small even when there are million of reads.
 *
 * @author Fabien Campagne
 *         Date: Jun 3, 2009
 *         Time: 3:22:49 PM
 */
public class ReadSet {
    final IntSet filter = new IntOpenHashSet();
    final Int2IntMap multiplicityMap;
    private int smallestStoredMultiplicity;

    public ReadSet() {
        multiplicityMap = new Int2IntOpenHashMap();
    }


    public void clear() {
        filter.clear();
        multiplicityMap.clear();
    }

    public void add(int readIndex, int multiplicity) {
        assert readIndex >= 0 : "read indices must be positive";
        filter.add(readIndex);
        if (multiplicity >= smallestStoredMultiplicity) {

            multiplicityMap.put(readIndex, multiplicity);
        }
    }

    public boolean contains(int readIndex) {
        return filter.contains(readIndex);

    }

    /**
     * Save the filter to disk.
     *
     * @param basename
     * @param suffix
     * @throws IOException
     */
    public void save(String basename, String suffix) throws IOException {
        String filename = basename + "-" + suffix + ".filter";
        FileOutputStream stream = new FileOutputStream(filename);
        OutputBitStream out = new OutputBitStream(stream);


        // sort in increasing order:
        IntList sorted = new IntArrayList();
        sorted.addAll(filter);
        Collections.sort(sorted);

        // append to output:
        out.writeGamma(sorted.size());
        int previous = -1;
        for (int readIndex : sorted) {
            out.writeDelta(readIndex - previous);
            out.writeGamma(multiplicityMap.get(readIndex));
            previous = readIndex;

        }
        out.close();
    }

    /**
     * Load the filter from disk.
     *
     * @param basename
     * @param suffix
     * @throws IOException
     */
    public void load(String basename, String suffix) throws IOException {
        String filename = basename + "-" + suffix + ".filter";
        load(new File(filename));

    }

    /**
     * Load a read set from disk.
     * @param file File to load.
     * @throws IOException If an error occurs reading the read set file.
     */
    public void load(File file) throws IOException {
        FileInputStream stream = new FileInputStream(file);
        InputBitStream in = new InputBitStream(stream);

        filter.clear();
        int numDeltas = in.readGamma();
        int previous = -1;
        for (int i = 0; i < numDeltas; i++) {
            int delta = in.readDelta();
            int multiplicity = in.readGamma();
            final int readIndex = previous + delta;
            filter.add(readIndex);
            multiplicityMap.put(readIndex, multiplicity);
            previous = readIndex;
        }
        in.close();
    }

    public void add(int readIndex) {
        add(readIndex, 0);
    }

    public int size() {
        return filter.size();
    }

    /**
     * Returns the multiplicity of the read, or zero is the read is not in this read set.
     * @param readIndex
     * @return
     */
    public int getMultiplicity(int readIndex) {
        if (filter.contains(readIndex)) {
            return multiplicityMap.get(readIndex);
        } else return 0;
    }

    /**
     * Add all the read indices with the given multiplicity.
     *
     * @param otherReadIndices Set of read indices to add.
     * @param multiplicity
     */
    public void addAll(IntSet otherReadIndices, int multiplicity) {
        for (int readIndex : otherReadIndices) {
            add(readIndex, multiplicity);
        }
    }

    /**
     * Set the default multiplicity value
     *
     * @param value
     */
    public void smallestStoredMultiplicity(int value) {
        multiplicityMap.defaultReturnValue(value);
        smallestStoredMultiplicity = value;
    }
}
