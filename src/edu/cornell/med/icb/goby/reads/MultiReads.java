package edu.cornell.med.icb.goby.reads;

import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import it.unimi.dsi.fastutil.ints.IntArrayList;

/**
 * @author Fabien Campagne
 *         Date: Jun 20, 2009
 *         Time: 10:51:11 AM
 */
public class MultiReads {
    int readIndex = -1;
    int currentByteIndex = 0;
    ByteArrayList bytes = new ByteArrayList();
    private IntArrayList ends = new IntArrayList();
    private IntArrayList starts = new IntArrayList();

    public void newRead() {
        // previous read ends at previous index.
        ends.set(readIndex, currentByteIndex - 1);
        readIndex += 1;
        starts.set(readIndex, currentByteIndex);
    }

    public void addByte(byte base) {
        currentByteIndex += 1;
        bytes.add(base);

    }
}
