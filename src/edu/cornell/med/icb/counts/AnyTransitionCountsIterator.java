package edu.cornell.med.icb.counts;

import java.io.Closeable;
import java.io.IOException;
import java.util.Arrays;

/**
 * Iterates through a set of count readers, returning a transition whenever one of the
 * underlying readers would report a transition at the position. A transition is never
 * triggered at positions where the underlying readers are all constant, making it
 * efficient to compare counts across readers.
 *
 * @author Fabien Campagne
 *         Date: Jun 13, 2009
 *         Time: 2:00:21 PM
 */
public class AnyTransitionCountsIterator implements Closeable, CountsReaderI {
    final CountsReaderI[] readers;
    final int position[];
    final int countUpToPosition[];
    private boolean nextTransitionLoaded;
    private final int numReaders;


    public AnyTransitionCountsIterator(CountsReaderI... countReader) {

        readers = countReader;
        numReaders = countReader.length;
        position = new int[numReaders];
        Arrays.fill(position, -1);
        countUpToPosition = new int[numReaders];
        currentPosition = -1;

    }

    public boolean hasNextTransition() throws IOException {
        if (nextTransitionLoaded) return true;
        // advance by at least one position:
        currentPosition++;
        // load count for each reader that transitions at this position:
        int countReadersFinished = 0;
        for (int i = 0; i < numReaders; i++) {
            if (position[i] >= currentPosition) continue;

            // passed the previous transition, must advance this reader.
            if (!readers[i].hasNextTransition()) {
                // reader has no more transitions.
                countReadersFinished++;
                position[i]++;
            } else {
                // load the next transition:
                readers[i].nextTransition();
                position[i] = readers[i].getPosition();
                //   currentPosition=currentPosition();
            }
        }

        // determine the new current position, since we may have advanced all the readers in a big jump:

        currentPosition = currentPosition();
        for (int i = 0; i < numReaders; i++) {
            // get the next count if we have reached the position of the delegate reader.
            if (position[i] == currentPosition) {
                countUpToPosition[i] = readers[i].getCount();
            }
        }
        if (countReadersFinished == numReaders) {
            // no more transitions in any reader.
            return false;
        }

        nextTransitionLoaded = true;
        return true;
    }

    /*
   Return current position and update currentLength as the minimum length of any transition at the position
   returned.
    */
    private int currentPosition() {
        int pos = Integer.MAX_VALUE;
        currentLength = Integer.MAX_VALUE;

        for (int i = 0; i < numReaders; i++) {

            pos = Math.min(pos, position[i]);
            if (pos == position[i]) {
                currentLength = Math.min(readers[i].getLength(), currentLength);
            }
        }
        return pos;
    }

    int currentPosition;
    private int currentLength;

    /**
     * Advance to the next transition. After this method has been called successfully, position,
     * length, deltaCount and currentCount are available through getters of this reader.
     *
     * @throws IOException
     */
    public void nextTransition() throws IOException {
        if (hasNextTransition()) {
            nextTransitionLoaded = false;

        } else {
            throw new IllegalStateException("next cannot be called when hasNext would return false.");
        }
    }

    /**
     * Return the sum of counts over the readers that have non zero counts at the current position.
     */
    public int getCount() {
        int count = 0;
        for (CountsReaderI reader : getReaders()) {
            if (reader.getPosition() == getPosition()) count += reader.getCount();
        }
        return count;
    }

    public final CountsReaderI[] getReaders() {
        return readers;
    }

    public final int getCount(int readerIndex) {
        return countUpToPosition[readerIndex];
    }

    public void close() throws IOException {
        for (CountsReaderI reader : readers) {
            reader.close();
        }
    }

    public void skipTo(int position) throws IOException {
        throw new UnsupportedOperationException("skipTo is not currently supported by this implementation.");
    }

    public int getLength() {
        return currentLength;
    }

    public int getPosition() {
        return currentPosition;
    }
}
