package edu.cornell.med.icb.counts;

import java.io.IOException;

/**
 * A facade over a countsReader that returns positions summed with an offset.
 * The position of count transitions returned is the position of the transition
 * returned by the delegate reader, plus the offset.
 *
 * @author Fabien Campagne
 *         Date: Jun 15, 2009
 *         Time: 6:41:52 PM
 */
public class OffsetCountsReader implements Cloneable, CountsReaderI {
    CountsReader delegate;
    private int offset;

    /**
     * {@inheritDoc}
     */
    public int getCount() {
        return delegate.getCount();
    }

    /**
     * {@inheritDoc}
     */
    public void close() throws IOException {
        delegate.close();
    }

    /**
     * {@inheritDoc}
     */
    public void skipTo(int position) throws IOException {
        delegate.skipTo(position - offset);
    }

    public int getLength() {
        return delegate.getLength();
    }

    /**
     * {@inheritDoc}
     */
    public int getPosition() {
        return delegate.getPosition() + offset;
    }

    /**
     * {@inheritDoc}
     */
    public boolean hasNextTransition() throws IOException {
        return delegate.hasNextTransition();
    }

    /**
     * {@inheritDoc}
     */
    public void nextTransition() throws IOException {
        delegate.nextTransition();
    }

    public OffsetCountsReader(CountsReader delegate, int offset) {
        this.delegate = delegate;
        this.offset = offset;
    }

    /**
     * Change the offset.
     *
     * @param offset New value that is being added to count transition positions.
     */
    public void setOffset(int offset) {
        this.offset = offset;
    }
}
