package edu.cornell.med.icb.counts;

import java.io.IOException;

/**
 * Minimum contract for all implementations that support iterating over counts.
 *
 * @author Fabien Campagne
 *         Date: Jun 15, 2009
 *         Time: 6:37:21 PM
 */
public interface CountsReaderI {
    /**
     * Return the position along the sequence where the count is observed.
     *
     * @return
     */
    int getPosition();

    /**
     * Determines if the reader has data about another transition.
     *
     * @return True when a call to nextTransition() will succeed, False otherwise.
     * @throws IOException
     */
    boolean hasNextTransition() throws IOException;

    /**
     * Advance to the next transition. After this method has been called successfully, position,
     * length, deltaCount and currentCount are available through getters of this reader.
     *
     * @throws IOException
     */
    void nextTransition() throws IOException;

    /**
     * Return the count at the given position.
     *
     * @return
     */
    int getCount();

    /**
     * {@inheritDoc}
     */
    void close() throws IOException;

    /**
     * Advance up to or past the specified position. The reader is advanced until the position returned by getPosition()
     * is at least equal, or greater to the specified position.
     *
     * @param position
     * @throws IOException
     */
    public void skipTo(int position) throws IOException;

    /**
     * The length of the region/peak where the count is observed.
     * @return
     */
    int getLength();
}
