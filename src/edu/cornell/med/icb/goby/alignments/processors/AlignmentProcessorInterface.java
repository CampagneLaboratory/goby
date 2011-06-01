/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.alignments.processors;

import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import edu.cornell.med.icb.goby.alignments.Alignments;

import java.io.IOException;
import java.io.PrintWriter;

/**
 * @author Fabien Campagne
 *         Date: May 12, 2011
 *         Time: 9:16:42 AM
 */
public interface AlignmentProcessorInterface {
    /**
     * Provides the next processed entry. The processor had the opportunity to modify the entry before returning it.
     *
     * @param targetIndex index of the reference sequence, similar to a call to AlignmentReaderI.skipTo()
     * @param position    position on the reference sequence, similar to a call to AlignmentReaderI.skipTo()
     * @return The next available alignment entry.
     * @throws IOException If an error occurred reading the alignment.
     */
    Alignments.AlignmentEntry nextRealignedEntry(int targetIndex, int position) throws IOException;

    /**
     * Provide this processor with a compressed genome. The genome must match the reference against which the alignment was
     * produced.
     *
     * @param genome The genome corresponding to the alignment
     */
    void setGenome(RandomAccessSequenceInterface genome);

    /**
     * Returns the number of entries actually modified by the processor. This is garanteed to be less or equal to
     * getProcessedCount().
     *
     * @return the number of entries modified by this processor.
     */
    int getModifiedCount();

    /**
     * Returns the number of entries processed by the processor. This is exactly the number of times  nextRealignedEntry
     * was called.
     *
     * @return the number of entries modified by this processor.
     */
    int getProcessedCount();
}
