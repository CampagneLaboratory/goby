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

import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.alignments.ConcatAlignmentReader;
import edu.cornell.med.icb.goby.alignments.ConcatSortedAlignmentReader;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;

import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: 5/31/11
 *         Time: 5:35 PM
 */
public class DummyProcessorUnsorted implements AlignmentProcessorInterface {

    private ConcatAlignmentReader reader;


    /**
     * Set the input reader.
     *
     * @param reader Input reader.
     */
    public DummyProcessorUnsorted(ConcatAlignmentReader reader) {
        this.reader = reader;

    }

    /**
     * Return the next available entry in the input reader (calls next()), at targetIndex and position. No processing is done.
     *
     * @param targetIndex
     * @param position
     * @return An entry, or null when hasNext() of the delegate is false.
     * @throws java.io.IOException
     */
    public Alignments.AlignmentEntry nextRealignedEntry(int targetIndex, int position) throws IOException {
        ++processedCount;
        if (reader.hasNext()) {
            return reader.next();
        }   else {
            return null;
        }
    }

    public void setGenome(RandomAccessSequenceInterface genome) {

    }

    private int processedCount;

    @Override
    public int getModifiedCount() {
        return 0;
    }

    @Override
    public int getProcessedCount() {
        return processedCount;
    }
}
