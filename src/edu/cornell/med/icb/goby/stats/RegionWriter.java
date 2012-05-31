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

package edu.cornell.med.icb.goby.stats;

import edu.cornell.med.icb.goby.algorithmic.data.GroupComparison;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;

import java.util.ArrayList;

/**
 * Interface for all region writers. Region writers write output files about genomic regions. Some writers
 * define regions with annotations, others discover regions de-novo. Regions writers obtain data from a provider,
 * that is typically given to the interface implementation at construction time. The writeRecord method obtains
 * data for a specific site from the provider, writes output when needed, then advance the provider to the next site
 * and waits to be called again. Interactions end when the close method is called.
 * @author Fabien Campagne
 *         Date: 2/19/12
 *         Time: 12:19 PM
 */
public interface RegionWriter {
    /**
     * Obtain information from the provider and write as needed to the output.
     */
    void writeRecord();

    /**
     * Write any remaining information to the output and close it.
     */
    void close();

    /**
     * Set the group comparisons.
     *
     * @param groupComparisons
     */
    void setGroupComparisons(ArrayList<GroupComparison> groupComparisons);

    /**
     * Set the genome. Can be used by the region writer to obtain the genomic context of sites.
     *
     * @param genome
     */
    void setGenome(RandomAccessSequenceInterface genome);

    /**
     * Set the annotation filename.
     *
     * @param annotationFilename
     */
    void setAnnotationFilename(String annotationFilename);

    /**
     * Set the sample index to group index array
     *
     * @param readerIndexToGroupIndex
     */
    void setSampleIndexToGroupIndex(int[] readerIndexToGroupIndex);

    /**
     * Indicate that the writer should combine all contexts (True) or write output where contexts are
     * separately described (False).
     *
     * @param aggregateAllContexts
     */
    void setAggregateAllContexts(boolean aggregateAllContexts);


}
