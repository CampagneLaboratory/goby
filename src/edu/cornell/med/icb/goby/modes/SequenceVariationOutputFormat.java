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

package edu.cornell.med.icb.goby.modes;

import edu.cornell.med.icb.goby.alignments.ContigHelper;
import edu.cornell.med.icb.goby.alignments.DiscoverVariantIterateSortedAlignments;
import edu.cornell.med.icb.goby.alignments.DiscoverVariantPositionData;
import edu.cornell.med.icb.goby.alignments.SampleCountInfo;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import edu.cornell.med.icb.goby.util.OutputInfo;

import java.io.File;

/**
 * Interface to support different kinds of outputs for discover sequence variants mode.
 * @author Fabien Campagne
 *         Date: Mar 20, 2011
 *         Time: 11:43:56 AM
 */
public interface SequenceVariationOutputFormat {
    /**
     * Define columns needed by this file format in the statsWriter.
     * @param statsWriter where output will be written.
     * @param mode The discover mode provides arguments that the format may need to define columns.
     */
    public void defineColumns(OutputInfo statsWriter, DiscoverSequenceVariantsMode mode);

    /**
     * Allocate storate needed to estimate the statistics the format outputs.
     * @param numberOfSamples Number of samples in the comparison
     * @param numberOfGroups Number of groups being compared.
     */
    public void allocateStorage(int numberOfSamples, int numberOfGroups);

    /**
     * Write statistics to the statsWriter for each position that should be reported.
     * @param iterator
     * @param sampleCounts
     * @param referenceIndex
     * @param position
     * @param list
     * @param groupIndexA
     * @param groupIndexB
     */
    public void writeRecord(DiscoverVariantIterateSortedAlignments iterator,
                            SampleCountInfo[] sampleCounts, int referenceIndex,
                            int position,
                            DiscoverVariantPositionData list,
                            int groupIndexA, int groupIndexB);


    void close();

    /**
     * Optionally give access to the genome to the file format writer.
     * @param genome
     */
    void setGenome(RandomAccessSequenceInterface genome);

    void setGenomeReferenceIndex(int index);
}
