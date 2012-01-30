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

package edu.cornell.med.icb.goby.modes;

import edu.cornell.med.icb.goby.algorithmic.data.GroupComparison;
import edu.cornell.med.icb.goby.algorithmic.data.MethylCountInfo;
import edu.cornell.med.icb.goby.alignments.DiscoverVariantIterateSortedAlignments;
import edu.cornell.med.icb.goby.alignments.DiscoverVariantPositionData;
import edu.cornell.med.icb.goby.alignments.ReadIndexStats;
import edu.cornell.med.icb.goby.alignments.SampleCountInfo;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import edu.cornell.med.icb.goby.stats.VCFWriter;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.PrintWriter;
import java.util.ArrayList;

/**
 * A Variant Call Format to estimate average methylation rates over parts of the genome that overlap annotations or other
 * regions.
 * Skeleton for Nyasha to wire her implementation in.
 *
 * @author Fabien Campagne
 *         Date: Jan 28 2012
 *         Time: 1:50:13 PM
 */
public class MethylationRegionsOutputFormat implements SequenceVariationOutputFormat {

    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(MethylationRegionsOutputFormat.class);

    VCFWriter statWriter;
    private int numberOfSamples;
    private int numberOfGroups;
    private String[] samples;
    private String[] groups;
    /**
     * Maps sampleIndex to group index.
     */
    int[] readerIndexToGroupIndex;

    private ArrayList<GroupComparison> groupComparisons = new ArrayList<GroupComparison>();

    private RandomAccessSequenceInterface genome;
    /**
       Where number of methylated and un-methylated bases are stored for each sample and
       group under investigation:
     */
    private MethylCountInfo mci;


    public void defineColumns(final PrintWriter writer, final DiscoverSequenceVariantsMode mode) {

        groups = mode.getGroups();
        samples = mode.getSamples();

        readerIndexToGroupIndex = mode.getReaderIndexToGroupIndex();
        final ObjectArrayList<ReadIndexStats> readIndexStats = mode.getReadIndexStats();
        this.statWriter = new VCFWriter(writer);
        groupComparisons = mode.getGroupComparisons();

        if (groups.length < 1) {
            System.err.println("Methylation format requires at least one group.");
            System.exit(1);
        }


    }

    @Override
    public void allocateStorage(final int numberOfSamples, final int numberOfGroups) {
        this.numberOfGroups = numberOfGroups;
        this.numberOfSamples = numberOfSamples;
        mci = new MethylCountInfo(numberOfSamples, numberOfGroups);

    }
    private String chromosome;
    private int referenceIndex;
    private int position;

    @Override
    public void writeRecord(final DiscoverVariantIterateSortedAlignments iterator,
                            final SampleCountInfo[] sampleCounts,
                            final int referenceIndex,
                            int position,
                            final DiscoverVariantPositionData list,
                            final int groupIndexA,
                            final int groupIndexB) {


        position = position + 1;

        final char refBase = sampleCounts[0].referenceBase;
        if (refBase != 'C' && refBase != 'G') {
            return;
        }
        MethylationRateVCFOutputFormat.fillMethylationCountArrays(sampleCounts, list, position, refBase, mci, readerIndexToGroupIndex);
        //TODO: hook the averaging writer in here (C and Cm are in mci):
        this.position=position;
        if (referenceIndex!=this.referenceIndex) {
            this.referenceIndex=referenceIndex;
            chromosome=genome.getReferenceName(referenceIndex);
        }
        statWriter.writeRecord();
    }

    public int getPosition() {
        return position;
    }

    public String getChromosome() {

        return chromosome;
    }

    @Override
    public void close() {
        //TODO: hook the averaging writer in here:
        statWriter.close();
    }

    @Override
    public void setGenome(RandomAccessSequenceInterface genome) {
        this.genome = genome;
    }


}