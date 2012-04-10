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
import edu.cornell.med.icb.goby.stats.AnnotationAveragingWriter;
import edu.cornell.med.icb.goby.stats.MethylCountProviderFromRegionsOutputFormat;
import edu.cornell.med.icb.goby.stats.RegionWriter;
import edu.cornell.med.icb.goby.util.dynoptions.DynamicOptionClient;
import edu.cornell.med.icb.goby.util.OutputInfo;
import edu.cornell.med.icb.goby.util.dynoptions.RegisterThis;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

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
public class MethylationRegionsOutputFormat implements SequenceVariationOutputFormat, MethylationFormat {

    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(MethylationRegionsOutputFormat.class);
    @RegisterThis
    public static DynamicOptionClient doc = new DynamicOptionClient(MethylationRegionsOutputFormat.class, "annotations:filename to a tab delimited annotation file:",
            "do-indel-rate:boolean, true value indicates that the indel rate should be output in the MR field:false",
            "de-novo-regions:boolean, true indicates that regions should be discovered de-novo, false indicates that regions are defined by annotations:false"
    );

    public static DynamicOptionClient doc() {
        return doc;
    }

    private String annotationFilename;


    private int numberOfSamples;
    private int numberOfGroups;
    private String[] samples;
    private String[] groups;
    /**
     * Maps sampleIndex to group index.
     */
    int[] readerIndexToGroupIndex;
    /**
     * The averaging writer that overlaps sites with annotations and writes averages:
     */
    private RegionWriter averagingWriter;
    private int minimumEventThreshold;
    private Boolean doIndels;
    private Boolean deNovoRegions;
    private int genomeReferenceIndex;

    public MethylationRegionsOutputFormat() {
        deNovoRegions = doc.getBoolean("de-novo-regions");
        assert !deNovoRegions : "de novo regions are not currently supported.";
        if (deNovoRegions == false) {
            final String annotations = doc.getString("annotations");

            if (annotations != null) {
                annotationFilename = annotations;
            } else {
                // there is no annotation file. Future, to use direct discovery of regions.
                annotationFilename = null;
            }
        }
        doIndels = doc.getBoolean("do-indel-rate");
    }

    private ArrayList<GroupComparison> groupComparisons = new ArrayList<GroupComparison>();

    private RandomAccessSequenceInterface genome;
    /**
     * Where number of methylated and un-methylated bases are stored for each sample and
     * group under investigation:
     */
    private MethylCountInfo mci;


    public String[] getGroups() {
        return groups;
    }

    public void defineColumns(final OutputInfo outputInfo, final DiscoverSequenceVariantsMode mode) {

        groups = mode.getGroups();
        samples = mode.getSamples();
        final MethylCountProviderFromRegionsOutputFormat provider = new MethylCountProviderFromRegionsOutputFormat(this);
        averagingWriter = new AnnotationAveragingWriter(outputInfo, provider);

        assert annotationFilename != null : "annotation filename must have been set";

        averagingWriter.setAnnotationFilename(annotationFilename);
        if (doIndels) {
            averagingWriter.setAggregateAllContexts(true);
        }
        readerIndexToGroupIndex = mode.getReaderIndexToGroupIndex();
        averagingWriter.setSampleIndexToGroupIndex(readerIndexToGroupIndex);
        final ObjectArrayList<ReadIndexStats> readIndexStats = mode.getReadIndexStats();

        groupComparisons = mode.getGroupComparisons();
        averagingWriter.setGroupComparisons(groupComparisons);
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
    private int referenceIndex = -1;
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
        averagingWriter.setGenome(genome);
        final char refBase = sampleCounts[0].referenceBase;
        if (refBase != 'C' && refBase != 'G') {
            return;
        }
        if (doIndels) {

            IndelCountOutputFormat.fillMethylationCountArrays(sampleCounts, list, position, refBase, mci, readerIndexToGroupIndex);
        } else {
            MethylationRateVCFOutputFormat.fillMethylationCountArrays(sampleCounts, list, position, refBase, mci, readerIndexToGroupIndex);
        }
        if (mci.eventCountAtSite < minimumEventThreshold) {
            return;
        }
        this.position = position;
        if (referenceIndex != this.referenceIndex) {
            this.referenceIndex = referenceIndex;
            chromosome = genome.getReferenceName(genomeReferenceIndex);
        }
        //   statWriter.writeRecord();
        averagingWriter.writeRecord();

    }

    public int getPosition() {
        return position;
    }

    public String getChromosome() {

        return chromosome;
    }

    @Override
    public void close() {


        averagingWriter.close();
    }

    @Override
    public void setGenome(RandomAccessSequenceInterface genome) {
        this.genome = genome;
    }

    @Override
    public void setGenomeReferenceIndex(int index) {
        genomeReferenceIndex = index;
    }


    public String[] getSamples() {
        return samples;
    }

    public MethylCountInfo getMci() {
        return mci;
    }

    @Override
    public void setMinimumEventThreshold(int minimumEventThreshold) {
        this.minimumEventThreshold = minimumEventThreshold;
    }


}
