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

package edu.cornell.med.icb.goby.algorithmic.algorithm.dmr;

import edu.cornell.med.icb.goby.algorithmic.data.GroupComparison;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import edu.cornell.med.icb.goby.stats.DummyObservationWriter;
import edu.cornell.med.icb.goby.stats.EmpiricalPValueEstimator;
import edu.cornell.med.icb.goby.stats.MethylCountProvider;
import edu.cornell.med.icb.goby.stats.RegionWriter;
import edu.cornell.med.icb.goby.stats.VCFWriter;
import edu.cornell.med.icb.goby.util.OutputInfo;
import edu.cornell.med.icb.goby.util.OutputInfoFromWriter;
import edu.cornell.med.icb.goby.util.dynoptions.DynamicOptionClient;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.output.NullWriter;
import org.apache.log4j.Logger;

import javax.swing.plaf.synth.Region;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Date;

/**
 * @Author: Nyasha Chambwe
 * @Date : 3/26/12
 * @Time: 2:26 PM
 */
public class RegionAveragingWriter extends VCFWriter implements RegionWriter {

    private ArrayList<GroupComparison> groupComparisons = new ArrayList<GroupComparison>();
    private RandomAccessSequenceInterface genome;
    private int[] sampleIndexToGroupIndex;
    private String[] contexts = {"CpG", "CpA", "CpC", "CpT", "CpN"};
    private boolean aggregateAllContexts;
    private Boolean estimateIntraGroupP;
    private Boolean estimateIntraGroupDifferences;
    private Boolean writeObservations;
    private Boolean writeCounts;
    private ObservationWriter obsWriter = new DummyObservationWriter();
    private boolean initialized;
    private MethylCountProvider provider;
    private OutputInfo outputInfo;



    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(RegionAveragingWriter.class);


    public static final DynamicOptionClient doc() {
        return doc;
    }

    public static final DynamicOptionClient doc = new DynamicOptionClient(RegionAveragingWriter.class,
            EmpiricalPValueEstimator.LOCAL_DYNAMIC_OPTIONS,
            "write-counts:boolean, when true write C and Cm for regions:false",
            "write-observations:boolean, when true write observations to disk: false",
            "contexts:string, coma delimited list of contexts for which to evaluate methylation rate. Contexts can be CpG, CpA,CpC,CpT,CpN. Default is CpG only:CpG"
    );

    public RegionAveragingWriter(final OutputInfo outputInfo, RandomAccessSequenceInterface genome,
                                 MethylCountProvider provider) {
        super(new NullWriter());
        final String contextString = doc.getString("contexts");
        final String[] contextTokens = contextString.split(",");
        if (contextTokens.length != 0) {
            LOG.info("registering user defined contexts: " + ObjectArrayList.wrap(contextTokens));
            contexts = contextTokens;
        }
        estimateIntraGroupDifferences = doc.getBoolean("estimate-intra-group-differences");
        estimateIntraGroupP = doc.getBoolean("estimate-empirical-P");
        writeCounts = doc.getBoolean("write-counts");
        writeObservations = doc.getBoolean("write-observations");

        if (estimateIntraGroupDifferences || estimateIntraGroupP) {
            String basename = FilenameUtils.removeExtension(outputInfo.getFilename());
            if (basename == null) {
                basename = Long.toString(new Date().getTime());
            }
            if (writeObservations) {
                String filename = basename + "-" + (estimateIntraGroupDifferences ? "null" : "test") + "-observations.tsv";
                try {
                    obsWriter = new ObservationWriter(new FileWriter(filename));
                    obsWriter.setHeaderIds(new String[]{"context", "chromosome", "start", "end", "annotation-id"});
                } catch (IOException e) {
                    LOG.error("Cannot open observation file for writing: " + filename);
                }

            }
        }

        this.provider = provider;
        this.outputInfo = outputInfo;
        if (!estimateIntraGroupDifferences) {
            //outputWriter = outputInfo.getPrintWriter();
        } else {
            //outputWriter = new NullWriter();
        }
        this.genome = genome;
        initialized = false;
        //processGroups = true;

    }

    public RegionAveragingWriter(OutputInfo outputInfo, MethylCountProvider provider) {
        this(outputInfo, null, provider);
    }

    public RegionAveragingWriter(final Writer writer, RandomAccessSequenceInterface genome, MethylCountProvider provider) {
        this(new OutputInfoFromWriter(writer), genome, provider);

    }

    @Override
    public void writeRecord() {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void close() {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void setGroupComparisons(ArrayList<GroupComparison> groupComparisons) {
        this.groupComparisons = groupComparisons;
    }

    @Override
    public void setGenome(RandomAccessSequenceInterface genome) {
        this.genome = genome;
    }

    @Override
    public void setAnnotationFilename(String annotationFilename) {
        //TODO
        // no need, trying to discover regions denovo either set annotaiton filename to null
        // remove method from Region Writer Interface
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void setSampleIndexToGroupIndex(int[] readerIndexToGroupIndex) {
        sampleIndexToGroupIndex = readerIndexToGroupIndex;
    }

    @Override
    public void setAggregateAllContexts(boolean aggregateAllContexts) {
        this.aggregateAllContexts = aggregateAllContexts;
        if (aggregateAllContexts) {
            contexts = new String[]{"ALL"};
        }
    }
}
