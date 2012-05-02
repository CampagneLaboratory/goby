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

import edu.cornell.med.icb.goby.R.GobyRengine;
import edu.cornell.med.icb.goby.algorithmic.algorithm.SortedAnnotations;
import edu.cornell.med.icb.goby.algorithmic.algorithm.dmr.EstimatedDistribution;
import edu.cornell.med.icb.goby.algorithmic.algorithm.dmr.ObservationWriter;
import edu.cornell.med.icb.goby.algorithmic.algorithm.dmr.Stat5StatisticAdaptor;
import edu.cornell.med.icb.goby.algorithmic.data.Annotation;
import edu.cornell.med.icb.goby.algorithmic.data.GroupComparison;
import edu.cornell.med.icb.goby.algorithmic.data.SamplePairEnumerator;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import edu.cornell.med.icb.goby.util.dynoptions.DynamicOptionClient;
import edu.cornell.med.icb.goby.util.OutputInfo;
import edu.cornell.med.icb.goby.util.OutputInfoFromWriter;
import edu.cornell.med.icb.goby.util.dynoptions.RegisterThis;
import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.io.output.NullWriter;
import org.apache.log4j.Logger;
import org.rosuda.JRI.Rengine;

import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Date;

/**
 * A VCF Writer that averages values of some fields over a set of annotations,
 * then writes the result to an IGV file format.
 *
 * @author Nyasha Chambwe
 * @author Fabien Campagne
 *         Date: 12/12/11
 *         Time: 1:39 PM
 */
public class AnnotationAveragingWriter extends VCFWriter implements RegionWriter {
    Writer outputWriter;
    private boolean initialized;
    private int numSamples;
    private RandomAccessSequenceInterface genome;
    String[] chosenFormatFields;
    private MethylCountProvider provider;
    private String annotationFilename = null;
    @RegisterThis
    public static final DynamicOptionClient doc = new DynamicOptionClient(AnnotationAveragingWriter.class,
            EmpiricalPValueEstimator.LOCAL_DYNAMIC_OPTIONS,
            "annotations:annotation filename:",
            "write-counts:boolean, when true write C and Cm for regions:false",
            "write-observations:boolean, when true write oservations to disk: false",
            "contexts:string, coma delimited list of contexts for which to evaluate methylation rate. Contexts can be CpG, CpA,CpC,CpT,CpN. Default is CpG only:CpG"
    );

    public static final DynamicOptionClient doc() {
        return doc;
    }

    private String[] groups;
    private int numGroups;
    private int[] sampleIndexToGroupIndex;
    private boolean fisherRInstalled;
    private SortedAnnotations annotations = new SortedAnnotations();
    private Int2ObjectMap<FormatFieldCounter> counterMap = new Int2ObjectAVLTreeMap<FormatFieldCounter>();

    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(AnnotationAveragingWriter.class);


    /**
     * An array that enumerates the sequence contexts under consideration
     */
    String[] contexts = {"CpG", "CpA", "CpC", "CpT", "CpN"};

    boolean processGroups;


    private IntList selectedFormatColumnIndices;
    private String[] samples;

    private ArrayList<GroupComparison> groupComparisons = new ArrayList<GroupComparison>();
    private boolean aggregateAllContexts;
    private Boolean estimateIntraGroupDifferences;
    private OutputInfo outputInfo;
    private Boolean writeCounts;
    private Boolean estimateIntraGroupP;
    private Boolean writeNumSites = true;
    final EmpiricalPValueEstimator empiricalPValueEstimator = new EmpiricalPValueEstimator();
    private final String[] identifiers = new String[5];
    private Boolean writeObservations;


    public AnnotationAveragingWriter(OutputInfo outputInfo, MethylCountProvider provider) {
        this(outputInfo, null, provider);
    }

    public AnnotationAveragingWriter(final Writer writer, RandomAccessSequenceInterface genome, MethylCountProvider provider) {
        this(new OutputInfoFromWriter(writer), genome, provider);

    }

    public void setWriteNumSites(boolean b) {
        writeNumSites = b;
    }

    public void setContexts(String[] strings) {
        contexts = strings;
    }


    public AnnotationAveragingWriter(final OutputInfo outputInfo, RandomAccessSequenceInterface genome, MethylCountProvider provider) {
        super(new NullWriter());
        String contextString = doc.getString("contexts");
        String[] contextTokens = contextString.split(",");
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
            outputWriter = outputInfo.getPrintWriter();
        } else {
            outputWriter = new NullWriter();
        }
        this.genome = genome;
        initialized = false;
        processGroups = true;

    }


    /**
     * Initialize variables for each declared sample. Write header to output. Executes exactly once per output file.
     */
    private void init() {
        if (!initialized) {
            initialized = true;
            samples = provider.getSamples();
            numSamples = samples.length;
            groups = provider.getGroups();

            if (groups == null) {
                processGroups = false;
                numGroups = 0;
            } else {
                if (groups.length < 1) {
                    System.err.println("Methylation format requires at least one group.");
                    System.exit(1);
                }
                numGroups = groups.length;
            }
            empiricalPValueEstimator.configure(contexts.length, doc);
            empiricalPValueEstimator.setGroupEnumerator(new SamplePairEnumerator(sampleIndexToGroupIndex, numSamples, numGroups, groupComparisons.size()));

            // load annotations
            annotations.setGenome(this.genome);
            try {
                assert annotationFilename != null : "annotation filename cannot be null";
                annotations.loadAnnotations(annotationFilename);
                LOG.info("annotations " + annotationFilename + " loaded.");
            } catch (IOException e) {
                LOG.warn("An error occurred loading the annotation file:  " + annotationFilename);
                return;
            }

            //write headers
            writeHeaders();

            //activate R 
            try {
                final Rengine rEngine = GobyRengine.getInstance().getRengine();
                fisherRInstalled = rEngine != null && rEngine.isAlive();
            } catch (java.lang.UnsatisfiedLinkError e) {
                System.out.println("Cannot initialize R");
                e.printStackTrace();
                throw e;
            }

            if (estimateIntraGroupDifferences) {
                empiricalPValueEstimator.setStatAdaptor(new Stat5StatisticAdaptor());
                empiricalPValueEstimator.setNullDistribution(new EstimatedDistribution(contexts.length, empiricalPValueEstimator.getStatAdaptor()));


            }
            if ((estimateIntraGroupDifferences || estimateIntraGroupP) && !(obsWriter instanceof DummyObservationWriter)) {
                empiricalPValueEstimator.getStatAdaptor().setObservationWriter(obsWriter);
            }
        }
    }

    private void writeHeaders() {
        try {
            //  IGV format - maintain fidelity
            outputWriter.append("Chromosome\tStart\tEnd\tFeature");

            if (writeCounts) {
                for (String context : contexts) {
                    for (String sample : samples) {
                        writeStatForSample(sample, context, "#C");
                    }
                }

                for (String context : contexts) {
                    for (String sample : samples) {
                        writeStatForSample(sample, context, "#Cm");
                    }
                }
            }

            for (String context : contexts) {
                for (String sample : samples) {
                    writeStatForSample(sample, context, "MR");
                }
            }

            if (groups != null) {

                if (writeCounts) {
                    for (String context : contexts) {
                        for (String group : groups) {
                            writeStatForSample(group, context, "#C");
                        }
                    }

                    for (String context : contexts) {
                        for (String group : groups) {
                            writeStatForSample(group, context, "#Cm");
                        }
                    }
                }
                for (String context : contexts) {
                    for (String group : groups) {
                        writeStatForSample(group, context, "MR");
                    }
                }
            }
            if (writeNumSites) {
                for (String context : contexts) {
                    for (String sample : samples) {
                        writeStatForSample(sample, context, "#Sites");
                    }
                }
                if (groups != null) {
                    for (String context : contexts) {
                        for (String group : groups) {
                            writeStatForSample(group, context, "#Sites");
                        }
                    }
                }
            }
            for (String context : contexts) {
                for (final GroupComparison comparison : groupComparisons) {

                    writeStatForGroupComparison(comparison, context, "fisherP");
                }
            }


            for (String context : contexts) {
                for (final GroupComparison comparison : groupComparisons) {

                    writeStatForGroupComparison(comparison, context, "deltaMR");
                }
            }
            if (estimateIntraGroupDifferences) {

                empiricalPValueEstimator.recordWithinGroupSamplePairs(groups);

            }
            if (estimateIntraGroupP) {
                for (final String context : contexts) {
                    for (final GroupComparison comparison : groupComparisons) {
                        empiricalPValueEstimator.recordBetweenGroupsSamplePairs(comparison);

                        writeStatForGroupComparison(comparison, context, "empiricalP");
                    }
                }
            }
            outputWriter.append('\n');
        } catch (IOException
                e) {
            throw new RuntimeException(e);
        }
    }

    private void writeStatForSample
            (String
                     trackName, String
                    context, String
                    statName) throws IOException {
        StringBuilder columnName = new StringBuilder();
        outputWriter.append('\t');
        columnName.append(statName + "[");
        columnName.append(trackName);
        columnName.append("]");
        if (!aggregateAllContexts) {
            columnName.append("[");
            columnName.append(context);
            columnName.append("]");
        }
        outputWriter.append(columnName.toString());
    }

    private void writeStatForGroupComparison
            (GroupComparison
                     comparison, String
                    context, String
                    statName) throws IOException {
        StringBuilder comparisonName = new StringBuilder();
        outputWriter.append('\t');
        comparisonName.append(statName + "[");
        comparisonName.append(comparison.nameGroup1);
        comparisonName.append("/");
        comparisonName.append(comparison.nameGroup2);
        comparisonName.append("]");
        if (!aggregateAllContexts) {
            comparisonName.append("[");
            comparisonName.append(context);
            comparisonName.append("]");
        }
        outputWriter.append(comparisonName.toString());
    }

    @Override
    public void writeRecord
            () {
        init();
        addValuesAndAverage();
        provider.next();
    }

    /**
     * Increment counters for methylated and non-methylated cytosines across
     * sequence contexts, samples and groups
     */
    private void addValuesAndAverage
    () {

        final String chromosome = provider.getChromosome().toString();
        final int pos = provider.getPosition();
        final int refIndex = genome.getReferenceIndex(chromosome);

        String currentContext = findGenomicContext(refIndex, pos);
        int contextIndex = codeIndex(currentContext);

        if (contextIndex == -1) {
            return;
        }
        if (annotations.hasOverlappingAnnotations(refIndex, pos)) {
            final IntAVLTreeSet validOverlappingAnnotations = annotations.getValidAnnotationIndices();

            // process each individual sample
            for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {
                for (int each : validOverlappingAnnotations) {
                    // increment counters for each annotation overlapping at this position
                    FormatFieldCounter cntr = counterMap.get(each);

                    if (cntr == null) {
                        cntr = new FormatFieldCounter(each, numSamples, numGroups, contexts);
                        counterMap.put(each, cntr);
                    }
                    cntr.incrementCounts(sampleIndex, sampleIndexToGroupIndex,
                            provider.getC(sampleIndex),
                            provider.getCm(sampleIndex), contextIndex);

                    if (LOG.isTraceEnabled()) {
                        LOG.debug("sample " + samples[sampleIndex] + " " + "position: " + pos);
                    }
                }
            }
        } else {
            if (LOG.isTraceEnabled()) {
                LOG.trace("Did not find overlapping annotations for " + chromosome + " : position: " + pos);
            }
        }


        final IntSet currentAnnotations = counterMap.keySet();

        for (final int annot : currentAnnotations) {
            buildAnnotationRecordForOutput(chromosome, pos, annot);
        }
    }

    private int codeIndex
            (String
                     currentContext) {
        int contextIndex = -1;
        for (int i = 0; i < contexts.length; i++) {
            if (currentContext.equals(contexts[i])) {
                contextIndex = i;
            }
        }
        if (aggregateAllContexts) {
            //all counts go to single "ALL" context.
            return 0;
        }
        if (contextIndex == -1) {
            if (LOG.isTraceEnabled()) LOG.trace("context was not recognized: " + currentContext);
        }
        return contextIndex;
    }

    private ObservationWriter obsWriter = new DummyObservationWriter();

    private void buildAnnotationRecordForOutput(String chromosome, int pos, int anno) {

        if (annotations.pastChosenAnnotation(anno, chromosome, pos)) {
            // this annotation is ready to be written
            Annotation annoOut = annotations.getAnnotation(anno);
            FormatFieldCounter counter = counterMap.get(anno);

            StringBuilder lineToOutput = new StringBuilder("");
            try {
                if (writeObservations) {
                    identifiers[0] = "context"; // will be filled below.
                    identifiers[1] = annoOut.getChromosome();
                    identifiers[2] = String.valueOf(annoOut.getStart());
                    identifiers[3] = String.valueOf(annoOut.getEnd());
                    identifiers[4] = annoOut.getId();
                    obsWriter.setElementIds(identifiers);
                }
                lineToOutput.append(annoOut.getChromosome());
                lineToOutput.append("\t");
                lineToOutput.append(String.valueOf(annoOut.getStart()));
                lineToOutput.append("\t");
                lineToOutput.append(String.valueOf(annoOut.getEnd())).append("\t");
                lineToOutput.append(annoOut.getId());

                if (writeCounts) {
                    for (int currentContext = 0; currentContext < contexts.length; currentContext++) {
                        for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {
                            lineToOutput.append("\t");
                            final int unMethylatedCCounterPerSample = counter.getUnmethylatedCCountPerSample(currentContext, sampleIndex);
                            lineToOutput.append(unMethylatedCCounterPerSample);
                        }
                    }

                    for (int currentContext = 0; currentContext < contexts.length; currentContext++) {
                        for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {
                            lineToOutput.append("\t");
                            final int methylatedCCounterPerSample = counter.getMethylatedCCountPerSample(currentContext, sampleIndex);
                            lineToOutput.append(methylatedCCounterPerSample);
                        }
                    }
                }

                for (int currentContext = 0; currentContext < contexts.length; currentContext++) {
                    for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {
                        lineToOutput.append("\t");
                        final double methylationRatePerSample = counter.getMethylationRatePerSample(currentContext, sampleIndex);
                        //     System.out.printf("context=%s sample=%s mr=%g %n", contexts[currentContext], samples[sampleIndex], methylationRatePerSample);
                        lineToOutput.append(formatDouble(methylationRatePerSample));
                    }
                }

                if (writeCounts) {
                    for (int currentContext = 0; currentContext < contexts.length; currentContext++) {
                        for (int groupIndex = 0; groupIndex < numGroups; groupIndex++) {
                            lineToOutput.append("\t");
                            final int unMethylatedCCounterPerGroup = counter.getUnmethylatedCcountPerGroup(currentContext, groupIndex);
                            lineToOutput.append(unMethylatedCCounterPerGroup);
                        }
                    }

                    for (int currentContext = 0; currentContext < contexts.length; currentContext++) {
                        for (int groupIndex = 0; groupIndex < numGroups; groupIndex++) {
                            lineToOutput.append("\t");
                            final int methylatedCCounterPerGroup = counter.getMethylatedCCountPerGroup(currentContext, groupIndex);
                            lineToOutput.append(methylatedCCounterPerGroup);
                        }
                    }
                }

                for (int currentContext = 0; currentContext < contexts.length; currentContext++) {

                    for (int groupIndex = 0; groupIndex < numGroups; groupIndex++) {

                        lineToOutput.append("\t");
                        lineToOutput.append(formatDouble(counter.getMethylationRatePerGroup(currentContext, groupIndex)));
                    }
                }
                if (writeNumSites) {
                    for (int currentContext = 0; currentContext < contexts.length; currentContext++) {

                        for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {

                            lineToOutput.append("\t");
                            final int numSitesPerSample = counter.getNumberOfSitesPerSample(currentContext, sampleIndex);
                            lineToOutput.append(numSitesPerSample);
                        }
                    }
                    for (int currentContext = 0; currentContext < contexts.length; currentContext++) {

                        for (int groupIndex = 0; groupIndex < numGroups; groupIndex++) {

                            lineToOutput.append("\t");
                            final int numSitesPerGroup = counter.getNumberOfSitesPerGroup(currentContext, groupIndex);
                            lineToOutput.append(numSitesPerGroup);
                        }
                    }
                }
                for (int currentContext = 0; currentContext < contexts.length; currentContext++) {

                    for (final GroupComparison comparison : groupComparisons) {
                        final int indexGroup1 = comparison.indexGroup1;
                        final int indexGroup2 = comparison.indexGroup2;
                        double fisherP = Double.NaN;

                        final boolean ok = checkCounts(counter, currentContext);
                        if (ok) {
                            fisherP = fisherRInstalled ? FisherExactRCalculator.getFisherPValue(
                                    counter.getUnmethylatedCcountPerGroup(currentContext, indexGroup1),
                                    counter.getMethylatedCCountPerGroup(currentContext, indexGroup1),
                                    counter.getUnmethylatedCcountPerGroup(currentContext, indexGroup2),
                                    counter.getMethylatedCCountPerGroup(currentContext, indexGroup2)) : Double.NaN;

                        } else {
                            LOG.error(String.format("An exception was caught evaluation the Fisher Exact test P-value. " +
                                    "Details are provided below%n" + "[[%s  %s] [%s   %s]]",
                                    counter.getUnmethylatedCcountPerGroup(currentContext, indexGroup1),
                                    counter.getMethylatedCCountPerGroup(currentContext, indexGroup1),
                                    counter.getUnmethylatedCcountPerGroup(currentContext, indexGroup2),
                                    counter.getMethylatedCCountPerGroup(currentContext, indexGroup2)
                            ));
                        }
                        lineToOutput.append("\t");
                        lineToOutput.append(formatDouble(fisherP));
                    }
                }
                for (int currentContext = 0; currentContext < contexts.length; currentContext++) {

                    for (final GroupComparison comparison : groupComparisons) {
                        final int indexGroup1 = comparison.indexGroup1;
                        final int indexGroup2 = comparison.indexGroup2;
                        final double deltaMR = Math.abs(counter.getMethylationRatePerGroup(currentContext, indexGroup1)
                                - counter.getMethylationRatePerGroup(currentContext, indexGroup2));

                        lineToOutput.append("\t");
                        lineToOutput.append(formatDouble(deltaMR));
                    }
                }
                if (estimateIntraGroupDifferences) {
                    obsWriter.setTypeOfPair(ObservationWriter.TypeOfPair.WITHIN_GROUP_PAIR);
                    for (int currentContext = 0; currentContext < contexts.length; currentContext++) {
                        identifiers[0] = contexts[currentContext];
                        for (final GroupComparison comparison : groupComparisons) {
                            obsWriter.setNullComparison(comparison.nameGroup1);
                            empiricalPValueEstimator.estimateNullDensity(currentContext, comparison.indexGroup1, counter);
                            obsWriter.setNullComparison(comparison.nameGroup2);
                            empiricalPValueEstimator.estimateNullDensity(currentContext, comparison.indexGroup2, counter);
                        }
                    }
                }
                if (estimateIntraGroupP) {
                    obsWriter.setTypeOfPair(ObservationWriter.TypeOfPair.BETWEEN_GROUP_PAIR);
                    for (int contextIndex = 0; contextIndex < contexts.length; contextIndex++) {
                        identifiers[0] = contexts[contextIndex];
                        for (final GroupComparison comparison : groupComparisons) {

                            obsWriter.setComparison(comparison);
                            final double p = empiricalPValueEstimator.estimateEmpiricalPValue(contextIndex, comparison, counter);
                            lineToOutput.append("\t");
                            lineToOutput.append(formatDouble(p));
                        }
                    }
                }
                outputWriter.append(lineToOutput.toString());
                outputWriter.append("\n");
                counterMap.remove(anno);
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }
    }


    /**
     * Format double, rendering NaN as empty string.
     *
     * @param value
     * @return
     */
    private String formatDouble
    (
            double value) {
        if (value != value) {
            // value is NaN
            return "";

        } else {
            return String.format("%g", value);
        }
    }

    private boolean checkCounts
            (FormatFieldCounter
                     tempCounter, int currentContext) {
        boolean ok = true;
        // detect if any count is negative (that's a bug)
        for (int indexGroup = 0; indexGroup < numGroups; indexGroup++) {

            if (tempCounter.getUnmethylatedCcountPerGroup(currentContext, indexGroup) < 0) {
                ok = false;
            }

            if (tempCounter.getMethylatedCCountPerGroup(currentContext, indexGroup) < 0) {
                ok = false;
            }
        }

        return ok;
    }


    @Override
    public void close
            () {
        if (!counterMap.isEmpty()) {
            IntSet currentAnnotations = counterMap.keySet();

            for (int anno : currentAnnotations) {
                buildAnnotationRecordForOutput(annotations.getAnnotationsLastChromosome(), Integer.MAX_VALUE, anno);
            }
        }
        outWriter.close();
        IOUtils.closeQuietly(outputWriter);

        if (estimateIntraGroupDifferences) {
            // when estimating intra-group differences, we  serialize the estimator to the output.
            try {
                LOG.debug("Storing density to filename: " + outputInfo.getFilename());
                EstimatedDistribution.store(empiricalPValueEstimator.getNullDistribution(), outputInfo.getFilename());
            } catch (IOException e) {
                LOG.error("Unable to write estimator to file", e);
            }
        }
        if (obsWriter != null) {
            obsWriter.close();

        }

    }


    private String findGenomicContext
            (
                    int referenceIndex,
                    int position) {
        int referenceLength = genome.getLength(referenceIndex);
        int zeroBasedPos = position - 1;
        char currentBase = genome.get(referenceIndex, zeroBasedPos);

        char nextBase = '?';
        String tempContext = new StringBuilder().append('C').append('p').toString();
        char concatBase = 'N';

        if (currentBase == 'C') {
            if (referenceLength == position) {
                return Character.toString(currentBase);
            }
            nextBase = genome.get(referenceIndex, (zeroBasedPos + 1));
            concatBase = nextBase;
        } else {
            if (currentBase == 'G') {
                if (zeroBasedPos == 0) {
                    return Character.toString(currentBase);
                }
                nextBase = genome.get(referenceIndex, (zeroBasedPos - 1));
                switch (nextBase) {
                    case 'C':
                        concatBase = 'G';
                        break;
                    case 'A':
                        concatBase = 'T';
                        break;
                    case 'T':
                        concatBase = 'A';
                        break;
                    case 'G':
                        concatBase = 'C';
                        break;
                }
            }
        }
        tempContext = tempContext.concat(Character.toString(concatBase));
        return tempContext;
    }


    @Override
    public void setGroupComparisons
            (ArrayList<GroupComparison> groupComparisons) {
        this.groupComparisons = groupComparisons;

    }

    @Override
    public void setGenome
            (RandomAccessSequenceInterface
                     genome) {
        this.genome = genome;
    }

    /**
     * Set the annotation filename.
     *
     * @param annotationFilename
     */
    @Override
    public void setAnnotationFilename
    (String
             annotationFilename) {
        this.annotationFilename = annotationFilename;
    }

    /**
     * Set the sample index to group index array
     *
     * @param readerIndexToGroupIndex
     */
    @Override
    public void setSampleIndexToGroupIndex
    (
            int[] readerIndexToGroupIndex) {
        sampleIndexToGroupIndex = readerIndexToGroupIndex;
    }


    @Override
    public void setAggregateAllContexts
            (
                    boolean aggregateAllContexts) {
        this.aggregateAllContexts = aggregateAllContexts;
        if (aggregateAllContexts) {
            contexts = new String[]{"ALL"};
        }
    }

    public void setWriteCounts
            (
                    boolean b) {
        writeCounts = b;
    }
}
