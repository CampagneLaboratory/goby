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
import edu.cornell.med.icb.goby.algorithmic.algorithm.dmr.DensityEstimator;
import edu.cornell.med.icb.goby.algorithmic.data.Annotation;
import edu.cornell.med.icb.goby.algorithmic.data.GroupComparison;
import edu.cornell.med.icb.goby.algorithmic.data.IntraGroupEnumerator;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import edu.cornell.med.icb.goby.util.DynamicOptionClient;
import edu.cornell.med.icb.goby.util.IOUtil;
import edu.cornell.med.icb.goby.util.OutputInfo;
import edu.cornell.med.icb.goby.util.OutputInfoFromWriter;
import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.apache.commons.io.IOUtils;
import org.apache.commons.io.output.NullWriter;
import org.apache.log4j.Logger;
import org.rosuda.JRI.Rengine;

import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;

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
    public static final DynamicOptionClient doc = new DynamicOptionClient(AnnotationAveragingWriter.class,
            "annotations:annotation filename:",
            "write-counts:boolean, when true write C and Cm for regions:false",
            "estimate-intra-group-differences: boolean, true indicates that pair-wise differences for sample in the same group should be tallied and written to the output. False indicates regular output.:false");
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


    public AnnotationAveragingWriter(OutputInfo outputInfo, MethylCountProvider provider) {
        this(outputInfo, null, provider);
    }

    public AnnotationAveragingWriter(final Writer writer, RandomAccessSequenceInterface genome, MethylCountProvider provider) {
        this(new OutputInfoFromWriter(writer), genome, provider);

    }

    public AnnotationAveragingWriter(final OutputInfo outputInfo, RandomAccessSequenceInterface genome, MethylCountProvider provider) {
        super(new NullWriter());
        estimateIntraGroupDifferences = doc.getBoolean("estimate-intra-group-differences");
        writeCounts = doc.getBoolean("write-counts");
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
            } else {
                if (groups.length < 1) {
                    System.err.println("Methylation format requires at least one group.");
                    System.exit(1);
                }
                numGroups = groups.length;
            }

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
                estimator = new DensityEstimator(contexts.length);
            }
        }
    }

    private void writeHeaders() {
        try {
            //  IGV format - maintain fidelity
            outputWriter.append("Chromosome\tStart\tEnd\tFeature");

            for (String context : contexts) {
                for (String sample : samples) {
                    if (writeCounts) {
                        writeCColumn(sample, context);
                        writeCmColumn(sample, context);
                        writeNumSitesColumn(sample, context);
                    }
                    writeRateColumn(sample, context);
                }
            }
            if (groups != null) {
                for (String context : contexts) {
                    for (String group : groups) {
                        if (writeCounts) {
                            writeCColumn(group, context);
                            writeCmColumn(group, context);
                            writeNumSitesColumn(group, context);
                        }
                        writeRateColumn(group, context);
                    }
                }
            }
            for (String context : contexts) {
                for (final GroupComparison comparison : groupComparisons) {

                    writeFisherColumn(comparison, context);
                }

            }


            for (String context : contexts) {
                for (final GroupComparison comparison : groupComparisons) {

                    writeDeltaMRColumn(comparison, context);
                }
            }
            if (estimateIntraGroupDifferences) {
                groupEnumerator = new IntraGroupEnumerator(sampleIndexToGroupIndex, numSamples, numGroups);

                for (final GroupComparison comparison : groupComparisons) {
                    groupEnumerator.recordPairForGroup(comparison.indexGroup1);
                    groupEnumerator.recordPairForGroup(comparison.indexGroup2);
                }
            }
            outputWriter.append('\n');
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private IntraGroupEnumerator groupEnumerator;

    private void writeCColumn(String trackName, String context) throws IOException {
        StringBuilder columnName = new StringBuilder();
        outputWriter.append('\t');
        columnName.append("#C[");
        columnName.append(trackName);
        columnName.append("]");
        if (!aggregateAllContexts) {
            columnName.append("[");
            columnName.append(context);
            columnName.append("]");
        }
        outputWriter.append(columnName.toString());

    }

    private void writeCmColumn(String trackName, String context) throws IOException {
        StringBuilder columnName = new StringBuilder();
        outputWriter.append('\t');
        columnName.append("#Cm[");
        columnName.append(trackName);
        columnName.append("]");
        if (!aggregateAllContexts) {
            columnName.append("[");
            columnName.append(context);
            columnName.append("]");
        }
        outputWriter.append(columnName.toString());
    }

    private void writeNumSitesColumn(String trackName, String context) throws IOException {
        StringBuilder columnName = new StringBuilder();
        outputWriter.append('\t');
        columnName.append("#Sites[");
        columnName.append(trackName);
        columnName.append("]");
        if (!aggregateAllContexts) {
            columnName.append("[");
            columnName.append(context);
            columnName.append("]");
        }
        outputWriter.append(columnName.toString());
    }

    private void writeFisherColumn(GroupComparison comparison, String context) throws IOException {
        StringBuilder comparisonName = new StringBuilder();
        outputWriter.append('\t');
        comparisonName.append("fisherP[");
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

    private void writeDeltaMRColumn(GroupComparison comparison, String context) throws IOException {
        StringBuilder comparisonName = new StringBuilder();
        outputWriter.append('\t');
        comparisonName.append("deltaMR[");
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

    private void writeRateColumn(String trackName, String context) throws IOException {
        StringBuilder columnName = new StringBuilder();
        outputWriter.append('\t');
        columnName.append("MR[");
        columnName.append(trackName);
        columnName.append("]");
        if (!aggregateAllContexts) {
            columnName.append("[");
            columnName.append(context);
            columnName.append("]");
        }
        outputWriter.append(columnName.toString());

    }

     @Override
    public void writeRecord() {
        init();
        addValuesAndAverage();
        provider.next();
    }
    /**
     * Increment counters for methylated and non-methylated cytosines across
     * sequence contexts, samples and groups
     */
    private void addValuesAndAverage() {

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

                    LOG.debug("sample " + samples[sampleIndex] + " " + "position: " + pos);
                }
            }
        } else {
            LOG.debug("Did not find overlapping annotations for " + chromosome + " : position: " + pos);
        }


        IntSet currentAnnotations = counterMap.keySet();

        for (final int annot : currentAnnotations) {
            buildAnnotationRecordForOutput(chromosome, pos, annot);
        }
    }

    private int codeIndex(String currentContext) {
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
            LOG.warn("context was not recognized: " + currentContext);
        }
        return contextIndex;
    }

    private void buildAnnotationRecordForOutput(String chromosome, int pos, int anno) {

        if (annotations.pastChosenAnnotation(anno, chromosome, pos)) {
            // this annotation is ready to be written
            Annotation annoOut = annotations.getAnnotation(anno);
            FormatFieldCounter counter = counterMap.get(anno);

            StringBuilder lineToOutput = new StringBuilder("");
            try {
                lineToOutput.append(annoOut.getChromosome());
                lineToOutput.append("\t");
                lineToOutput.append(String.valueOf(annoOut.getStart()));
                lineToOutput.append("\t");
                lineToOutput.append(String.valueOf(annoOut.getEnd())).append("\t");
                lineToOutput.append(annoOut.getId());
                for (int currentContext = 0; currentContext < contexts.length; currentContext++) {

                    for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {
                        if (writeCounts) {
                            lineToOutput.append("\t");
                            final int unMethylatedCCounterPerSample = counter.getUnmethylatedCCountPerSample(currentContext, sampleIndex);
                            lineToOutput.append(unMethylatedCCounterPerSample);

                            lineToOutput.append("\t");
                            final int methylatedCCounterPerSample = counter.getMethylatedCCountPerSample(currentContext, sampleIndex);
                            lineToOutput.append(methylatedCCounterPerSample);

                            lineToOutput.append("\t");
                            final int numSitesPerSample = counter.getNumberOfSitesPerSample(currentContext, sampleIndex);
                            lineToOutput.append(numSitesPerSample);
                        }

                        lineToOutput.append("\t");
                        final double methylationRatePerSample = counter.getMethylationRatePerSample(currentContext, sampleIndex);
                        //     System.out.printf("context=%s sample=%s mr=%g %n", contexts[currentContext], samples[sampleIndex], methylationRatePerSample);
                        lineToOutput.append(formatDouble(methylationRatePerSample));
                    }
                }

                for (int currentContext = 0; currentContext < contexts.length; currentContext++) {

                    for (int groupIndex = 0; groupIndex < numGroups; groupIndex++) {
                        if (writeCounts) {
                            lineToOutput.append("\t");
                            final int unMethylatedCCounterPerGroup = counter.getUnmethylatedCcountPerGroup(currentContext, groupIndex);
                            lineToOutput.append(unMethylatedCCounterPerGroup);

                            lineToOutput.append("\t");
                            final int methylatedCCounterPerGroup = counter.getMethylatedCCountPerGroup(currentContext, groupIndex);
                            lineToOutput.append(methylatedCCounterPerGroup);


                            lineToOutput.append("\t");
                            final int numSitesPerGroup = counter.getNumberOfSitesPerGroup(currentContext, groupIndex);
                            lineToOutput.append(numSitesPerGroup);

                        }
                        lineToOutput.append("\t");
                        lineToOutput.append(formatDouble(counter.getMethylationRatePerGroup(currentContext, groupIndex)));
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
                        double deltaMR = Double.NaN;

                        deltaMR = Math.abs(counter.getMethylationRatePerGroup(currentContext, indexGroup1) - counter.getMethylationRatePerGroup(currentContext, indexGroup2));

                        lineToOutput.append("\t");
                        lineToOutput.append(formatDouble(deltaMR));
                    }
                }
                if (estimateIntraGroupDifferences) {
                    for (int currentContext = 0; currentContext < contexts.length; currentContext++) {

                        for (final GroupComparison comparison : groupComparisons) {

                            collectWithinGroupEstimates(currentContext, comparison.indexGroup1, counter);
                            collectWithinGroupEstimates(currentContext, comparison.indexGroup2, counter);
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

    private void collectWithinGroupEstimates(final int contextIndex, final int groupIndex,
                                             final FormatFieldCounter counter) {
        if (estimateIntraGroupDifferences) {
            // enumerate sample pairs that belong to the group of interest:
            final ObjectArrayList<SamplePair> pairs = groupEnumerator.getPairs(groupIndex);
            for (final SamplePair next : pairs) {
                final int Cma = counter.getMethylatedCCountPerSample(contextIndex, next.sampleIndexA);
                final int Ca = counter.getUnmethylatedCCountPerSample(contextIndex, next.sampleIndexA);
                final int Cmb = counter.getMethylatedCCountPerSample(contextIndex, next.sampleIndexB);
                final int Cb = counter.getUnmethylatedCCountPerSample(contextIndex, next.sampleIndexB);
                if ((Cma + Ca) == 0 || (Cmb + Cb) == 0) {
                    System.out.printf("Zero in one intra-group sample for %d %d %d %d samplexIndexA=%d sampleIndexB=%d %n",
                            Cma, Ca, Cmb, Cb,next.sampleIndexA, next.sampleIndexB);
                } else {
                    estimator.observe(contextIndex, Cma, Ca, Cmb, Cb);
                }
            }
        }
    }

    DensityEstimator estimator;

    /**
     * Format double, rendering NaN as empty string.
     *
     * @param value
     * @return
     */
    private String formatDouble(double value) {
        if (value != value) {
            // value is NaN
            return "";

        } else {
            return String.format("%g", value);
        }
    }

    private boolean checkCounts(FormatFieldCounter tempCounter, int currentContext) {
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
    public void close() {
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
                DensityEstimator.store(estimator, outputInfo.getFilename());
            } catch (IOException e) {
                LOG.error("Unable to write estimator to file", e);
            }
        }
    }


    private String findGenomicContext(int referenceIndex, int position) {
        int zeroBasedPos = position - 1;
        char currentBase = genome.get(referenceIndex, zeroBasedPos);
        int referenceLength = genome.getLength(referenceIndex);
        char nextBase = '?';
        String tempContext = new StringBuilder().append('C').append('p').toString();
        char concatBase = '?';

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
    public void setGroupComparisons(ArrayList<GroupComparison> groupComparisons) {
        this.groupComparisons = groupComparisons;
    }

    @Override
    public void setGenome(RandomAccessSequenceInterface genome) {
        this.genome = genome;
    }

    /**
     * Set the annotation filename.
     *
     * @param annotationFilename
     */
    @Override
    public void setAnnotationFilename(String annotationFilename) {
        this.annotationFilename = annotationFilename;
    }

    /**
     * Set the sample index to group index array
     *
     * @param readerIndexToGroupIndex
     */
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

    public void setWriteCounts(boolean b) {
        writeCounts = b;
    }
}
