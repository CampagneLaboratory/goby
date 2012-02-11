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
import edu.cornell.med.icb.goby.algorithmic.data.Annotation;
import edu.cornell.med.icb.goby.algorithmic.data.GroupComparison;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import edu.cornell.med.icb.goby.util.DynamicOptionClient;
import it.unimi.dsi.fastutil.ints.*;
import org.apache.commons.io.output.NullWriter;
import org.apache.commons.lang.ArrayUtils;
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
public class VCFAveragingWriter extends VCFWriter {
    Writer outputWriter;
    private boolean initialized;
    private int numSamples;
    private RandomAccessSequenceInterface genome;
    String[] chosenFormatFields;
    private MethylCountProvider provider;
    private String annotationFilename = null;
    public static final DynamicOptionClient doc = new DynamicOptionClient(VCFAveragingWriter.class, "annotations:annotation filename:");
    private String[] groups;
    private int numGroups;
    private int[] sampleIndexToGroupIndex;
    private boolean fisherRInstalled;
    private SortedAnnotations annotations = new SortedAnnotations();
    private Int2ObjectMap<FormatFieldCounter> counterMap = new Int2ObjectAVLTreeMap<FormatFieldCounter>();

    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(VCFAveragingWriter.class);


    /**
     * An array that enumerates the sequence contexts under consideration
     */
    String[] contexts = {"CpG", "CpA", "CpC", "CpT", "CpN"};

    boolean processGroups;


    private IntList selectedFormatColumnIndices;
    private String[] samples;

    private ArrayList<GroupComparison> groupComparisons = new ArrayList<GroupComparison>();
    private boolean aggregateAllContexts;


    public VCFAveragingWriter(final Writer writer, MethylCountProvider provider) {
        this(writer, null, provider);
    }

    public VCFAveragingWriter(final Writer writer, RandomAccessSequenceInterface genome, MethylCountProvider provider) {
        super(new NullWriter());
        this.provider = provider;
        outputWriter = writer;
        this.genome = genome;
        initialized = false;
        processGroups = true;
    }


    @Override
    public void writeRecord() {
        init();
        addValuesAndAverage();
        provider.next();
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
        }
    }

    private void writeHeaders() {
        try {
            //  IGV format - maintain fidelity
            outputWriter.append("Chromosome\tStart\tEnd\tFeature");
            int i = 1;
            int j = 1;
            String[] outputTracks = (String[]) ArrayUtils.addAll(samples, groups);

            for (String trackName : outputTracks) {
                for (String context : contexts) {
                    j = writeRateColumn(i, j, outputTracks, trackName, context);
                }

                i++;
            }
            i = 1;
            j = 1;

            for (final GroupComparison comparison : groupComparisons) {
                for (String context : contexts) {
                    j = writeFisherColumn(i, j, comparison, context);
                }

                i++;
            }
            outputWriter.append('\n');
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private int writeFisherColumn(int i, int j, GroupComparison comparison, String context) throws IOException {
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

        j++;
        return j;
    }

    private int writeRateColumn(int i, int j, String[] outputTracks, String trackName, String context) throws IOException {
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

        j++;
        return j;
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
                    FormatFieldCounter cntr;
                    if (counterMap.containsKey(each)) {
                        cntr = counterMap.get(each);
                    } else {
                        cntr = new FormatFieldCounter(each, numSamples, numGroups, contexts.length);
                        counterMap.put(each, cntr);
                    }
                    cntr.incrementCounts(sampleIndex, sampleIndexToGroupIndex,
                            provider.getC(sampleIndex),
                            provider.getCm(sampleIndex), contextIndex, processGroups);

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
            FormatFieldCounter temp = counterMap.get(anno);

            StringBuilder lineToOutput = new StringBuilder("");
            try {
                lineToOutput.append(annoOut.getChromosome());
                lineToOutput.append("\t");
                lineToOutput.append(String.valueOf(annoOut.getStart()));
                lineToOutput.append("\t");
                lineToOutput.append(String.valueOf(annoOut.getEnd())).append("\t");
                lineToOutput.append(annoOut.getId());

                for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {

                    for (int currentContext = 0; currentContext < contexts.length; currentContext++) {

                        lineToOutput.append("\t");
                        temp.calculateSampleMethylationRate(sampleIndex, currentContext);
                        lineToOutput.append(String.format("%g", temp.getMethylationRatePerSample(currentContext, sampleIndex)));
                    }

                }


                for (int groupIndex = 0; groupIndex < numGroups; groupIndex++) {

                    for (int currentContext = 0; currentContext < contexts.length; currentContext++) {

                        lineToOutput.append("\t");

                        temp.calculateGroupMethylationRate(groupIndex, currentContext);
                        lineToOutput.append(String.format("%g", temp.getMethylationRatePerGroup(currentContext, groupIndex)));
                    }

                }
                for (final GroupComparison comparison : groupComparisons) {
                    final int indexGroup1 = comparison.indexGroup1;
                    final int indexGroup2 = comparison.indexGroup2;
                    double fisherP = Double.NaN;

                    for (int currentContext = 0; currentContext < contexts.length; currentContext++) {
                        final boolean ok = checkCounts(temp, currentContext);
                        if (ok) {
                            fisherP = fisherRInstalled ? FisherExactRCalculator.getFisherPValue(
                                    temp.getUnmethylatedCcounterPerGroup(currentContext, indexGroup1),
                                    temp.getMethylatedCCounterPerGroup(currentContext, indexGroup1),
                                    temp.getUnmethylatedCcounterPerGroup(currentContext, indexGroup2),
                                    temp.getMethylatedCCounterPerGroup(currentContext, indexGroup2)) : Double.NaN;

                        } else {
                            LOG.error(String.format("An exception was caught evaluation the Fisher Exact test P-value. " +
                                    "Details are provided below%n" + "[[%s  %s] [%s   %s]]",
                                    temp.getUnmethylatedCcounterPerGroup(currentContext, indexGroup1),
                                    temp.getMethylatedCCounterPerGroup(currentContext, indexGroup1),
                                    temp.getUnmethylatedCcounterPerGroup(currentContext, indexGroup2),
                                    temp.getMethylatedCCounterPerGroup(currentContext, indexGroup2)
                            ));
                        }
                        lineToOutput.append("\t");
                        lineToOutput.append(String.format("%g", fisherP));
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

    private boolean checkCounts(FormatFieldCounter tempCounter, int currentContext) {
        boolean ok = true;
        // detect if any count is negative (that's a bug)
        for (int indexGroup = 0; indexGroup < numGroups; indexGroup++) {

            if (tempCounter.getUnmethylatedCcounterPerGroup(currentContext, indexGroup) < 0) {
                ok = false;
            }

            if (tempCounter.getMethylatedCCounterPerGroup(currentContext, indexGroup) < 0) {
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

            outWriter.close();
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


    public void setGroupComparisons(ArrayList<GroupComparison> groupComparisons) {
        this.groupComparisons = groupComparisons;
    }

    public void setGenome(RandomAccessSequenceInterface genome) {
        this.genome = genome;
    }

    /**
     * Set the annotation filename.
     *
     * @param annotationFilename
     */
    public void setAnnotationFilename(String annotationFilename) {
        this.annotationFilename = annotationFilename;
    }

    /**
     * Set the sample index to group index array
     *
     * @param readerIndexToGroupIndex
     */
    public void setSampleIndexToGroupIndex(int[] readerIndexToGroupIndex) {
        sampleIndexToGroupIndex = readerIndexToGroupIndex;
    }


    public void setAggregateAllContexts(boolean aggregateAllContexts) {
        this.aggregateAllContexts = aggregateAllContexts;
        if (aggregateAllContexts) {
            contexts = new String[]{"ALL"};
        }
    }
}
