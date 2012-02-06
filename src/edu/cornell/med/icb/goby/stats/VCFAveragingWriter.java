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

import edu.cornell.med.icb.goby.algorithmic.algorithm.SortedAnnotations;
import edu.cornell.med.icb.goby.algorithmic.data.Annotation;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import edu.cornell.med.icb.goby.util.DynamicOptionClient;
import it.unimi.dsi.fastutil.ints.*;
import org.apache.commons.io.output.NullWriter;
import org.apache.commons.lang.ArrayUtils;

import java.io.IOException;
import java.io.Writer;
import java.lang.reflect.Array;

/**
 * A VCF Writer that averages values of some fields over a set of annotations,
 * then writes the result to IGV file format.
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
    private String chromosome;
    private RandomAccessSequenceInterface genome;
    String[] chosenFormatFields;
    private MethylCountProvider provider;
    private String annotationFilename = null;
    public static final DynamicOptionClient doc = new DynamicOptionClient(VCFAveragingWriter.class, "annotations:annotation filename:");
    private String[] groups;
    private int numGroups;
    private int[] sampleIndexToGroupIndex;
    private boolean processGroups;

    public void setGenome(RandomAccessSequenceInterface genome) {
        this.genome = genome;
    }

    public VCFAveragingWriter(final Writer writer, MethylCountProvider provider) {
        this(writer, null, provider);
    }

    public VCFAveragingWriter(final Writer writer, RandomAccessSequenceInterface genome, MethylCountProvider provider) {
        super(new NullWriter());
        this.provider = provider;
        outputWriter = writer;
        this.genome = genome;
        initialized = false;
        processGroups=false;
    }

    /**
     * Set the annotation filename.
     *
     * @param annotationFilename
     */
    public void setAnnotationFilename(String annotationFilename) {
        this.annotationFilename = annotationFilename;
    }



    private IntList selectedFormatColumnIndices;
    private String[] samples;


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
            groups = provider.getGroups();

            numSamples = samples.length;
           if(groups != null){
               processGroups=true;
            numGroups= groups.length;
           }
            // load annotations
            annotations.setGenome(this.genome);

            try {
                assert annotationFilename != null : "annotation filename cannot be null";
                annotations.loadAnnotations(annotationFilename);
                System.out.println("annotations loaded");
            } catch (IOException e) {
                System.err.println("An error occurred loading the annotation file. ");
                System.exit(1);
            }

            //write header
            try {
                //  IGV format - maintain fidelity
                outputWriter.append("Chromosome\tStart\tEnd\tFeature\t");
                int i = 0;

                String[] outputTracks= (String[]) ArrayUtils.addAll(samples, groups);
                for (final String trackName : outputTracks) {
                    StringBuilder columnName = new StringBuilder();
                    columnName.append("MR[");
                    columnName.append(trackName);
                    columnName.append("]");
                    outputWriter.append(columnName.toString());
                    i++;
                    if (i != outputTracks.length) {
                        outputWriter.append('\t');
                    }
                }
                outputWriter.append("\n");
                outputWriter.flush();
            } catch (IOException e) {
                throw new RuntimeException(e);
            }

        }
    }

    private SortedAnnotations annotations = new SortedAnnotations();
    private Int2ObjectMap<FormatFieldCounter> counterMap = new Int2ObjectAVLTreeMap<FormatFieldCounter>();

    private void addValuesAndAverage() {

        final String chromosome = provider.getChromosome().toString();
        final int pos = provider.getPosition();
        final int refIndex = genome.getReferenceIndex(chromosome);

        if (annotations.hasOverlappingAnnotations(refIndex, pos)) {
            System.out.println(chromosome + " position: " + pos + " has overlapping annotations");
            final IntAVLTreeSet validOverlappingAnnotations = annotations.getValidAnnotationIndices();

            // process each individual sample
            for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {
                for (int each : validOverlappingAnnotations) {
                    // increment counters for each annotation overlapping at this position
                    FormatFieldCounter cntr;
                    if (counterMap.containsKey(each)) {
                        cntr = counterMap.get(each);
                    } else {
                        cntr = new FormatFieldCounter(each, numSamples, numGroups);
                        counterMap.put(each, cntr);
                    }

                    cntr.unmethylatedCCounterPerSample[sampleIndex] += provider.getC(sampleIndex);
                    cntr.methylatedCCounterPerSample[sampleIndex] += provider.getCm(sampleIndex);
                    if(processGroups){
                    cntr.unmethylatedCcounterPerGroup[sampleIndexToGroupIndex[sampleIndex]] += provider.getC(sampleIndex);
                    cntr.methylatedCCounterPerGroup[sampleIndexToGroupIndex[sampleIndex]] += provider.getCm(sampleIndex);
                    }
                    cntr.numberOfSites[sampleIndex] += 1;
                    System.out.println("sample " + samples[sampleIndex] + " " + "position: " + pos + " " + cntr.toString(sampleIndex));
                }
            }



        } else {
            System.out.println("Did not find overlapping annotations for " + chromosome + " : position: " + pos);
        }


        IntSet currentAnnotations = counterMap.keySet();

        for (int anno : currentAnnotations) {
            outputMethylationRate(chromosome, pos, anno);
        }
    }

    private void outputMethylationRate(String chromosome, int pos, int anno) {

        if (annotations.pastChosenAnnotation(anno, chromosome, pos)) {
            // this annotation is ready to be written
            Annotation annoOut = annotations.getAnnotation(anno);
            FormatFieldCounter temp = counterMap.get(anno);
            StringBuilder lineToOutput= new StringBuilder("");
            try {
                lineToOutput.append(annoOut.getChromosome());
                lineToOutput.append("\t");
                lineToOutput.append(String.valueOf(annoOut.getStart()));
                lineToOutput.append("\t");
                lineToOutput.append(String.valueOf(annoOut.getEnd())).append("\t");
                lineToOutput.append(annoOut.getId());
                lineToOutput.append("\t");
                for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {
                    temp.CalculateSampleMethylationRate(sampleIndex);
                    lineToOutput.append(String.format("%g", temp.getMethylationRatePerSample()[sampleIndex]));
                    if ((sampleIndex != (numSamples - 1)) || numGroups > 0 ) {
                        lineToOutput.append("\t");
                    }
                }
                for(int groupIndex=0; groupIndex < numGroups; groupIndex++){
                    temp.CalculateGroupMethylationRate(groupIndex);
                    lineToOutput.append(String.format("%g", temp.getMethylationRatePerGroup()[groupIndex]));
                    if (groupIndex != (numGroups - 1)) {
                        lineToOutput.append("\t");
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

    @Override
    public void close() {
        if (!counterMap.isEmpty()) {
            IntSet currentAnnotations = counterMap.keySet();

            for (int anno : currentAnnotations) {
                outputMethylationRate(annotations.getAnnotationsLastChromosome(), Integer.MAX_VALUE, anno);
            }

            outWriter.close();
        }
    }
/*
    private int getFormatFieldValue(int orderInColumnIndices, int sampleIndex) {
        int formatFieldIndex = selectedFormatColumnIndices.getInt(orderInColumnIndices);
        return Integer.parseInt(getSampleValue(formatFieldIndex, sampleIndex).toString());
    }
*/


    /**
     * Select the FORMAT fields whose values will be averaged per sample.
     *
     * @param selectedFormatFieldNames names for the FORMAT fields whose values will be averaged per sample and annotation.
     */
   /** public void selectFormatFields(final String[] selectedFormatFieldNames) {
        selectedFormatColumnIndices = new IntArrayList();
        for (final String fieldName : selectedFormatFieldNames) {
            selectedFormatColumnIndices.add(formatTypeToFormatFieldIndex.get(fieldName));
        }
    }
    * @param readerIndexToGroupIndex
    */
    public void setSampleIndexToGroupIndex(int[] readerIndexToGroupIndex) {
        sampleIndexToGroupIndex= readerIndexToGroupIndex;
    }
}
