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

import com.sun.corba.se.spi.logging.LogWrapperFactory;
import edu.cornell.med.icb.goby.algorithmic.algorithm.LogGCCorrectionWeight;
import edu.cornell.med.icb.goby.algorithmic.algorithm.SortedAnnotations;
import edu.cornell.med.icb.goby.algorithmic.data.Annotation;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import it.unimi.dsi.fastutil.ints.Int2ObjectAVLTreeMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.IntAVLTreeSet;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntSet;
import org.apache.commons.io.output.NullWriter;
import sun.rmi.runtime.Log;

import java.io.IOException;
import java.io.Writer;

/**
 * A VCF Writer that averages values of some fields over a set of annotations,
 * then writes the result to IGV file format.
 *
 * @author Fabien Campagne
 *         Date: 12/12/11
 *         Time: 1:39 PM
 */
public class VCFAveragingWriter extends VCFWriter {
    Writer outputWriter;
    private boolean initialized;
    private int numSamples;
    private int numFields;
    private String chromosome;
    private RandomAccessSequenceInterface genome;
    String[] chosenFormatFields;
    private MethylCountProvider provider;

    public VCFAveragingWriter(final Writer writer, RandomAccessSequenceInterface genome, MethylCountProvider provider) {
        super(new NullWriter());
        this.provider = provider;
        outputWriter = writer;
        this.genome = genome;
        initialized = false;
    }


    /**
     * Select the FORMAT fields whose values will be averaged per sample.
     *
     * @param selectedFormatFieldNames names for the FORMAT fields whose values will be averaged per sample and annotation.
     */
    public void selectFormatFields(final String[] selectedFormatFieldNames) {
        selectedFormatColumnIndices = new IntArrayList();
        for (final String fieldName : selectedFormatFieldNames) {
            selectedFormatColumnIndices.add(formatTypeToFormatFieldIndex.get(fieldName));
        }
    }

    private IntList selectedFormatColumnIndices;
    private String[] samples;


    @Override
    public void writeRecord() {
        init();
        provider.feedFrom(this);
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
            numFields = 2;

            // load annotations
            annotations.setGenome(this.genome);
            String filename = "/Users/nyasha/Dropbox/Professional/DMR-caller/testVCFAveragingWriter.txt";
            try {
                annotations.loadAnnotations(filename);
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
                for (final String sample : samples) {
                    StringBuilder columnName = new StringBuilder();
                    columnName.append("MR[");
                    columnName.append(sample);
                    columnName.append("]");
                    outputWriter.append(columnName.toString());
                    i++;
                    if (i != samples.length) {
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

    SortedAnnotations annotations = new SortedAnnotations();
    Int2ObjectMap<FormatFieldCounter> counterMap = new Int2ObjectAVLTreeMap<FormatFieldCounter>();

    private void addValuesAndAverage() {

        final String chromosome = provider.getChromosome().toString();
        final int pos = provider.getPosition();
        final int refIndex = genome.getReferenceIndex(chromosome);

        if (annotations.hasOverlappingAnnotations(refIndex, pos)) {
            System.out.println(chromosome + " position: " + pos + " has overlapping annotations");
            final IntAVLTreeSet overlappers = annotations.getValidAnnotationIndices();
            for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {
                for (int each : overlappers) {
                    // increment counters for each annotation overlapping at this position
                    FormatFieldCounter cntr;
                    if (counterMap.containsKey(each)) {
                        cntr = counterMap.get(each);
                    } else {
                        cntr = new FormatFieldCounter(each, numSamples);
                        counterMap.put(each, cntr);
                    }

                    cntr.unmethylatedCounter[sampleIndex] += provider.getC(sampleIndex, provider.getIndex());
                    cntr.methylatedCounter[sampleIndex] += provider.getCm(sampleIndex, provider.getIndex());
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

            try {
                outputWriter.append(annoOut.getChromosome());
                outputWriter.append("\t");
                outputWriter.append(String.valueOf(annoOut.getStart()));
                outputWriter.append("\t");
                outputWriter.append(String.valueOf(annoOut.getEnd())).append("\t");
                outputWriter.append(annoOut.getId());
                outputWriter.append("\t");
                for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {
                    temp.CalculateSampleMethylationRate(sampleIndex);
                    outputWriter.append(String.format("%g", temp.getMethylationRate()[sampleIndex]));
                    if (sampleIndex != (numSamples - 1)) {
                        outputWriter.append("\t");
                    }
                }
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
}
