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

import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.readers.vcf.VCFParser;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceCache;
import edu.cornell.med.icb.goby.xml.MethylStats;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.log4j.Logger;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Marshaller;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;


/**
 * Calculates statistics about alignments produced from bisulfite converted reads.
 * Estimates the number of CpG observed in the alignment(s).
 */
public class MethylStatsMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "methyl-stats";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION =
            "Calculates statistics about alignments produced from bisulfite converted reads.";
    /**
     * Basename of the counts archives to analyze.
     */
    private String[] inputBasenames;
    /**
     * Output filename where to write the stats.
     */
    private String statsOuputFilename;

    private int numThreads;
    private double[] percentiles = new double[]{0.9, .75, .5, .1, .01};
    private int[] depths = new int[]{5, 10, 15, 20, 30};
    private RandomAccessSequenceCache genome;
    private int[] fragmentLengthBins = new int[]{1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250, 300, 350, 400, 450, 500};
    private static final Logger LOG = Logger.getLogger(MethylStatsMode.class);
    private String inputFilename;


    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }

    /**
     * Configure.
     *
     * @param args command line arguments
     * @return this object for chaining
     * @throws java.io.IOException error parsing
     * @throws com.martiansoftware.jsap.JSAPException
     *                             error parsing
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args) throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        inputFilename = jsapResult.getString("input");
        statsOuputFilename = jsapResult.getString("output");

        final String genomeBasename = jsapResult.getString("genome");
        if (genomeBasename != null) {
            genome = new RandomAccessSequenceCache();
            try {
                genome.load(genomeBasename);
            } catch (ClassNotFoundException e) {
                System.err.println("Cannot load genome. An exception occured. Details may be provided below.");
                e.printStackTrace();
                System.exit(1);
            }
        }
        fragmentLengthBins = stringToInts(jsapResult.getString("fragment-lengths"));
        depths = stringToInts(jsapResult.getString("depths"));
        return this;
    }

    private double[] stringToDoubles(String percentiles) {
        String[] values = percentiles.split(",");
        DoubleArrayList result = new DoubleArrayList();
        for (String value : values) {
            result.add(Double.parseDouble(value) / 100.0d);
        }
        return result.toDoubleArray();
    }

    private int[] stringToInts(String depths) {
        String[] values = depths.split(",");
        IntArrayList result = new IntArrayList();
        for (String value : values) {
            result.add(Integer.parseInt(value.trim()));
        }
        return result.toIntArray();
    }

    int positionOfNextCpG;
    int refIndexOfNextCpG;

    /**
     * Run the mode.
     */
    @Override
    public void execute() {
        final PrintWriter output;
        try {
            output = statsOuputFilename.equals("-") ? new PrintWriter(System.out) : new PrintWriter(new FileWriter(statsOuputFilename));
            final MethylStats backgroundStats = new MethylStats(depths, fragmentLengthBins);
            // First process the genome to find out the background distribution of CpGs:
            ProgressLogger pg = new ProgressLogger(LOG);

            pg.priority = org.apache.log4j.Level.INFO;
            pg.itemsName = "bases";
            pg.displayFreeMemory = true;
            pg.start("counting genome CpG sites ");
            for (int sequenceIndex = 0; sequenceIndex < genome.numberOfSequences(); sequenceIndex++) {
                final int genomeSequenceSize = genome.getSequenceSize(sequenceIndex);
                char previousBase = '\0';
                char currentBase = '\0';
                int positionCpG1 = 0;
                int positionCpG2 = 0;
                int fragmentLength = -1;
                for (int position = 0; position < genomeSequenceSize; position++) {

                    currentBase = genome.get(sequenceIndex, position);
                    pg.lightUpdate();
                    if (currentBase == 'G' && previousBase == 'C') {
                        // found a new CpG.
                        positionCpG1 = positionCpG2;
                        positionCpG2 = position;
                        fragmentLength = positionCpG2 - positionCpG1;
                        if (fragmentLength > 0) {
                            backgroundStats.genomeHasCpG(fragmentLength);

                        }
                    }
                    previousBase = currentBase;
                }
            }
            pg.done();
            System.out.printf("Found %d CpG sites in genome.%n", backgroundStats.getNumberCpGsInGenome());
            String VcfFilename = inputFilename;
            VCFParser vcfParser = new VCFParser(VcfFilename);
            try {
                vcfParser.readHeader();

                String samples[] = vcfParser.getColumnNamesUsingFormat();
                int[] sampleDepthGlobalFieldIndex = new int[samples.length];

                int i = 0;
                MethylStats[] methylStats = new MethylStats[samples.length];

                for (String sample : samples) {
                    sampleDepthGlobalFieldIndex[i] = vcfParser.getGlobalFieldIndex(sample, "GB");

                    methylStats[i] = backgroundStats.copy();
                    methylStats[i].sampleId = sample;
                    i++;
                }
                int positionGlobalFieldIndex = vcfParser.getGlobalFieldIndex("POS", "VALUE");
                int chromosomeGlobalFieldIndex = vcfParser.getGlobalFieldIndex("CHROM", "VALUE");


                MethylStats stats = backgroundStats.copy();
                pg = new ProgressLogger(LOG);
                pg.priority = org.apache.log4j.Level.INFO;
                pg.itemsName = "sites";
                pg.displayFreeMemory = true;
                pg.start();
                int numSamples = samples.length;


                while (vcfParser.hasNextDataLine()) {
                    String reference = vcfParser.getFieldValue(chromosomeGlobalFieldIndex).toString();
                    int referenceIndex = genome.getReferenceIndex(reference);
                    // VCF positions are 1-based, but Goby genome positions are 0-based, adjust here:
                    int sitePosition = vcfParser.getFieldNumValues(positionGlobalFieldIndex) - 1;
                    pg.lightUpdate();
                    if (isCpG(referenceIndex, sitePosition)) {
                        final int fragmentLength = calculateFragmentLength(reference, sitePosition);

                        for (i = 0; i < numSamples; i++) {
                            final int depthInSample = vcfParser.getFieldNumValues(sampleDepthGlobalFieldIndex[i]);

                            methylStats[i].observedInSample(depthInSample, fragmentLength);
                        }
                    }
                    vcfParser.next();
                }
                pg.stop();
                pg.done();
                final JAXBContext jc = JAXBContext.newInstance("edu.cornell.med.icb.goby.xml");
                final Marshaller m = jc.createMarshaller();
                i = 0;
                for (String sample : samples) {
                    m.marshal(methylStats[i], output);
                    i++;
                }
                output.flush();
            } catch (VCFParser.SyntaxException e) {
                System.err.println("exception reading VCF file: " + e);
                e.printStackTrace();
            }

        } catch (IOException e) {
            System.err.println("An error occured opening the output file. ");
            System.exit(1);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private boolean isCpG(final int referenceIndex, final int sitePosition) {
        if (sitePosition + 1 > referenceSequenceSize) {
            return false;
        }
        if (genome.get(referenceIndex, sitePosition) == 'C' &&
                genome.get(referenceIndex, sitePosition + 1) == 'G') {
            return true;
        }
        return false;
    }

    private int calculateFragmentLength(String reference, int sitePosition) {
        int referenceIndex = genome.getReferenceIndex(reference);
        if (referenceIndex != refIndexOfNextCpG || sitePosition >= this.positionOfNextCpG) {
            // we have change reference sequence or are past the next site.
            // need to update next CpG position
            referenceSequenceSize = genome.getSequenceSize(referenceIndex);

            char previousBase = '\0';
            char currentBase = '\0';
            int positionCpG1 = 0;
            int positionCpG2 = 0;
            int fragmentLength = -1;
            for (int pos = sitePosition + 2; pos < referenceSequenceSize; pos++) {
                currentBase = genome.get(referenceIndex, pos);
                if (currentBase == 'G' && previousBase == 'C') {
                    // found a new CpG.
                    positionCpG1 = positionCpG2;
                    positionCpG2 = pos;
                    fragmentLength = positionCpG2 - positionCpG1;
                    if (fragmentLength > 0) {
                        this.positionOfNextCpG = positionCpG2;
                        return positionOfNextCpG;

                    }
                }

                previousBase = currentBase;
            }


        }
        return positionOfNextCpG;

    }


    int referenceSequenceSize;

    private double sum(long[] array) {
        double sum = 0;
        int o = 0;
        for (long value : array) {

            sum += value;

        }
        return sum;
    }

    private double divide
            (
                    double v1,
                    double v2) {
        return v1 / v2;
    }

    private double sum
            (LongArrayList
                     array, int offset) {
        double sum = 0;
        int o = 0;
        for (long value : array) {
            if (o >= offset) {
                sum += value;
            }
            ++o;
        }
        return sum;
    }


    /**
     * Main method.
     *
     * @param args command line args.
     * @throws com.martiansoftware.jsap.JSAPException
     *                             error parsing
     * @throws java.io.IOException error parsing or executing.
     */
    public static void main
    (
            final String[] args) throws JSAPException, IOException {
        new MethylStatsMode().configure(args).execute();
    }
}