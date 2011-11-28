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
import javax.xml.bind.JAXBException;
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
    private int[] depths = new int[]{5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 300, 500};
    private RandomAccessSequenceCache genome;
    private int[] fragmentLengthBins = new int[]{1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250, 300, 350, 400, 450, 500};
    private static final Logger LOG = Logger.getLogger(MethylStatsMode.class);
    private String inputFilename;
    private static final boolean QUICK = false;

    private boolean doFragments = true;
    private int[] baseCallGlobalFieldIndex;

    private int[] convertedCystosineGlobalFieldIndex;
    private int[] unconvertedCystosineGlobalFieldIndex;


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
        doFragments = jsapResult.getBoolean("fragments");
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

    int positionOfNextCpG = -1;
    int refIndexOfNextCpG = -1;
    int[] sampleDepthGlobalFieldIndex;
    int[] methylationRateGlobalFieldIndex;

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
            if (doFragments) {
                System.out.printf("Pre-processing genome.%n");
                scanOneStrand(backgroundStats, '+');
                scanOneStrand(backgroundStats, '-');
                System.out.printf("Found %d CpG sites in genome.%n", backgroundStats.getNumberCpGsInGenome());

            } else {
                System.out.println("Fragment analysis is not activated (activate with --fragments/-f).");
            }
            String VcfFilename = inputFilename;
            VCFParser vcfParser = new VCFParser(VcfFilename);
            try {
                vcfParser.readHeader();

                String samples[] = vcfParser.getColumnNamesUsingFormat();
                sampleDepthGlobalFieldIndex = new int[samples.length];
                methylationRateGlobalFieldIndex = new int[samples.length];
                convertedCystosineGlobalFieldIndex = new int[samples.length];
                unconvertedCystosineGlobalFieldIndex = new int[samples.length];
                baseCallGlobalFieldIndex = new int[samples.length];

                int i = 0;
                MethylStats[] methylStats = new MethylStats[samples.length];

                for (String sample : samples) {
                    sampleDepthGlobalFieldIndex[i] = vcfParser.getGlobalFieldIndex(sample, "GB");
                    methylationRateGlobalFieldIndex[i] = vcfParser.getGlobalFieldIndex(sample, "MR");
                    convertedCystosineGlobalFieldIndex[i] = vcfParser.getGlobalFieldIndex(sample, "C");
                    unconvertedCystosineGlobalFieldIndex[i] = vcfParser.getGlobalFieldIndex(sample, "Cm");
                    if (convertedCystosineGlobalFieldIndex[i] == -1 || unconvertedCystosineGlobalFieldIndex[i] == -1) {
                        System.err.println("Fatal: the vcf file must contain the FORMAT fields Cm and C.");
                        System.exit(1);
                    }
                    baseCallGlobalFieldIndex[i] = vcfParser.getGlobalFieldIndex(sample, "BC");

                    methylStats[i] = backgroundStats.copy();
                    methylStats[i].sampleId = sample;
                    i++;
                }
                int positionGlobalFieldIndex = vcfParser.getGlobalFieldIndex("POS", "VALUE");
                int chromosomeGlobalFieldIndex = vcfParser.getGlobalFieldIndex("CHROM", "VALUE");
                int strandGlobalFieldIndex = vcfParser.getGlobalFieldIndex("INFO", "Strand");


                MethylStats stats = backgroundStats.copy();
                ProgressLogger pg = new ProgressLogger(LOG);
                pg.priority = org.apache.log4j.Level.INFO;
                pg.itemsName = "sites";
                pg.displayFreeMemory = true;
                pg.start();
                final int numSamples = samples.length;


                while (vcfParser.hasNextDataLine()) {
                    pg.lightUpdate();
                    String reference = vcfParser.getFieldValue(chromosomeGlobalFieldIndex).toString();
                    if (ignoreRef(reference)) {
                        vcfParser.next();
                        continue;
                    }
                    int referenceIndex = genome.getReferenceIndex(reference);

                    if (referenceIndex != refIndexOfNextCpG) {
                        // we have change reference sequence or are past the next site.
                        // need to update next CpG position
                        if (referenceIndex > genome.numberOfSequences()) {
                            System.out.printf("alignment reference %s does not exist in genome.", reference);
                            System.exit(1);
                        }
                        referenceSequenceSize = genome.getSequenceSize(referenceIndex);
                        refIndexOfNextCpG = referenceIndex;
                    }
                    // VCF positions are 1-based, but Goby genome positions are 0-based, adjust here:
                    int sitePosition = Integer.parseInt(vcfParser.getFieldValue(positionGlobalFieldIndex).toString()) - 1;
                    char strand = vcfParser.getFieldValue(strandGlobalFieldIndex).charAt(0);
                    if (doFragments) {
                        if (isCpG(referenceIndex, sitePosition, strand)) {

                            final int fragmentLength = calculateFragmentLength(reference, sitePosition, strand);
                            if (fragmentLength > 0) {
                                for (i = 0; i < numSamples; i++) {
                                    final int depthInSample = Integer.parseInt(vcfParser.getFieldValue(sampleDepthGlobalFieldIndex[i]).toString());
                                    if (depthInSample > 10) {
                                        methylStats[i].observedInSample(depthInSample, fragmentLength);
                                    }
                                }
                            }
                        }
                    }

                    updateCpXs(reference, referenceIndex, sitePosition, strand, methylStats, vcfParser, numSamples);

                    vcfParser.next();
                }
                pg.stop();
                pg.done();
//                writeXml(output, samples, methylStats);
                writeTab(output, samples, methylStats);
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

    private void updateCpXs(String reference, int referenceIndex, int sitePosition, char strand,
                            MethylStats[] methylStats, VCFParser vcfParser, int numSamples) {
        if (sitePosition + 1 >= referenceSequenceSize) {
            return;
        }
        if (strand == '-' && sitePosition < 1) { // there is no previous base to check
            return;
        }
        final char firstBase = genome.get(referenceIndex, sitePosition);
        final char secondBase = genome.get(referenceIndex, sitePosition + (strand == '+' ? 1 : -1));
        if (isC(firstBase, strand)) {


            for (int i = 0; i < numSamples; i++) {
                final int depthInSample = Integer.parseInt(vcfParser.getFieldValue(sampleDepthGlobalFieldIndex[i]).toString());
                final CharSequence baseCalls = vcfParser.getFieldValue(baseCallGlobalFieldIndex[i]);
                if ("ignore".equals(baseCalls)) continue;
                final float mr = Integer.parseInt(vcfParser.getFieldValue(methylationRateGlobalFieldIndex[i]).toString());
                final int numCm = Integer.parseInt(vcfParser.getFieldValue(unconvertedCystosineGlobalFieldIndex[i]).toString());
                final int numCConverted = Integer.parseInt(vcfParser.getFieldValue(convertedCystosineGlobalFieldIndex[i]).toString());

                final MethylStats stats = methylStats[i];
                if (base(secondBase, strand) != 'G') {
                    stats.numConvertedNotCpGContext += numCConverted;
                    stats.numNotCpGContext += numCConverted + numCm;
                    /*   System.out.printf("non CpG: pos=%d strand=%c numCm=%d numC=%d mr=%g depth=%d %c%c %n", sitePosition,
             strand,
             numCm, numCConverted, mr,
             depthInSample,
             firstBase, secondBase);      */
                }
                if (depthInSample < 10 || mr < 10) {
                    // discard positions if less than 10 bases observed methylated or less than 10% methylation. We do this to try to avoid
                    // sequencing errors.
                    continue;
                }

                long[] mCpXfreqs = stats.getMethylCpXFreqs();

                switch (base(secondBase, strand)) {
                    case 'C':
                        mCpXfreqs[MethylStats.CPC] += numCm;
                        // final CharSequence baseCalls = vcfParser.getFieldValue(baseCallGlobalFieldIndex[i]);
                        //   System.out.printf("ref: %s position: %d strand %c %c %c baseCalls=%s %n", reference, sitePosition+1, strand, firstBase, secondBase,baseCalls);
                        break;
                    case 'A':
                        mCpXfreqs[MethylStats.CPA] += numCm;

                        break;
                    case 'T':
                        mCpXfreqs[MethylStats.CPT] += numCm;
                        break;
                    case 'G':
                        stats.numCTpG += depthInSample;
                        mCpXfreqs[MethylStats.CPG] += numCm;
                        break;
                }
            }


        }
    }

    private char base(char base, char strand) {
        if (strand == '+') {
            return base;
        } else {
            switch (base) {
                case 'A':
                    return 'T';
                case 'C':
                    return 'G';
                case 'T':
                    return 'A';
                case 'G':
                    return 'C';
                default:
                    return base;
            }
        }

    }

    private boolean isC(final char firstBase, final char strand) {
        return strand == '+' && firstBase == 'C' || strand == '-' && firstBase == 'G';
    }

    private void scanOneStrand(MethylStats backgroundStats, char strand) {
        ProgressLogger pg = new ProgressLogger(LOG);

        pg.priority = org.apache.log4j.Level.INFO;
        pg.itemsName = "bases";
        pg.displayFreeMemory = true;
        pg.start(String.format("counting genome CpG sites on %c strand.", strand));
        for (int sequenceIndex = 0; sequenceIndex < genome.numberOfSequences(); sequenceIndex++) {
            if (QUICK && genome.getReferenceIndex("1") != sequenceIndex) {
                continue;
            }

            final int genomeSequenceSize = genome.getSequenceSize(sequenceIndex);
            char previousBase = '\0';
            char currentBase = '\0';
            int positionCpG1 = 0;
            int positionCpG2 = 0;
            int fragmentLength = -1;
            for (int position = strand == '+' ? 0 :
                    genomeSequenceSize - 1; strand == '+' ? position < genomeSequenceSize :
                    position >= 0; position += strand == '+' ? 1 : -1) {

                currentBase = genome.get(sequenceIndex, position);
                pg.lightUpdate();
                if (isCpG(previousBase, currentBase, strand)) {
                    // found a new CpG.
                    positionCpG1 = positionCpG2;
                    positionCpG2 = position;
                    fragmentLength = strand == '+' ?
                            positionCpG2 - positionCpG1 :
                            positionCpG1 - positionCpG2;
                    if (fragmentLength > 0) {
                        backgroundStats.genomeHasCpG(fragmentLength);

                    }
                }

                previousBase = currentBase;
            }
        }
        pg.done();
    }

    private boolean isCpG(char previousBase, char currentBase, char strand) {
        if (strand == '+') {
            return previousBase == 'C' && currentBase == 'G';
        } else {
            return previousBase == 'G' && currentBase == 'C';
        }
    }

    private void writeTab(PrintWriter output, String[] samples, MethylStats[] methylStats) {
        int sampleIndex = 0;

        if (doFragments) {
            for (String sample : samples) {
                MethylStats methylStat = methylStats[sampleIndex];
                int fragBinIndex = 0;
                output.printf("%s\tnumCpGsObserved\t%d%n", sample, methylStats[sampleIndex].getNumberCpGsObserved());
                for (final int fragmentLengthStart : fragmentLengthBins) {

                    output.printf("%s\tnumberCpGsPerFragmentBinObserved\t%s\t%d%n",
                            sample,
                            getRange(methylStat.getFragmentLengthBins(), fragBinIndex),
                            methylStat.getNumberCpGsPerFragmentBinObserved()[fragBinIndex]);
                    fragBinIndex++;
                }
                fragBinIndex = 0;
                for (final int fragmentLengthStart : fragmentLengthBins) {
                    output.printf("%s\tnormalizedCpGsPerFragmentBinObserved\t%s\t%3.3g%n", sample,
                            getRange(methylStat.getFragmentLengthBins(), fragBinIndex),
                            divide(methylStat.getNumberCpGsPerFragmentBinObserved()[fragBinIndex],
                                    methylStat.getNumberCpGsPerFragmentBinGenome()[fragBinIndex])
                    );
                    fragBinIndex++;
                }
                double sum = 0;
                int depthIndex = 0;
                for (final int depth : depths) {
                    sum += methylStat.getNumberCpGsPerDepth()[depthIndex];
                    depthIndex++;
                }
                depthIndex = 0;
                for (final int depth : depths) {
                    output.printf("%s\tnormalizedDepth\t%s\t%3.3g%n", sample,
                            getRange(methylStat.getDepths(), depthIndex),
                            divide(methylStat.getNumberCpGsPerDepth()[depthIndex],
                                    sum)
                    );
                    depthIndex++;
                }
            }
        }
        for (String sample : samples) {
            MethylStats methylStat = methylStats[sampleIndex];
            double sumCmpX = 0;
            for (int j = MethylStats.CPMIN; j < MethylStats.CPMAX; j++) {
                sumCmpX += methylStat.getMethylCpXFreqs()[j];
            }
            for (int j = MethylStats.CPMIN; j < MethylStats.CPMAX; j++) {

                output.printf("%s\tCPX\t%s\t%3.4g%n", sample,

                        getCpString(j),
                        100.0 * divide(methylStat.getMethylCpXFreqs()[j],
                                sumCmpX)
                );
            }

            output.printf("%s\tconversionRateNonCpGContext\t%d\t%d\t%3.5g%n",
                    sample,
                    methylStat.numConvertedNotCpGContext, methylStat.numNotCpGContext,
                    100.0 * divide(methylStat.numConvertedNotCpGContext,
                            methylStat.numNotCpGContext)
            );

            sampleIndex++;
        }
    }


    private String getCpString(final int cpIndex) {
        switch (cpIndex) {
            case MethylStats.CPA:
                return "CpA";
            case MethylStats.CPG:
                return "CpG";
            case MethylStats.CPC:
                return "CpC";
            case MethylStats.CPT:
                return "CpT";
        }
        return "???";
    }

    private String getRange(final int[] array, final int index) {
        return String.format("[%d-%s[",
                array[index],
                index + 1 < array.length ?
                        Integer.toString(array[index + 1]) :
                        "inf");

    }

    private boolean ignoreRef(String reference) {
        if (!QUICK) {
            return false;
        } else {
            return (!reference.equals("1"));
        }
    }

    private void writeXml(PrintWriter output, String[] samples, MethylStats[] methylStats) throws JAXBException {
        int i;
        final JAXBContext jc = JAXBContext.newInstance("edu.cornell.med.icb.goby.xml");
        final Marshaller m = jc.createMarshaller();
        i = 0;
        for (String sample : samples) {
            m.marshal(methylStats[i], output);
            i++;
        }
    }

    private boolean isCpG(final int referenceIndex, final int sitePosition, char strand) {
        if (sitePosition + 1 > referenceSequenceSize) {
            return false;
        }
        char firstBase = genome.get(referenceIndex, sitePosition);
        char secondBase = genome.get(referenceIndex, sitePosition + 1);
        /*  if (firstBase != 'N' && secondBase != 'N') {
           System.out.printf("%c %c %n", firstBase, secondBase);
       } */
        return isCpG(firstBase, secondBase, strand);

    }

    private int calculateFragmentLength(String reference, int sitePosition, final char strand) {
        int referenceIndex = genome.getReferenceIndex(reference);

        char previousBase = '\0';
        char currentBase = '\0';
        int positionCpG1 = sitePosition;
        int positionCpG2 = 0;
        int fragmentLength = -1;
        int direction = strand == '+' ? 1 : -1;
        for (int pos = sitePosition + 2 * direction; strand == '+' ? pos < referenceSequenceSize : pos >= 0;
             pos += direction) {
            currentBase = genome.get(referenceIndex, pos);
            if (isCpG(previousBase, currentBase, strand)) {
                // found a new CpG.
                positionCpG1 = sitePosition;
                positionCpG2 = pos;

                fragmentLength = strand == '+' ? positionCpG2 - positionCpG1 : positionCpG1 - positionCpG2;

                if (fragmentLength > 0) {
                    return fragmentLength;
                } else return -1;
            }
            previousBase = currentBase;
        }
        return -1;
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