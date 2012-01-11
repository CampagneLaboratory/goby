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
import edu.cornell.med.icb.goby.stats.SampleStats;
import edu.cornell.med.icb.goby.util.GrepReader;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.*;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;


/**
 * Compare two VCF files. In contrast to compare-vcf (vcf-tools), this mode
 * <LI/>Can handle non diploid genotype calls.
 * <LI/>Can compare specific pairs of samples identified by their name (or a keyword matching one sample name) on
 * the command line.
 *
 * @author Fabien Campagne
 */
public class VCFCompareMode extends AbstractGobyMode {
    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(VCFCompareMode.class);

    /**
     * The input files.
     */
    private List<File> inputFiles;

    /**
     * The output filename.
     */
    private String outputFilename;


    /**
     * The mode name.
     */
    private static final String MODE_NAME = "vcf-compare";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Compare genotype calls across VCF files.";

    private int[] chromosomeFieldIndex;
    private int[] positionFieldIndex;

    private String[] genotypeColumnSet;
    private ObjectArrayList<VCFLine>[] lines;
    private Object2IntMap<VCFPosition>[] indices;
    private int[] refFieldIndex;
    private int[] altFieldIndex;


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
    public AbstractCommandLineMode configure(final String[] args)
            throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);
        setInputFilenames(jsapResult.getStringArray("input"));
        outputFilename = jsapResult.getString("output");
        this.genotypeColumnSet = jsapResult.getStringArray("column");

        return this;
    }

    /**
     * Add an input file.
     *
     * @param inputFile the input file to add.
     */
    public synchronized void addInputFile(final File inputFile) {
        if (inputFiles == null) {
            inputFiles = new LinkedList<File>();
        }
        this.inputFiles.add(inputFile);
    }

    /**
     * Clear the input files list.
     */
    public synchronized void clearInputFiles() {
        if (inputFiles != null) {
            inputFiles.clear();
        }
    }

    /**
     * Set the input filenames.
     *
     * @param inputFilenames the input filename
     */
    public synchronized void setInputFilenames(final String[] inputFilenames) {
        clearInputFiles();
        for (final String inputFilname : inputFilenames) {
            addInputFile(new File(inputFilname));
        }
    }

    /**
     * Get the input filenames.
     *
     * @return the input filenames
     */
    public synchronized String[] getInputFilenames() {
        if (inputFiles == null) {
            return new String[0];
        }
        final String[] array = new String[inputFiles.size()];
        int i = 0;
        for (final File inputFile : inputFiles) {
            array[i++] = inputFile.toString();
        }
        return array;
    }

    Object2IntMap<String> sampleIndexToGenotypeIndex;

    /**
     * Compare VCF files.
     *
     * @throws java.io.IOException
     */
    @Override
    public void execute() throws IOException {
        if (inputFiles == null || inputFiles.size() == 0) {
            throw new IOException("--input not specified");
        }
        /*  if (StringUtils.isBlank(outputFilename)) {
            throw new IOException("--output not specified");
        }*/
        final int numInputFiles = inputFiles.size();
        final VCFParser[] parsers = new VCFParser[numInputFiles];
        int parserIndex = 0;
        chromosomeFieldIndex = new int[numInputFiles];
        positionFieldIndex = new int[numInputFiles];
        refFieldIndex = new int[numInputFiles];
        altFieldIndex = new int[numInputFiles];

        scanInput(parsers, parserIndex);
        // brute force, load everything in memory:
        loadGenotypes(numInputFiles, parsers);

        System.out.println("Determining common positions");
        ObjectSet<VCFPosition> commonPositions = new ObjectAVLTreeSet<VCFPosition>();

        commonPositions.addAll(reduce(lines[0]));
        System.out.printf("# position parser[0] %d%n", commonPositions.size());

        for (parserIndex = 1; parserIndex < numInputFiles; parserIndex++) {

            commonPositions.retainAll(reduce(lines[parserIndex]));
        }
        DoubleIndexedIdentifier reverseIdentifiers = new DoubleIndexedIdentifier(identifiers);
        for (int i = 0; i < numInputFiles; i++) {
            dumpPositionsUniqueTo(i, inputFiles, lines, commonPositions, reverseIdentifiers);
        }
        System.out.println("Sorting..");
        ObjectArrayList<VCFPosition> sortedPositions = new ObjectArrayList<VCFPosition>();
        sortedPositions.addAll(commonPositions);
        Collections.sort(sortedPositions);
        System.out.println("Done sorting.");

        printStats(numInputFiles, commonPositions, sortedPositions, reverseIdentifiers);
        System.out.flush();
        System.exit(0);
    }

    private void loadGenotypes(int numInputFiles, VCFParser[] parsers) throws IOException {
        int parserIndex;
        final IntArrayList[] keepGlobalGenotypeIndex = new IntArrayList[numInputFiles];
        sampleIndexToGenotypeIndex = new Object2IntAVLTreeMap<String>();
        lines = new ObjectArrayList[numInputFiles];
        indices = new Object2IntMap[numInputFiles];
        IntArrayList[] allGenotypeFieldIndices = new IntArrayList[numInputFiles];
        for (parserIndex = 0; parserIndex < numInputFiles; parserIndex++) {
            keepGlobalGenotypeIndex[parserIndex] = new IntArrayList(numInputFiles);
            lines[parserIndex] = new ObjectArrayList<VCFLine>(10000000);
            indices[parserIndex] = new Object2IntAVLTreeMap<VCFPosition>();
            allGenotypeFieldIndices[parserIndex] = new IntArrayList(parsers[parserIndex].getColumnNamesUsingFormat().length);
            for (String sample : genotypeColumnSet) {
                allGenotypeFieldIndices[parserIndex].add(parsers[parserIndex].getGlobalFieldIndex(getSampleColumn(parsers[parserIndex], sample), "GT"));
            }
            sampleIndexToGenotypeIndex = new Object2IntArrayMap<String>();
            sampleIndexToGenotypeIndex.defaultReturnValue(-1);

            for (String keepColumn : genotypeColumnSet) {

                String sampleId = getSampleColumn(parsers[parserIndex], keepColumn);
                final int index = parsers[parserIndex].getGlobalFieldIndex(sampleId, "GT");
                if (index != -1) {
                    keepGlobalGenotypeIndex[parserIndex].add(index);

                    sampleIndexToGenotypeIndex.put(sampleId, Arrays.binarySearch(genotypeColumnSet, keepColumn));
                }

            }
            int[] indicesToKeep = keepGlobalGenotypeIndex[parserIndex].toIntArray();
            System.out.printf("Loading %s..%n", inputFiles.get(parserIndex).getName());
            int index = 0;
            ProgressLogger pg = new ProgressLogger(LOG);
            pg.displayFreeMemory = true;
            pg.priority = Level.INFO;
            pg.start();
            try {

                int earlyStopCount = Integer.MAX_VALUE; //30000;
                //int earlyStopCount = 300;
                int count = 0;
                while (parsers[parserIndex].hasNextDataLine()) {
                    final VCFLine line = new VCFLine();

                    final String chr = parsers[parserIndex].getFieldValue(chromosomeFieldIndex[parserIndex]).toString();
                    String ref = parsers[parserIndex].getStringFieldValue(refFieldIndex[parserIndex]);
                    String alts = parsers[parserIndex].getStringFieldValue(altFieldIndex[parserIndex]);
                    line.pos.chromosome = identifiers.registerIdentifier(new MutableString(chr));
                    line.pos.position = Integer.parseInt(parsers[parserIndex].
                            getFieldValue(positionFieldIndex[parserIndex]).toString());
                    line.ref=ref;
                    // keep this line since there is a variant somewhere on it.
                    int sampleIndex = 0;
                    for (final int fieldIndex : indicesToKeep) {
                        final String genotypeCode = parsers[parserIndex].getFieldValue(fieldIndex).toString();
                        if (genotypeCode == null) {
                            // skip position, not typed in some sample.
                            parsers[parserIndex].next();
                            pg.lightUpdate();
                            continue;
                        }
                        final String genotype = convertCode(genotypeCode, parsers[parserIndex], ref, alts);

                        line.genotypes.add(genotype);

                    }
                    lines[parserIndex].add(line);
                    indices[parserIndex].put(line.pos, index++);

                    parsers[parserIndex].next();
                    pg.lightUpdate();
                    /* if (count++ > earlyStopCount) {
                       break;
                   } */
                }
            } catch (Exception e) {
                e.printStackTrace();
            } finally {
                pg.done();
                parsers[parserIndex].close();
            }
            System.out.println("Line.size() " + lines[parserIndex].size());
        }
    }

    private void scanInput(VCFParser[] parsers, int parserIndex) throws IOException {
        for (final File inputFile : inputFiles) {
            // eliminate lines from the header that try to define some field in ALT:
            final GrepReader filter = new GrepReader(inputFile.getPath(), "^##ALT=");

            parsers[parserIndex] = new VCFParser(filter);
            try {
                parsers[parserIndex].readHeader();
            } catch (VCFParser.SyntaxException e) {
                e.printStackTrace();
                System.exit(1);
            }

            if (genotypeColumnSet.length == 0) {
                // if no columns selected, used the list of columns from the very first command line input file:
                genotypeColumnSet = parsers[parserIndex].getColumnNamesUsingFormat();
                System.err.println("Using default columns from first file: " + ObjectArrayList.wrap(genotypeColumnSet));
            }
            chromosomeFieldIndex[parserIndex] = parsers[parserIndex].getGlobalFieldIndex("CHROM", "VALUE");
            positionFieldIndex[parserIndex] = parsers[parserIndex].getGlobalFieldIndex("POS", "VALUE");
            refFieldIndex[parserIndex] = parsers[parserIndex].getGlobalFieldIndex("REF", "VALUE");
            altFieldIndex[parserIndex] = parsers[parserIndex].getGlobalFieldIndex("ALT", "VALUE");
            assert refFieldIndex[parserIndex] != -1 : "REF field must be found.";
            assert altFieldIndex[parserIndex] != -1 : "ALT field must be found.";
            parserIndex++;
        }
    }

    private void printStats(int numInputFiles,
                            ObjectSet<VCFPosition> commonPositions,
                            ObjectArrayList<VCFPosition> sortedPositions, DoubleIndexedIdentifier reverseIdentifiers) throws IOException {
        final VCFLine alignedLines[] = new VCFLine[numInputFiles];
        int parserIndex = 0;

        ObjectSet<String> distinctGenotypes = new ObjectArraySet(numInputFiles);
        ObjectList<String> sampleGenotypes = new ObjectArrayList(numInputFiles);

        SampleStats[] sampleStats = new SampleStats[genotypeColumnSet.length];
        for (VCFPosition pos : sortedPositions) {
            for (parserIndex = 0; parserIndex < numInputFiles; parserIndex++) {
                alignedLines[parserIndex] = lines[parserIndex].get(indices[parserIndex].get(pos));
                alignedLines[parserIndex].ref=lines[parserIndex].get(indices[parserIndex].get(pos)).ref;
            }

            int sampleIndex = 0;
            for (String sample : genotypeColumnSet) {
                SampleStats sampleStat = sampleStats[sampleIndex];
                if (sampleStat == null) {
                    sampleStat = new edu.cornell.med.icb.goby.stats.SampleStats(numInputFiles);
                    sampleStats[sampleIndex] = sampleStat;
                    sampleStat.sampleId = sample;
                }
                distinctGenotypes.clear();
                sampleGenotypes.clear();
                String ref = null;
                for (final VCFLine line : alignedLines) {
                    final String genotype = line.genotypes.get(sampleIndex);
                    ref=line.ref;
                    distinctGenotypes.add(genotype);
                    sampleGenotypes.add(genotype);
                }
                for (int fileIndex = 0; fileIndex < numInputFiles; fileIndex++) {
                    sampleStat.observeTransitionToTransversions(fileIndex, sampleGenotypes, ref);
                }
                if (distinctGenotypes.size() > 1) {
                    sampleStat.numGenotypeDisagreements++;

                    if (distinctGenotypes.contains("")) {
                        sampleStat.numGenotypeNotCalled++;
                    } else {
                        /*System.out.printf("%s\t%d\tdisagreement: %s %n",
                                reverseIdentifiers.getId(alignedLines[0].pos.chromosome).toString(),
                                alignedLines[0].pos.position,
                                ObjectArrayList.wrap(distinctGenotypes.toArray()));
                          */
                        sampleStat.analyze(distinctGenotypes, sampleGenotypes);
                    }
                } else {
                    sampleStat.numGenotypeAgreements++;
                }
                sampleIndex++;
            }
        }
        System.out.printf("# common positions across files: %d (overlap with larger set: %g %%) " +
                "(overlap with smaller set: %g%%) %n", commonPositions.size(),
                fraction(commonPositions.size(), maxSize(lines)), fraction(commonPositions.size(), minSize(lines)));
        for (final SampleStats sampleStat : sampleStats) {
            System.out.println("Sample: " + sampleStat.sampleId);
            final int sumErrors = sampleStat.numGenotypeNotCalled + sampleStat.numGenotypeDisagreements;
            System.out.printf("Among the common positions, %d positions (%g %%) had the same genotype, while %d positions (%g %%) had some disagreements (failure to call a genotype in other method, failure to call one or more alleles, or different genotype called: hard error). \n" +
                    "Among the differences, %g %% were failures to call any genotype %g %% were failures to call one allele, %g %% to call two, and %g %% to call more than two. %g %% sites had differences in genotypes that could not be explained by a failure to call an allele (e.g., G/G vs G/T when the reference is A/A)%n",
                    sampleStat.numGenotypeAgreements, fractionCumul(sampleStat.numGenotypeAgreements, sampleStat.numGenotypeDisagreements),
                    sampleStat.numGenotypeDisagreements, fractionCumul(sampleStat.numGenotypeDisagreements, sampleStat.numGenotypeAgreements),

                    fractionCumul(sampleStat.numGenotypeNotCalled, sumErrors),
                    fractionCumul(sampleStat.missedOneAlleles, sumErrors),
                    fractionCumul(sampleStat.missedTwoAlleles, sumErrors),
                    fractionCumul(sampleStat.missedMoreThanTwoAlleles, sumErrors),
                    fractionCumul(sampleStat.numHadDifferentAllele, sumErrors)
            );
        }
        if (outputFilename != null) {
            PrintWriter out = new PrintWriter(new FileWriter(outputFilename));
            out.write(SampleStats.header(numInputFiles));
            for (final SampleStats sampleStat : sampleStats) {

                out.write(sampleStat.toString());
            }
            out.close();
        }


    }

    private void dumpPositionsUniqueTo(int i, List<File> inputFiles, ObjectArrayList<VCFLine>[] lines, ObjectSet<VCFPosition> commonPositions, DoubleIndexedIdentifier reverseIdentifiers) {
        /*
   System.out.printf("Unique to %s ---------:%n", inputFiles.get(i).getName());
   for (VCFLine line : lines[i]) {
       if (!commonPositions.contains(line.pos)) {
           System.out.printf("%s GT=%s %n", line.pos.toString(reverseIdentifiers), line.genotypes);

       }
   }
   System.out.printf("<---------%n");
        */
    }


    MutableString buffer = new MutableString();
    ObjectArrayList<String> list = new ObjectArrayList<String>();

    private String convertCode(String genotypeCode, VCFParser parser, String ref, String alts) {
        String[] tokens = genotypeCode.split("[/|]");
        buffer.setLength(0);
        String[] altArray = alts.split(",");
        list.clear();
        for (String token : tokens) {
            try {
                if (".".equals(token)) {
                    return "";
                }
                int index = Integer.parseInt(token);
                if (index == 0) {

                    list.add("ref");

                } else {
                    list.add(altArray[index - 1]);

                }
            } catch (NumberFormatException e) {
                LOG.info("genotype could not be parsed: " + genotypeCode);
                return "genotype-error:" + genotypeCode;
            }


        }
        Collections.sort(list);
        for (String allele : list) {
            buffer.append(allele);
            buffer.append('/');
        }
        return buffer.toString();
    }

    private String getSampleColumn(VCFParser parser, String keepColumn) {
        String[] samples = parser.getColumnNamesUsingFormat();
        ObjectArraySet<String> set = new ObjectArraySet<String>();
        for (String sample : samples) {
            if (sample.contains(keepColumn)) {
                set.add(sample);
            }
        }
        if (set.size() == 1) {
            return set.iterator().next();
        } else {
            System.err.println("Several columns match the sample keyword provided: " + keepColumn + " " + set);
            throw new RuntimeException("Several columns match the sample keyword provided: " + set);
        }
    }

    private double maxSize(final ObjectArrayList<VCFLine>[] lines) {
        int max = 0;
        for (final ObjectArrayList<VCFLine> line : lines) {
            max = Math.max(line.size(), max);
        }
        return max;
    }

    private double minSize(final ObjectArrayList<VCFLine>[] lines) {
        int min = Integer.MAX_VALUE;
        for (final ObjectArrayList<VCFLine> line : lines) {
            min = Math.min(line.size(), min);
        }
        return min;
    }

    private double fraction(double a, double b) {
        return (a / (b)) * 100;
    }

    private double fractionCumul(double a, double b) {
        return (a / (a + b)) * 100;
    }

    private ObjectOpenHashSet<VCFPosition> reduce(ObjectArrayList<VCFLine> lines) {
        ObjectOpenHashSet<VCFPosition> positions = new ObjectOpenHashSet<VCFPosition>();
        for (VCFLine line : lines) {
            positions.add(line.pos);

        }
        return positions;
    }


    /**
     * Main mode for splitting compact reads files from a start position
     * to and end position.
     *
     * @param args command line arguments
     * @throws java.io.IOException IO error
     * @throws com.martiansoftware.jsap.JSAPException
     *                             command line parsing error.
     */
    public static void main(final String[] args) throws IOException, JSAPException {
        new VCFCompareMode().configure(args).execute();
    }

    private class VCFLine {
        VCFPosition
                pos = new VCFPosition();
        String ref;
        String alt;
        /**
         * One genotype per sample kept:
         */
        ObjectArrayList<String> genotypes = new ObjectArrayList<String>(2);
    }

    static IndexedIdentifier identifiers = new IndexedIdentifier();

    private class VCFPosition implements Comparable {

        int chromosome;
        int position;

        @Override
        public String toString() {
            return String.format("%s\t%d", chromosome, position);
        }

        public String toString(DoubleIndexedIdentifier reverse) {
            return String.format("%s\t%d", reverse.getId(chromosome), position);
        }

        @Override
        public boolean equals(Object o) {
            VCFPosition other = (VCFPosition) o;
            return position == other.position && chromosome == other.chromosome;
        }

        @Override
        public int hashCode() {
            int hash = 1;
            hash = hash * 31 + position;
            hash = hash * 31 + chromosome;
            return hash;
        }

        public int compareTo(Object o) {
            VCFPosition other = (VCFPosition) o;
            int result = chromosome - other.chromosome;
            if (result != 0) return result;
            else {
                return position - other.position;
            }
        }
    }

}