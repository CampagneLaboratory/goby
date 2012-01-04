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
import edu.cornell.med.icb.goby.util.GrepReader;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.*;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;
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
        if (StringUtils.isBlank(outputFilename)) {
            throw new IOException("--output not specified");
        }
        final int numInputFiles = inputFiles.size();
        VCFParser parsers[] = new VCFParser[numInputFiles];
        int parserIndex = 0;
        chromosomeFieldIndex = new int[numInputFiles];
        positionFieldIndex = new int[numInputFiles];
        refFieldIndex = new int[numInputFiles];
        altFieldIndex = new int[numInputFiles];
        for (File inputFile : inputFiles) {
            // eliminate lines from the header that try to define some field in ALT:
            final GrepReader filter = new GrepReader(inputFile.getPath(), "^##ALT=");

            parsers[parserIndex] = new VCFParser(filter);
            try {
                parsers[parserIndex].readHeader();
            } catch (VCFParser.SyntaxException e) {
                e.printStackTrace();
                System.exit(1);
            }
            chromosomeFieldIndex[parserIndex] = parsers[parserIndex].getGlobalFieldIndex("CHROM", "VALUE");
            positionFieldIndex[parserIndex] = parsers[parserIndex].getGlobalFieldIndex("POS", "VALUE");
            refFieldIndex[parserIndex] = parsers[parserIndex].getGlobalFieldIndex("REF", "VALUE");
            altFieldIndex[parserIndex] = parsers[parserIndex].getGlobalFieldIndex("ALT", "VALUE");
            assert refFieldIndex[parserIndex] != -1 : "REF field must be found.";
            assert altFieldIndex[parserIndex] != -1 : "ALT field must be found.";
            parserIndex++;
        }
        // brute force, load everything in memory:
        final IntSet[] keepGlobalGenotypeIndex = new IntSet[numInputFiles];
        lines = new ObjectArrayList[numInputFiles];
        indices = new Object2IntMap[numInputFiles];
        for (parserIndex = 0; parserIndex < numInputFiles; parserIndex++) {
            keepGlobalGenotypeIndex[parserIndex] = new IntArraySet();
            lines[parserIndex] = new ObjectArrayList<VCFLine>(10000000);
            indices[parserIndex] = new Object2IntAVLTreeMap<VCFPosition>();

            for (String keepColumn : genotypeColumnSet) {

                final int index = parsers[parserIndex].getGlobalFieldIndex(getSampleColumn(parsers[parserIndex], keepColumn), "GT");
                if (index != -1) {
                    keepGlobalGenotypeIndex[parserIndex].add(index);
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
                // load at most three million lines:
                int earlyStopCount = Integer.MAX_VALUE; //30000;
                int count = 0;
                while (parsers[parserIndex].hasNextDataLine()) {
                    final VCFLine line = new VCFLine();

                    final String chr = parsers[parserIndex].getFieldValue(chromosomeFieldIndex[parserIndex]).toString();
                    String ref = parsers[parserIndex].getStringFieldValue(refFieldIndex[parserIndex]);
                    String alts = parsers[parserIndex].getStringFieldValue(altFieldIndex[parserIndex]);
                    line.pos.chromosome = identifiers.registerIdentifier(new MutableString(chr));
                    line.pos.position = Integer.parseInt(parsers[parserIndex].
                            getFieldValue(positionFieldIndex[parserIndex]).toString());
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
                    /*  if (count++ > earlyStopCount) {
                        break;
                    }*/
                }
            } catch (Exception e) {
                e.printStackTrace();
            } finally {
                pg.done();
                parsers[parserIndex].close();
            }
            System.out.println("Line.size() " + lines[parserIndex].size());
        }
        System.out.println("Determining common positions");
        ObjectSet<VCFPosition> commonPositions = new ObjectAVLTreeSet<VCFPosition>();

        commonPositions.addAll(reduce(lines[0]));
        System.out.printf("# position parser[0] %d%n", commonPositions.size());

        for (parserIndex = 1; parserIndex < numInputFiles; parserIndex++) {

            commonPositions.retainAll(reduce(lines[parserIndex]));
        }
        System.out.printf("# common positions across files: %d (%g %%)%n", commonPositions.size(),
                fraction(commonPositions.size(), maxSize(lines)));
        System.out.println("Sorting..");
        ObjectArrayList<VCFPosition> sortedPositions = new ObjectArrayList<VCFPosition>();
        sortedPositions.addAll(commonPositions);
        Collections.sort(sortedPositions);
        System.out.println("Done sorting.");
        final VCFLine alignedLines[] = new VCFLine[numInputFiles];
        parserIndex = 0;
        int numGenotypeAgreements = 0;
        int numGenotypeDisagreements = 0;
        ObjectSet<String> distinctGenotypes = new ObjectArraySet(numInputFiles);

        for (VCFPosition pos : sortedPositions) {
            for (parserIndex = 0; parserIndex < numInputFiles; parserIndex++) {
                alignedLines[parserIndex] = lines[parserIndex].get(indices[parserIndex].get(pos));
            }
            distinctGenotypes.clear();
            for (final VCFLine line : alignedLines) {
                distinctGenotypes.addAll(line.genotypes);
            }
            if (distinctGenotypes.size() > 1) {
                numGenotypeDisagreements++;
            } else {
                numGenotypeAgreements++;
            }
        }


        System.out.printf("Among the common positions, %d positions (%g %%) had the same genotype, while %d positions (%g %%) had different genotypes.",
                numGenotypeAgreements, fraction(numGenotypeAgreements, numGenotypeDisagreements),
                numGenotypeDisagreements, fraction(numGenotypeDisagreements, numGenotypeAgreements));
        System.exit(0);
    }

    MutableString buffer = new MutableString();

    private String convertCode(String genotypeCode, VCFParser parser, String ref, String alts) {
        String[] tokens = genotypeCode.split("[|/]");
        buffer.setLength(0);
        String[] altArray = alts.split(",");
        for (String token : tokens) {
            try {
                if (".".equals(token)) {
                    return "";
                }
                int index = Integer.parseInt(token);
                if (index == 0) {

                    buffer.append(ref);
                    buffer.append(',');
                } else {

                    buffer.append(altArray[index - 1]);
                    buffer.append(',');
                }
            } catch (NumberFormatException e) {
                LOG.info("genotype could not be parsed: " + genotypeCode);
                return "genotype-error";
            }
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

    private double fraction(double a, double b) {
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
        VCFPosition pos = new VCFPosition();
        String ref;
        String alt;
        ObjectArrayList<String> genotypes = new ObjectArrayList<String>(4);
    }

    static IndexedIdentifier identifiers = new IndexedIdentifier();

    private class VCFPosition implements Comparable {

        int chromosome;
        int position;

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