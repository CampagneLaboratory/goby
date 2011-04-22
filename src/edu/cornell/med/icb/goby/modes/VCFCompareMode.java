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
import edu.cornell.med.icb.goby.reads.Reads;
import edu.cornell.med.icb.goby.reads.ReadsReader;
import edu.cornell.med.icb.goby.reads.ReadsWriter;
import edu.cornell.med.icb.goby.reads.ReadSet;
import edu.cornell.med.icb.goby.reads.MessageChunksWriter;
import edu.cornell.med.icb.goby.readers.vcf.VCFParser;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import edu.rit.numeric.ListXYSeries;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.logging.ProgressLogger;
import it.unimi.dsi.fastutil.objects.*;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.chars.Char2ShortOpenHashMap;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.FileInputStream;
import java.util.LinkedList;
import java.util.List;
import java.util.Collections;
import java.nio.channels.FileChannel;

/**
 * Compare two VCF files. In contrast to compare-vcf (vcf-tools), this mode
 * <LI/>Can handle non diployd genotype calls.
 * <LI/>Can compare specific pairs of samples identified by their name on the command line.
 *
 * @author Fabien Campagne
 */
public class VCFCompareMode extends AbstractGobyMode {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(VCFCompareMode.class);

    /**
     * The input files.
     */
    private List<File> inputFiles;

    /**
     * The output filename.
     */
    private String outputFilename;

    /**
     * sequences per chunck in the written file.
     */
    private int sequencePerChunk = 10000;

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
    private Object2IntMap<VCFPosition> indices = new Object2IntArrayMap<VCFPosition>();


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
        for (File inputFile : inputFiles) {

            parsers[parserIndex] = new VCFParser(inputFile.getPath());
            try {
                parsers[parserIndex].readHeader();
            } catch (VCFParser.SyntaxException e) {
                e.printStackTrace();
                System.exit(1);
            }
            chromosomeFieldIndex[parserIndex] = parsers[parserIndex].getGlobalFieldIndex("CHROM", "VALUE");
            positionFieldIndex[parserIndex] = parsers[parserIndex].getGlobalFieldIndex("POS", "VALUE");
            parserIndex++;
        }
        // brute force, load everything in memory:
        IntSet[] keepGlobalGenotypeIndex = new IntSet[numInputFiles];
        lines = new ObjectArrayList[numInputFiles];
        for (parserIndex = 0; parserIndex < numInputFiles; parserIndex++) {
            keepGlobalGenotypeIndex[parserIndex] = new IntArraySet();
            lines[parserIndex] = new ObjectArrayList<VCFLine>();
            for (String keepColumn : genotypeColumnSet) {

                final int index = parsers[parserIndex].getGlobalFieldIndex(keepColumn, "GT");
                if (index != -1) {
                    keepGlobalGenotypeIndex[parserIndex].add(index);
                }

            }

            System.out.printf("Loading %s..%n", inputFiles.get(parserIndex).getName());
            int index = 0;
            while (parsers[parserIndex].hasNextDataLine()) {
                VCFLine line = new VCFLine();

                final String chr = parsers[parserIndex].getFieldValue(chromosomeFieldIndex[parserIndex]).toString();
                line.pos.chromosome = identifiers.registerIdentifier(new MutableString(chr));
                line.pos.position = Integer.parseInt(parsers[parserIndex].
                        getFieldValue(positionFieldIndex[parserIndex]).toString());
                for (int fieldIndex : keepGlobalGenotypeIndex[parserIndex]) {
                    line.genotypes.add(parsers[parserIndex].getFieldValue(fieldIndex).toString());
                }
                lines[parserIndex].add(line);
                indices.put(line.pos, index++);
                parsers[parserIndex].next();
            }
            parsers[parserIndex].close();
            System.out.println("Line.size()" + lines[parserIndex].size());
        }
        ObjectOpenHashSet<VCFPosition> commonPositions = new ObjectOpenHashSet<VCFPosition>();

        commonPositions.addAll(reduce(lines[0]));
        System.out.printf("# position parser[0] %d%n", commonPositions.size());


        for (parserIndex = 1; parserIndex < numInputFiles; parserIndex++) {

            commonPositions.retainAll(reduce(lines[parserIndex]));
        }
        System.out.printf("# common positions across files: %d%n", commonPositions.size());
        ObjectArrayList<VCFPosition> sortedPositions = new ObjectArrayList<VCFPosition>();
        sortedPositions.addAll(commonPositions);
        Collections.sort(sortedPositions);
        VCFLine alignedLines[] = new VCFLine[numInputFiles];
        parserIndex = 0;
        for (VCFPosition pos : sortedPositions) {
            for (parserIndex = 0; parserIndex < numInputFiles; parserIndex++) {
                alignedLines[parserIndex] = lines[parserIndex].get(indices.get(pos));
            }
            
        }
        System.exit(0);
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
        ObjectArrayList<String> genotypes = new ObjectArrayList<String>();
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