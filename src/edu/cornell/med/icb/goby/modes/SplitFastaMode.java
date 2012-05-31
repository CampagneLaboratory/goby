/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
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
import edu.cornell.med.icb.goby.readers.FastXEntry;
import edu.cornell.med.icb.goby.readers.FastXReader;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import org.apache.commons.io.FilenameUtils;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Split a FASTA/FASTQ file into multiple parts. Primarily used for splitting FASTA
 * into multiple chunks for Eland with a maximum number of entries per file and
 * one read length per file.
 * @author Kevin Dorff
 */
public class SplitFastaMode extends AbstractGobyMode {

    /** The mode name. */
    private static final String MODE_NAME = "split-fasta";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Split a FASTA/FASTQ file into multiple parts. "
            + "Primarily used for splitting FASTA into multiple chunks for Eland with a maximum "
            + "number of entries per file and one read length per file.";


    /** Default value. */
    private static final int DEFAULT_SPLIT_MAX_LENGTH = 32;

    /**
     * The number of entries to write to a single file in splitfastx, this is
     * needed due to a limitation of Eland.
     */
    private static final int MAX_READS_PER_FILE_DEFAULT = 16777000;

    /** The input file. */
    private String inputFile;

    /** Max read length for reading / writing files. */
    private int fastxSplitMaxLength;

    /** The number of read lengths for a file, if 1 each read length will get it's own file.
     * If 50, a file will be generated for read lengths 1-50, another for 51-100, etc. */
    private int splitReadsMod;

    private int maxReadsPerFile;

    /** Map to override help / default values. */
    private static final Map<String, String> HELP_VALUES;
    static {
        HELP_VALUES = new HashMap<String, String>();
        HELP_VALUES.put("[SPLIT_MAX_LENGTH]",
                Integer.toString(DEFAULT_SPLIT_MAX_LENGTH));
        HELP_VALUES.put("[MAX_READS_PER_FILE_DEFAULT]",
                Integer.toString(MAX_READS_PER_FILE_DEFAULT));
    }

    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }

    /**
     * Configure the mode arguements.
     * @param args the arguments
     * @return this object for chaining
     * @throws IOException error parsing arguments
     * @throws JSAPException error parsing arguments
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args) throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args, HELP_VALUES);

        inputFile = new File(jsapResult.getString("input")).getCanonicalPath();
        fastxSplitMaxLength = jsapResult.getInt("split-max-length");
        splitReadsMod = jsapResult.getInt("split-reads-mod");
        maxReadsPerFile = jsapResult.getInt("max-reads-per-file");
        return this;
    }

    /**
     * Split a fasta / fastq file by (a) readlength and (b) the maximum number of
     * entries per file. This will output the files that are written to stdout
     * @throws IOException error reading / writing files.
     */
    @Override
    public void execute() throws IOException {
        final FastXReader reader = new FastXReader(inputFile);
        final Int2ObjectMap<PrintStream> outputMap = new Int2ObjectOpenHashMap<PrintStream>();
        final Int2IntMap entriesPerReadLen = new Int2IntOpenHashMap();
        final Int2IntMap filesPerReadLen = new Int2IntOpenHashMap();
        final List<String> removeExt = Arrays.asList(
                "gz", "fa", "mpfa", "fna", "fsa", "fas", "fasta",
                "fq", "mpfq", "fnq", "fsq", "fas", "fastq");
        String inputName = FilenameUtils.getName(inputFile);
        while (true) {
            // Remove the unwanted extensions from the file name
            final String ext = FilenameUtils.getExtension(inputName);
            if (!removeExt.contains(ext)) {
                break;
            }
            inputName = FilenameUtils.getBaseName(inputName);
        }
        final String outputFilenameTemplate =
                FilenameUtils.getFullPath(inputFile) + inputName
                + "._READLENGTH_._PART_." + reader.getFileType();
        final NumberFormat nf3 = NumberFormat.getInstance();
        nf3.setMinimumIntegerDigits(3);
        final NumberFormat nf2 = NumberFormat.getInstance();
        nf2.setMinimumIntegerDigits(2);
        for (final FastXEntry entry : reader) {
            final int readLen = Math.min(fastxSplitMaxLength,
                    roundReadLen(entry.getReadLength(), splitReadsMod));
            PrintStream out = outputMap.get(readLen);
            if (out == null) {
                filesPerReadLen.put(readLen, 1);
                entriesPerReadLen.put(readLen, 0);
                String outputFilename = outputFilenameTemplate.replaceAll(
                        "_READLENGTH_", nf3.format(readLen));
                outputFilename = outputFilename.replaceAll(
                        "_PART_", nf2.format(1));
                System.out.println(outputFilename);
                out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFilename)));
                outputMap.put(readLen, out);
            }
            int numEntries = entriesPerReadLen.get(readLen);
            if (numEntries == maxReadsPerFile) {
                out.close();
                numEntries = 0;
                int numFiles = filesPerReadLen.get(readLen);
                numFiles++;
                filesPerReadLen.put(readLen, numFiles);
                String outputFilename = outputFilenameTemplate.replaceAll(
                        "_READLENGTH_", nf3.format(readLen));
                outputFilename = outputFilename.replaceAll(
                        "_PART_", nf2.format(numFiles));
                System.out.println(outputFilename);
                out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFilename)));
                outputMap.put(readLen, out);
            }
            out.println(entry.getEntry());
            entriesPerReadLen.put(readLen, numEntries + 1);
        }
        for (final PrintStream out : outputMap.values()) {
            out.close();
        }
        outputMap.clear();
        reader.close();
    }

    /**
     * Round the readLen based on the value in mod. If mod==1
     * then this will return i. Otherwise it will use mod to
     * round up. If mod==50, values 1..50 will return 50,
     * values 51..100 will return 100, etc.
     * @param i the readLen to round
     * @param mod the mod to use
     * @return the rounded readLen
     */
    public int roundReadLen(final int i, final int mod) {
        if (mod == 1) {
            return i;
        }
        final int j = i / mod;
        final int k = i % mod;
        if (k == 0) {
            return (j * mod);
        } else {
            return ((j + 1) * mod);
        }
    }

    /**
     * Main method.
     * @param args command line args.
     * @throws JSAPException error parsing
     * @throws IOException error parsing or executing.
     */
    public static void main(final String[] args) throws JSAPException, IOException {
        new SplitFastaMode().configure(args).execute();
    }
}
