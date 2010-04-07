/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                         Weill Medical College of Cornell University
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
import edu.cornell.med.icb.goby.reads.QualityEncoding;
import edu.cornell.med.icb.goby.reads.ReadsWriter;
import edu.cornell.med.icb.goby.util.DoInParallel;
import edu.cornell.med.icb.goby.util.FileExtensionHelper;
import it.unimi.dsi.fastutil.io.FastBufferedOutputStream;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * Converts a <a href="http://en.wikipedia.org/wiki/FASTA_format">FASTA</a>
 * or <a href="http://en.wikipedia.org/wiki/FASTQ_format">FASTQ</a> file to the
 * protocol buffer file format described by Reads.proto.
 *
 * @author Fabien Campagne
 *         Date: Apr 28, 2009
 *         Time: 6:03:56 PM
 */
public class FastaToCompactMode extends AbstractGobyMode {
    /**
     * Used to log informational and debug messages.
     */
    private static final Log LOG = LogFactory.getLog(FastaToCompactMode.class);

    /**
     * The mode name.
     */
    private static final String MODE_NAME = "fasta-to-compact";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Converts a FASTA/FASTQ file to the "
            + "Goby \"compact-reads\" file format.";

    /**
     * The files to convert to compact reads.
     */
    private String[] inputFilenames;

    /**
     * The output file or basename of the output files if there is more than one input file.
     */
    private String outputFile;

    /**
     * Include descriptions in the compact output.
     */
    private boolean includeDescriptions;

    /**
     * Include identifiers in the compact output.
     */
    private boolean includeIdentifiers;

    /**
     * The number of sequences that will be written in each compressed chunk. Th default is
     * suitable for very many short sequences but should be reduced to a few sequences per
     * chunk if each sequence is very large.
     */
    private int sequencePerChunk = 10000;

    /**
     * Exclude sequence data from the compact output.
     */
    private boolean excludeSequences;

    /**
     * Exclude quality scores from the compact output.
     */
    private boolean excludeQuality;
    private boolean verboseQualityScores;

    private boolean parallel = true;

    private QualityEncoding qualityEncoding;

    /**
     * Returns the mode name defined by subclasses.
     *
     * @return The name of the mode
     */
    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    /**
     * Returns the mode description defined by subclasses.
     *
     * @return A description of the mode
     */
    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }

    /**
     * Configure.
     *
     * @param args command line arguments
     * @return this object for chaining
     * @throws IOException   error parsing
     * @throws JSAPException error parsing
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args) throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);
        parallel = jsapResult.getBoolean("parallel", false);
        inputFilenames = jsapResult.getStringArray("input");
        includeDescriptions = jsapResult.getBoolean("include-descriptions");
        includeIdentifiers = jsapResult.getBoolean("include-identifiers");
        excludeSequences = jsapResult.getBoolean("exclude-sequences");
        excludeQuality = jsapResult.getBoolean("exclude-quality");
        verboseQualityScores = jsapResult.getBoolean("verbose-quality-scores");
        qualityEncoding = QualityEncoding.valueOf(jsapResult.getString("quality-encoding").toUpperCase());
        if (inputFilenames.length == 1) {
            outputFile = jsapResult.getString("output");
        }
        sequencePerChunk = jsapResult.getInt("sequence-per-chunk");
        return this;
    }

    /**
     * Perform the conversion fasta -> compact-reads on one or more files.
     *
     * @throws IOException if the input/output files cannot be read/written
     */
    @Override
    public void execute() throws IOException {
        //    final int numToProcess = inputFilenames.length;
        //  int numProcessed = 0;
        //   for (final String inputFilename : inputFilenames) {
        //       numProcessed = processOneFile(numToProcess, numProcessed, inputFilename);
        //  }

        try {
            final DoInParallel loop = new DoInParallel() {
                @Override
                public void action(final DoInParallel forDataAccess, final String inputBasename, final int loopIndex) {

                    try {
                        debugStart(inputBasename);
                        processOneFile(loopIndex, inputFilenames.length, inputBasename);
                        debugEnd(inputBasename);
                    } catch (IOException e) {
                        LOG.error(e);
                    }
                }
            };
            System.out.println("parallel: " + parallel);
            loop.execute(parallel, inputFilenames);
        } catch (Exception e) {
            LOG.error(e);
        }
    }


    private void processOneFile(final int loopIndex, final int length, final String inputFilename) throws IOException {
        final String outputFilename;
        if (loopIndex == 0 && StringUtils.isNotBlank(outputFile)) {
            outputFilename = outputFile;
        } else {
            outputFilename = stripFastxExtensions(inputFilename) + ".compact-reads";
        }
        System.out.println("Creating file " + outputFilename);
        final File output = new File(outputFilename);
        final File readsFile = new File(inputFilename);
        if (!output.exists() || FileUtils.isFileNewer(readsFile, output) || output.length() == 0) {
            convert(loopIndex, length, inputFilename, outputFilename);
        }

    }

    private void convert(final int loopIndex, final int length, final String inputFilename, final String outputFilename) throws IOException {
        System.out.printf("Converting [%d/%d] %s to %s%n",
                loopIndex, length, inputFilename, outputFilename);

        // Create directory for output file if it doesn't already exist
        final String outputPath = FilenameUtils.getFullPath(outputFilename);
        if (StringUtils.isNotBlank(outputPath)) {
            FileUtils.forceMkdir(new File(outputPath));
        }
        final ReadsWriter writer = new ReadsWriter(new FastBufferedOutputStream(new FileOutputStream(outputFilename)));
        try {

            writer.setNumEntriesPerChunk(sequencePerChunk);
            for (final FastXEntry entry : new FastXReader(inputFilename)) {
                if (includeDescriptions) {
                    writer.setDescription(entry.getEntryHeader());
                }
                if (includeIdentifiers) {
                    final MutableString description = entry.getEntryHeader();
                    final String identifier = description.toString().split("[\\s]")[0];
                    writer.setIdentifier(identifier);
                }
                if (!excludeSequences) {
                    writer.setSequence(entry.getSequence());
                } else {
                    writer.setSequence("");
                }
                if (!excludeQuality) {
                    writer.setQualityScores(convertQualityScores(entry.getQuality()));
                }
                writer.appendEntry();
            }
        } finally {

            writer.close();
            writer.printStats(System.out);
        }

    }

    private byte[] convertQualityScores(final MutableString quality) {
        // Only sanger and illumina encoding are supported at this time
        if (qualityEncoding == QualityEncoding.SOLEXA) {
            throw new UnsupportedOperationException("SOLEXA encoding is not supported "
                    + "at this time for lack of clear documentation.");
        } else if (qualityEncoding != QualityEncoding.SANGER
                && qualityEncoding != QualityEncoding.ILLUMINA) {
            throw new UnsupportedOperationException("Unknown encoding: " + qualityEncoding);
        }

        final int size = quality.length();
        final byte[] qualityScoreBuffer = new byte[size];

        if (verboseQualityScores) {
            System.out.println(quality);
        }

        for (int position = 0; position < size; position++) {
            qualityScoreBuffer[position] =
                    qualityEncoding.asciiEncodingToQualityScore(quality.charAt(position));
            checkRange(qualityScoreBuffer[position]);
            if (verboseQualityScores) {
                System.out.print(qualityScoreBuffer[position]);
                System.out.print(" ");
            }
        }
        if (verboseQualityScores) {
            System.out.println();
        }

        return qualityScoreBuffer;
    }

    private void checkRange(final byte qualityDecoded) {
        if (!(qualityDecoded >= 0 && qualityDecoded <= 40)) {
            System.err.println("Phred quality scores must be within 0 and 40. The value decoded "
                    + "was " + qualityDecoded + " You may have selected an incorrect encoding.");
            System.exit(10);
        }
    }

    /**
     * Get the filename including path WITHOUT fastx extensions (including .gz if it is there).
     *
     * @param name the full path to the file in question
     * @return the full path to file without the fastx/gz extensions or the same name if
     * those extensions weren't found.
     * @see edu.cornell.med.icb.goby.util.FileExtensionHelper#FASTX_FILE_EXTS
     */
    private static String stripFastxExtensions(final String name) {
        final String filename = FilenameUtils.getName(name);
        for (final String ext : FileExtensionHelper.FASTX_FILE_EXTS) {
            if (filename.endsWith(ext)) {
                return FilenameUtils.getFullPath(name)
                        + filename.substring(0, filename.lastIndexOf(ext));
            }
        }
        return name;
    }

    public static void main(final String[] args) throws IOException, JSAPException {
        new FastaToCompactMode().configure(args).execute();
    }
}
