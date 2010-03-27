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
import edu.cornell.med.icb.goby.reads.ReadsWriter;
import edu.cornell.med.icb.goby.util.FileExtensionHelper;
import edu.cornell.med.icb.goby.util.DoInParallel;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.StringUtils;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * Converts a FASTA file to the protocol buffer file format described by Reads.proto.
 *
 * @author Fabien Campagne
 *         Date: Apr 28, 2009
 *         Time: 6:03:56 PM
 */
public class FastaToCompactMode extends AbstractGobyMode {
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
    private byte[] qualityScoreBuffer;
    private boolean parallel = true;

    // see http://en.wikipedia.org/wiki/FASTQ_format for a description of these various encodings
    /*

    Sanger format can encode a Phred quality score from 0 to 93 using ASCII 33 to 126 (although in raw read data the Phred quality score rarely exceeds 60, higher scores are possible in assemblies or read maps).
    Illumina 1.3+ format can encode a Phred quality score from 0 to 62 using ASCII 64 to 126 (although in raw read data Phred scores from 0 to 40 only are expected).
    Solexa/Illumina 1.0 format can encode a Solexa/Illumina quality score from -5 to 62 using ASCII 59 to 126 (although in raw read data Solexa scores from -5 to 40 only are expected)
  <PRE>
  SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
  ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
  ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
  |                         |    |        |                              |                     |
 33                        59   64       73                            104                   126
       </PRE>
 S - Sanger       Phred+33,  41 values  (0, 40)
 I - Illumina 1.3 Phred+64,  41 values  (0, 40)
 X - Solexa       Solexa+64, 68 values (-5, 62)
     */
    public enum QualityEncoding {
        SANGER, ILLUMINA, SOLEXA
    }

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
        outputFile = jsapResult.getString("output");
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
            DoInParallel loop = new DoInParallel() {
                @Override
                public void action(DoInParallel forDataAccess, String inputBasename, int loopIndex) {

                    try {
                        debugStart(inputBasename);
                        processOneFile(loopIndex, inputFilenames.length,  inputBasename);
                        debugEnd(inputBasename);
                    } catch (IOException e) {
                        LOG.error(e);
                    }
                }
            };
            System.out.println("parallel: "+parallel);
            loop.execute(parallel, inputFilenames);
        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }


    private void processOneFile(int loopIndex, int length, String inputFilename) throws IOException {
        final String outputFilename;
        if (loopIndex == 1 && StringUtils.isNotBlank(outputFile)) {
            outputFilename = outputFile;
        } else {
            outputFilename = stripFastxExtensions(inputFilename) + ".compact-reads";
        }
        final File output = new File(outputFilename);
        final File readsFile = new File(inputFilename);
        if (!output.exists() || FileUtils.isFileNewer(readsFile, output) || output.length()==0) {
             convert(loopIndex, length, inputFilename, outputFilename);
        }

    }

    private void convert(int loopIndex, int length, String inputFilename, String outputFilename) throws IOException {
        System.out.printf("Converting [%d/%d] %s to %s%n",
                loopIndex, length,inputFilename, outputFilename);

        // Create directory for output file if it doesn't already exist
        final String outputPath = FilenameUtils.getFullPath(outputFilename);
        if (StringUtils.isNotBlank(outputPath)) {
            FileUtils.forceMkdir(new File(outputPath));
        }
        final ReadsWriter writer = new ReadsWriter(new FileOutputStream(outputFilename));
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
        final int size = quality.length();
        if (qualityScoreBuffer == null || size != qualityScoreBuffer.length) {
            qualityScoreBuffer = new byte[size];
        }
        if (verboseQualityScores) {
            System.out.println(quality);
        }

        switch (qualityEncoding) {
            case SANGER:
                for (int position = 0; position < size; position++) {
                    qualityScoreBuffer[position] = (byte) (quality.charAt(position) - 33);
                    checkRange(qualityScoreBuffer[position]);
                    if (verboseQualityScores) {
                        System.out.print(qualityScoreBuffer[position]);
                        System.out.print(" ");
                    }
                }
                if (verboseQualityScores) {
                    System.out.println("");
                }
                break;
            case ILLUMINA:
                for (int position = 0; position < size; position++) {
                    qualityScoreBuffer[position] = (byte) (quality.charAt(position) - 64);
                    checkRange(qualityScoreBuffer[position]);

                    if (verboseQualityScores) {
                        System.out.print(qualityScoreBuffer[position]);
                        System.out.print(" ");
                    }
                }
                if (verboseQualityScores) {
                    System.out.println("");
                }
                break;
            case SOLEXA:
                throw new UnsupportedOperationException("SOLEXA encoding is not supported at this time for lack of clear documentation.");
                /* for (int position = 0; position < size; position++) {
                   final int code = quality.charAt(position) - 64;
                   double phred=Math.log(1+Math.pow(10,code/10))/Math.log(10);
                   qualityScoreBuffer[position] = (byte) phred;
               break;
               } */

            default:
                throw new UnsupportedOperationException("Unknown encoding: " + qualityEncoding);
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
     *         those extensions weren't found.
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
