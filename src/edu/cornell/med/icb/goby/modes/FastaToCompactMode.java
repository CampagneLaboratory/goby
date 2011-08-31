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
import edu.cornell.med.icb.goby.reads.QualityEncoding;
import edu.cornell.med.icb.goby.reads.ReadCodec;
import edu.cornell.med.icb.goby.reads.ReadsWriter;
import edu.cornell.med.icb.goby.stats.NormalizationMethod;
import edu.cornell.med.icb.goby.util.DoInParallel;
import edu.cornell.med.icb.goby.util.FileExtensionHelper;
import it.unimi.dsi.fastutil.io.FastBufferedOutputStream;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.omg.IOP.Codec;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.FileReader;
import java.util.Properties;
import java.util.List;
import java.util.LinkedList;
import java.util.ServiceLoader;

/**
 * Converts a <a href="http://en.wikipedia.org/wiki/FASTA_format">FASTA</a>
 * or <a href="http://en.wikipedia.org/wiki/FASTQ_format">FASTQ</a> file to the compact reads
 * format. Compact reads are in the chunked protocol buffer file format described by Reads.proto.
 * Since Goby 1.7, this mode can load paired-end runs into single compact files.
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
    private static final String MODE_DESCRIPTION = "Converts FASTA/FASTQ files to the "
            + "Goby \"compact-reads\" file format. Since Goby 1.7, this mode can load paired-end "
            + "runs into single compact files.";

    /**
     * The files to convert to compact reads.
     */
    private String[] inputFilenames;
    // The list version is used by the API with addInputFilename()
    private List<String> inputFilenamesList = null;

    /**
     * The output file or basename of the output files if there is more than one input file.
     */
    private String outputFilename;

    /**
     * Include descriptions in the compact output.
     */
    private boolean includeDescriptions = false;

    /**
     * Include identifiers in the compact output.
     */
    private boolean includeIdentifiers = false;

    /**
     * The number of sequences that will be written in each compressed chunk. Th default is
     * suitable for very many short sequences but should be reduced to a few sequences per
     * chunk if each sequence is very large.
     */
    private int sequencePerChunk = 10000;

    /**
     * Exclude sequence data from the compact output.
     */
    private boolean excludeSequences = false;

    /**
     * Exclude quality scores from the compact output.
     */
    private boolean excludeQuality = false;
    private boolean verboseQualityScores = false;

    private boolean parallel = false;

    private QualityEncoding qualityEncoding = QualityEncoding.ILLUMINA;
    private boolean processPairs = false;
    private String pairIndicator1 = null;
    private String pairIndicator2 = null;
    private File keyValuePairsFilename;
    private Properties keyValueProps;

    private boolean apiMode = true;
    private int numThreads;

    private ReadCodec codec;

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
     * Get the list of files (reads/alignments) to process.
     *
     * @return The list of files.
     */
    public String[] getInputFilenames() {
        return inputFilenames;
    }

    /**
     * Add the specified filename to the list of files to process.
     *
     * @param inputFilename The file to process
     */
    public synchronized void addInputFilename(final String inputFilename) {
        if (inputFilenamesList == null) {
            inputFilenamesList = new LinkedList<String>();
        }
        inputFilenamesList.add(inputFilename);
        inputFilenames = inputFilenamesList.toArray(new String[inputFilenamesList.size()]);
    }

    public synchronized void clearInputFilenames() {
        if (inputFilenamesList == null) {
            inputFilenamesList = new LinkedList<String>();
        } else {
            inputFilenamesList.clear();
        }
        inputFilenames = new String[0];
    }

    /**
     * Get the output filename. Only used when inputFilenames.length == 1.
     *
     * @return the output filename
     */
    public String getOutputFilename() {
        return outputFilename;
    }

    /**
     * Set the output filename. Only used when inputFilenames.length == 1.
     *
     * @param outputFilename the output filename
     */
    public void setOutputFilename(final String outputFilename) {
        this.outputFilename = outputFilename;
    }

    /**
     * Get the quality encoding scale used for the input fastq file.
     *
     * @return the quality encoding scale used for the input fastq file
     */
    public QualityEncoding getQualityEncoding() {
        return qualityEncoding;
    }

    /**
     * Set the quality encoding scale to be used for the input fastq file.
     * Acceptable values are "Illumina", "Sanger", and "Solexa".
     *
     * @param qualityEncoding the quality encoding scale to be used for the input fastq file
     */
    public void setQualityEncoding(final QualityEncoding qualityEncoding) {
        this.qualityEncoding = qualityEncoding;
    }

    /**
     * Set the quality encoding scale to be used for the input fastq file.
     * Acceptable values are "Illumina", "Sanger", and "Solexa".
     *
     * @param qualityEncoding the quality encoding scale to be used for the input fastq file
     */
    public void setQualityEncoding(final String qualityEncoding) {
        this.qualityEncoding = QualityEncoding.valueOf(qualityEncoding.toUpperCase());
    }

    private static final ServiceLoader<ReadCodec> codecLoader = ServiceLoader.load(ReadCodec.class);

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
        this.apiMode = false;
        final JSAPResult jsapResult = parseJsapArguments(args);
        parallel = jsapResult.getBoolean("parallel", false);
        inputFilenames = jsapResult.getStringArray("input");
        includeDescriptions = jsapResult.getBoolean("include-descriptions");
        includeIdentifiers = jsapResult.getBoolean("include-identifiers");
        excludeSequences = jsapResult.getBoolean("exclude-sequences");
        excludeQuality = jsapResult.getBoolean("exclude-quality");
        verboseQualityScores = jsapResult.getBoolean("verbose-quality-scores");
        qualityEncoding = QualityEncoding.valueOf(jsapResult.getString("quality-encoding").toUpperCase());
        numThreads = jsapResult.getInt("num-threads");
        if (inputFilenames.length == 1) {
            outputFilename = jsapResult.getString("output");
        }
        sequencePerChunk = jsapResult.getInt("sequence-per-chunk");
        processPairs = jsapResult.getBoolean("paired-end");
        final String tokens = jsapResult.getString("pair-indicator");

        String codecName = jsapResult.getString("codec");

        if (codecName != null) {
            codecLoader.reload();
            for (final ReadCodec c : codecLoader) {

                if (c.name().equals(codecName)) {
                    LOG.info("Will use read codec " + c.name());
                    codec = c;
                    break;
                }
            }
        }
        if (processPairs && tokens != null) {
            final String[] tmp = tokens.split("[,]");
            if (tmp.length != 2) {
                System.err.println("Pair indicator argument must have exactly two tokens, separated by " +
                        "comma. Offending syntax: " + tokens);
                System.exit(1);
            }
            pairIndicator1 = tmp[0];
            pairIndicator2 = tmp[1];
        }
        keyValueProps = new Properties();
        keyValuePairsFilename = jsapResult.getFile("key-value-pairs");

        if (keyValuePairsFilename != null) {
            if (!keyValuePairsFilename.exists()) {
                System.err.println("Key value pair file cannot be found. Please check filename and try again.");
                System.exit(1);
            } else {

                keyValueProps.load(new FileReader(keyValuePairsFilename));

            }
        }
        String keys[] = jsapResult.getStringArray("key");
        String values[] = jsapResult.getStringArray("value");
        if (keys.length != values.length) {
            System.out.println("Key and value arguments must be paired.");
            System.exit(1);
        }
        int length = keys.length;
        for (int i = 0; i < length; i++) {
            String key = keys[i];
            String value = values[i];
            // override same keys with the command line values:
            keyValueProps.put(key, value);
            System.out.printf("defining or overriding key=%s value=%s %n", key, value);
        }
        return this;
    }

    /**
     * Perform the conversion fasta -> compact-reads on one or more files.
     *
     * @throws IOException if the input/output files cannot be read/written
     */
    @Override
    public void execute() throws IOException {


        try {
            if (apiMode) {
                // Force parallel to false if in apiMode.
                parallel = false;
            }
            final DoInParallel loop = new DoInParallel(numThreads) {
                @Override
                public void action(final DoInParallel forDataAccess, final String inputBasename, final int loopIndex) {
                    try {
                        debugStart(inputBasename);
                        processOneFile(loopIndex, inputFilenames.length, inputBasename, keyValueProps);
                        debugEnd(inputBasename);
                    } catch (IOException e) {
                        LOG.error("Error processing index " + loopIndex + ", " + inputBasename, e);
                    }
                }
            };
            System.out.println("parallel: " + parallel);
            loop.execute(parallel, inputFilenames);
        } catch (IllegalArgumentException e) {
            if (apiMode) {
                throw e;
            } else {
                LOG.error("Error processing", e);
            }
        } catch (Exception e) {
            LOG.error("Error processing", e);
        }
    }


    private void processOneFile(final int loopIndex, final int length, final String inputFilename,
                                Properties keyValueProps) throws IOException {
        String outputFilename;
        if (loopIndex == 0 && StringUtils.isNotBlank(this.outputFilename)) {
            outputFilename = this.outputFilename;
        } else {
            outputFilename = stripFastxExtensions(inputFilename) + ".compact-reads";

        }
        if (processPairs) {
            // remove _1 from the destination compact filename.
            outputFilename = outputFilename.replace(pairIndicator1, "");
        }

        final File output = new File(outputFilename);
        final File readsFile = new File(inputFilename);
        if (!output.exists() || FileUtils.isFileNewer(readsFile, output) || output.length() == 0) {
            System.out.println("Creating file " + outputFilename);
            convert(loopIndex, length, inputFilename, outputFilename, keyValueProps);
        } else {
            System.out.printf("Skipping file that already exists %s%n", outputFilename);
        }
    }

    private void convert(final int loopIndex, final int length, final String inputFilename, final String outputFilename, Properties keyValueProps) throws IOException {
        System.out.printf("Converting [%d/%d] %s to %s%n",
                loopIndex + 1, length, inputFilename, outputFilename);

        // Create directory for output file if it doesn't already exist
        final String outputPath = FilenameUtils.getFullPath(outputFilename);
        if (StringUtils.isNotBlank(outputPath)) {
            FileUtils.forceMkdir(new File(outputPath));
        }
        final ReadsWriter writer = new ReadsWriter(new FastBufferedOutputStream(new FileOutputStream(outputFilename)));
        if (codec != null) {
            writer.setCodec(codec);
        }
        try {

            writer.setNumEntriesPerChunk(sequencePerChunk);
            FastXReader pairReader = null;
            if (processPairs) {
                final String pairInputFilename = inputFilename.replace(pairIndicator1, pairIndicator2);
                LOG.info(String.format("Located paired-end input files (%s,%s)", inputFilename, pairInputFilename));
                pairReader = new FastXReader(pairInputFilename);

            }
            FastXEntry pairEntry = null;
            int entryIndex = 0;
            writer.setMetaData(keyValueProps);

            for (final FastXEntry entry : new FastXReader(inputFilename)) {
                if (pairReader != null) {
                    pairEntry = pairReader.next();
                    if (pairEntry == null) {
                        System.err.println("Cannot find matching sequence for " + entry.getEntryHeader());
                    }
                }
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
                    if (pairEntry != null) {
                        writer.setPairSequence(pairEntry.getSequence());
                    }
                } else {
                    writer.setSequence("");
                }
                if (!excludeQuality) {
                    writer.setQualityScores(convertQualityScores(qualityEncoding, entry.getQuality(),
                            verboseQualityScores, apiMode));
                    if (pairEntry != null) {
                        writer.setQualityScoresPair(convertQualityScores(qualityEncoding, pairEntry.getQuality(),
                                verboseQualityScores, apiMode));
                    }
                }

                writer.appendEntry();
                entryIndex++;
            }
        } finally {
            writer.close();
            writer.printStats(System.out);
        }
    }

    static public byte[] convertQualityScores(QualityEncoding qualityEncoding, final CharSequence quality, boolean verboseQualityScores) {
        return convertQualityScores(qualityEncoding, quality, verboseQualityScores, false);
    }

    static public byte[] convertQualityScores(QualityEncoding qualityEncoding, final CharSequence quality, boolean verboseQualityScores, boolean apiMode) {
        // Only Solexa, Sanger and Illumina encoding are supported at this time

        final int size = quality.length();
        final byte[] qualityScoreBuffer = new byte[size];

        if (verboseQualityScores) {
            System.out.println(quality);
        }

        for (int position = 0; position < size; position++) {
            qualityScoreBuffer[position] =
                    qualityEncoding.asciiEncodingToPhredQualityScore(quality.charAt(position));

            if (!qualityEncoding.isWithinValidRange(qualityScoreBuffer[position])) {
                final String message = "Phred quality scores must be within specific ranges for specfic encodings. " +
                        "The value decoded was " + qualityScoreBuffer[position] +
                        " and outside of the valid range for " + qualityEncoding +
                        " You may have selected an incorrect encoding.";
                if (apiMode) {
                    throw new IllegalArgumentException(message);
                } else {
                    System.err.println(message);
                    System.exit(10);
                }
            }
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
