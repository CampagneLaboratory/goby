/*
 * Copyright (C) 2009 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 * WEILL MEDICAL COLLEGE OF CORNELL UNIVERSITY MAKES NO REPRESENTATIONS
 * ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY PURPOSE. IT IS PROVIDED
 * "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY. THE WEILL MEDICAL COLLEGE
 * OF CORNELL UNIVERSITY SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY
 * THE USERS OF THIS SOFTWARE.
 */

package edu.cornell.med.icb.goby.modes;

import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.fastxReaders.FastXEntry;
import edu.cornell.med.icb.goby.fastxReaders.FastXReader;
import edu.cornell.med.icb.goby.reads.ReadsWriter;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.StringUtils;

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


    private String[] inputFilenames;
    private boolean pushDescription;
    private boolean pushIdentifier;
    private String outputFile;

    /**
     * The mode name.
     */
    public static final String MODE_NAME = "fasta-to-compact";
    public static final String MODE_DESCRIPTION = "Converts a FASTA file to the protocol buffer file format described by Reads.proto.";

    private int sequencePerChunk = 10000;
    private boolean excludeSequences;

    public String getModeName() {
        return MODE_NAME;
    }

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
        inputFilenames = jsapResult.getStringArray("input");
        pushDescription = jsapResult.getBoolean("include-descriptions");
        pushIdentifier = jsapResult.getBoolean("include-identifiers");
        excludeSequences = jsapResult.getBoolean("exclude-sequences");
        outputFile = jsapResult.getString("output");
        sequencePerChunk = jsapResult.getInt("sequence-per-chunk");
        return this;
    }

    /**
     * Perform the conversion fasta -> compact-reads on one or more files.
     *
     * @throws IOException
     */
    @Override
    public void execute() throws IOException {
        final int numToProcess = inputFilenames.length;
        int numProcessed = 0;
        for (final String inputFilename : inputFilenames) {

            final String outputFilename;
            if (numToProcess == 1 && StringUtils.isNotBlank(outputFile)) {
                outputFilename = outputFile;
            } else {
                outputFilename = stripFastxExtensions(inputFilename) + ".compact-reads";
            }
            System.out.printf("Converting [%d/%d] %s to %s%n",
                    ++numProcessed, numToProcess, inputFilename, outputFilename);

            final ReadsWriter writer = new ReadsWriter(new FileOutputStream(outputFilename));
            writer.setNumEntriesPerChunk(sequencePerChunk);
            for (final FastXEntry entry : new FastXReader(inputFilename)) {
                if (pushDescription) {
                    writer.setDescription(entry.getEntryHeader());
                }
                if (pushIdentifier) {
                    final MutableString description = entry.getEntryHeader();
                    final String identifier = description.toString().split("[\\s]")[0];
                    writer.setIdentifier(identifier);
                }
                if (!excludeSequences){
                    writer.setSequence(entry.getSequence());
                } else {
                    writer.setSequence("");
                }
                writer.appendEntry();
            }

            writer.close();
            writer.printStats(System.out);
        }
    }

    /**
     * Get the filename including path WITHOUT fastx extensions (including .gz if it is there).
     *
     * @param name the full path to the file in question
     * @return the filename without the fastx/gz extensions or the same name of those extensions
     *         weren't found.
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
