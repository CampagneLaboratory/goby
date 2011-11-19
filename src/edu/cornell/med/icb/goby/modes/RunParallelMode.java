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
import edu.cornell.med.icb.goby.GobyVersion;
import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.alignments.ReferenceLocation;
import edu.cornell.med.icb.goby.alignments.UpgradeTo1_9_6;
import edu.cornell.med.icb.goby.util.DoInParallel;
import edu.cornell.med.icb.goby.util.LoggingOutputStream;
import edu.rit.pj.ParallelRegion;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.exec.CommandLine;
import org.apache.commons.exec.DefaultExecutor;
import org.apache.commons.exec.PumpStreamHandler;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.apache.log4j.Level;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

/**
 * Converts a compact alignment to plain text.
 *
 * @author Fabien Campagne
 */
public class RunParallelMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "run-parallel";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Run some command in parallel for parts of a compact-reads file.";
    private int numParts;
    private String processPartCommand;
    private String input;
    private String output;


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

        input = jsapResult.getString("input");
        output = jsapResult.getString("output");
        processPartCommand = jsapResult.getString("process-part-command");
        numParts = jsapResult.getInt("num-parts");

        return this;


    }

    class Slice {
        long startOffset;
        long endOffset;
    }

    public void execute() throws IOException {
        final Slice slices[] = new Slice[numParts];
        File file = new File(input);
        if (!(file.isFile() && file.exists() && file.canRead())) {
            System.err.println("Input file cannot be read: " + input);
            System.exit(1);
        }
        int i = 0;
        for (final Slice slice : slices) {

            slices[i++] = new Slice();
        }
        final long fileLength = file.length();
        final long sliceLength = fileLength / numParts;
        long currentOffset = 0;

        for (final Slice slice : slices) {

            slice.startOffset = currentOffset;
            slice.endOffset = currentOffset + sliceLength;
            currentOffset = slice.endOffset;
        }

        final ObjectOpenHashSet<String> allOutputs = new ObjectOpenHashSet<String>();
        final ObjectOpenHashSet<String> allFastq = new ObjectOpenHashSet<String>();

        final DoInParallel loop = new DoInParallel(numParts) {
            @Override
            public void action(final DoInParallel forDataAccess, final String inputBasename, final int loopIndex) {
                try {

                    CompactToFastaMode ctfm = new CompactToFastaMode();
                    ctfm.setInputFilename(input);
                    ctfm.setOutputFormat(CompactToFastaMode.OutputFormat.FASTQ);
                    ctfm.setStartPosition(slices[loopIndex].startOffset);
                    ctfm.setEndPosition(slices[loopIndex].endOffset);

                    String s = FilenameUtils.getBaseName(FilenameUtils.removeExtension(input));
                    String fastqFilename = s + "-" + loopIndex + ".fq";
                    allFastq.add(fastqFilename);
                    File tmp1 = File.createTempFile(s, "-tmp");
                    File output = File.createTempFile(s, "-out");
                    ctfm.setOutputFilename(fastqFilename);
                    LOG.info(String.format("Extracting FASTQ for slice [%d-%d]%n",slices[loopIndex].startOffset,slices[loopIndex].endOffset));
                    ctfm.execute();
                    if (loopIndex>0) {
                       // wait a bit to give the first thread the time to load the database and establish shared memory pool
                       sleep(60);
                    }
                    String transformedCommand = processPartCommand.replaceAll("%read.fastq%", fastqFilename);
                    transformedCommand = transformedCommand.replaceAll("%tmp1%", tmp1.getName());
                    String outputFilename = output.getName();
                    transformedCommand = transformedCommand.replaceAll("%output%", outputFilename);
                    final DefaultExecutor executor = new DefaultExecutor();
                    OutputStream logStream = null;
                    try {
                        logStream = new LoggingOutputStream(getClass(), Level.INFO, "");
                        executor.setStreamHandler(new PumpStreamHandler(logStream));
                        LOG.info("About to execute: " + transformedCommand);
                        final int exitValue = executor.execute(CommandLine.parse(transformedCommand));
                        LOG.info("Exit value = " + exitValue);
                    } finally {
                        IOUtils.closeQuietly(logStream);
                        // remove the fastq file
                        new File(fastqFilename).delete();
                    }

                    if (new File(outputFilename + ".header").exists()) {
                        // found output alignment:
                        System.out.println("found output file: " + outputFilename);
                        allOutputs.add(outputFilename);
                    } else {
                        System.out.println("Warning: did not find output alignment: " + outputFilename);
                    }

                } catch (IOException e) {
                    LOG.error("Error processing index " + loopIndex + ", " + inputBasename, e);
                }
            }
        };
        String[] parts = new String[numParts];


        try {
            loop.execute(true, parts);
        } catch (Exception e) {
            System.err.println("An error occurred executing a parallel command: ");
            e.printStackTrace();
        }


        final ConcatenateAlignmentMode concat = new ConcatenateAlignmentMode();
        concat.setInputFileNames(allOutputs.toArray(new String[allOutputs.size()]));
        concat.setOutputFilename(output);
        concat.setAdjustQueryIndices(false);
        concat.setAdjustSampleIndices(false);
        concat.execute();

    }

    /**
     * Sleep for the specified number of seconds.
     * @param seconds
     */
    private void sleep(int seconds) {
        try {
            Thread.sleep(1000 *seconds);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }


    /**
     * Main method.
     *
     * @param args command line args.
     * @throws com.martiansoftware.jsap.JSAPException
     *                             error parsing
     * @throws java.io.IOException error parsing or executing.
     */

    public static void main(final String[] args) throws JSAPException, IOException {
        new RunParallelMode().configure(args).execute();
    }

}