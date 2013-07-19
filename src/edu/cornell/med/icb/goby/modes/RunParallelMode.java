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
import edu.cornell.med.icb.goby.util.DoInParallel;
import edu.cornell.med.icb.goby.util.IsDone;
import edu.cornell.med.icb.goby.util.LoggingOutputStream;
import edu.cornell.med.icb.goby.util.StreamSignal;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import org.apache.commons.exec.CommandLine;
import org.apache.commons.exec.DefaultExecutor;
import org.apache.commons.exec.PumpStreamHandler;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.apache.log4j.Level;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.HashMap;
import java.util.Map;

/**
 * Run some command in parallel for parts of a compact-reads file.
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
    private String[] command;


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
        command = jsapResult.getStringArray("process-part-command");
        paired = jsapResult.getBoolean("paired");

        StringBuffer sb = new StringBuffer();
        for (String arg : command) {
            sb.append(arg);
            sb.append(" ");
            if (!paired && args.equals("%pair.fastq%")) {
                System.err.println("%pair.fastq% requires the --paired argument.");
                System.exit(1);
            }
        }
        processPartCommand = sb.toString();
        numParts = jsapResult.getInt("num-parts");

        return this;


    }

    class Slice {
        long startOffset;
        long endOffset;
    }

    private boolean paired = false;

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
            IsDone done = new IsDone();

            @Override
            public void action(final DoInParallel forDataAccess, final String inputBasename, final int loopIndex) {
                try {

                    CompactToFastaMode ctfm = new CompactToFastaMode();
                    ctfm.setInputFilename(input);
                    ctfm.setOutputFormat(CompactToFastaMode.OutputFormat.FASTQ);
                    ctfm.setStartPosition(slices[loopIndex].startOffset);
                    ctfm.setEndPosition(slices[loopIndex].endOffset);
                    String s = FilenameUtils.getBaseName(FilenameUtils.removeExtension(input)) + "-" + Integer.toString(loopIndex);
                    String fastqFilename = s + "-input.fq";
                    String fastqPairedFilename = s + "-pair-input.fq";
                    allFastq.add(fastqFilename);
                    allFastq.add(fastqPairedFilename);
                    File tmp1 = new File(s + "-tmp");
                    tmp1.deleteOnExit();
                    File output = new File(s + "-out");
                    output.deleteOnExit();
                    ctfm.setOutputFilename(fastqFilename);
                    if (paired) {
                        ctfm.setOutputPairFilename(fastqPairedFilename);
                    }
                    LOG.info(String.format("Extracting FASTQ for slice [%d-%d] loopIndex=%d %n",
                            slices[loopIndex].startOffset,
                            slices[loopIndex].endOffset, loopIndex));
                    ctfm.execute();
                    if (loopIndex > 0) {
                        while (!done.isDone()) {
                            // wait a bit to give the first thread the time to load the database and establish shared memory pool
                            //   System.out.println("sleep 5 thread "+loopIndex);
                            sleep(5);
                        }
                        System.out.println("Thread " + loopIndex + " can now start.");
                    }
                    final Map<String, String> replacements = new HashMap<String, String>();

                    final String outputFilename = output.getName();

                    replacements.put("%read.fastq%", fastqFilename);
                    if (paired) {
                        replacements.put("%pair.fastq%", fastqPairedFilename);
                    }
                    replacements.put("%tmp1%", tmp1.getName());
                    replacements.put("%output%", outputFilename);
                    final String transformedCommand = transform(processPartCommand, replacements);
                    final DefaultExecutor executor = new DefaultExecutor();
                    OutputStream logStream = null;
                    try {
                        logStream = new LoggingOutputStream(getClass(), Level.INFO, "");
                        executor.setStreamHandler(new PumpStreamHandler(new StreamSignal(done, "scanning", logStream)));

                        final CommandLine parse = CommandLine.parse(transformedCommand, replacements);
                        LOG.info("About to execute: " + parse);
                        final int exitValue = executor.execute(parse);
                        LOG.info("Exit value = " + exitValue);
                        if (new File(outputFilename + ".header").exists()) {
                            // found output alignment:
                            System.out.println("found output file: " + outputFilename);
                            allOutputs.add(outputFilename + ".header");
                        } else {
                            System.out.println("Warning: did not find output alignment: " + outputFilename);
                        }
                    } finally {
                        IOUtils.closeQuietly(logStream);
                        // remove the fastq file
                        new File(fastqFilename).delete();
                    }


                } catch (IOException e) {
                    LOG.error("Error processing index " + loopIndex + ", " + inputBasename, e);
                }
            }
        };
        String[] parts = new String[numParts];

        for (int j = 0; j < numParts; j++) {
            parts[j] = Integer.toString(j);
        }
        try {
            loop.execute(true, parts);
        } catch (Exception e) {
            System.err.println("An error occurred executing a parallel command: ");
            e.printStackTrace();
        }

        System.out.printf("Preparing to concatenate %d outputs..%n", allOutputs.size());
        final ConcatenateAlignmentMode concat = new ConcatenateAlignmentMode();
        concat.setInputFileNames(allOutputs.toArray(new String[allOutputs.size()]));
        concat.setOutputFilename(output);
        concat.setAdjustQueryIndices(false);
        concat.setAdjustSampleIndices(false);
        concat.execute();

    }

    private String transform(String processPartCommand, Map<String, String> replacements) {
        for (String key : replacements.keySet()) {
            processPartCommand = processPartCommand.replaceAll(key, replacements.get(key));
        }
        return processPartCommand;
    }

    /**
     * Sleep for the specified number of seconds.
     *
     * @param seconds
     */
    private void sleep(int seconds) {
        try {
            Thread.sleep(1000 * seconds);
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