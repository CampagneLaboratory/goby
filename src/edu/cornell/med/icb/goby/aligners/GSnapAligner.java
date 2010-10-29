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

package edu.cornell.med.icb.goby.aligners;

import org.apache.commons.configuration.Configuration;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.exec.CommandLine;
import org.apache.commons.exec.DefaultExecutor;
import org.apache.commons.exec.PumpStreamHandler;
import org.apache.commons.lang.time.StopWatch;
import org.apache.commons.lang.SystemUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.log4j.Level;

import java.io.*;

import edu.cornell.med.icb.goby.modes.CompactToFastaMode;
import edu.cornell.med.icb.goby.modes.AbstractAlignmentToCompactMode;
import edu.cornell.med.icb.goby.util.LoggingOutputStream;
import edu.cornell.med.icb.goby.config.GobyConfiguration;
import edu.mssm.crover.cli.CLI;

/**
 * Wrapper around a modified GSNAP that supports Goby compact formats.
 *
 * @author Fabien Campagne
 *         Date: Oct 26, 2010
 *         Time: 4:00:02 PM
 */
public class GSnapAligner extends AbstractAligner {


    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(BWAAligner.class);
    private int numThreads;

    public CompactToFastaMode getReferenceCompactToFastaConverter() {
        return null;
    }

    public CompactToFastaMode getReadsCompactToFastaConverter() {
        return null;
    }

    public AbstractAlignmentToCompactMode getNativeAlignmentToCompactMode(String outputBasename) {
        return null;
    }

    public void setAlignerOptions(String options) {
        alignerOptions = "";
        if (options == null) {
            return;
        }
        final String[] opts = options.split("[,=]");

        // seedLength parameter
        if (CLI.isKeywordGiven(opts, "t")) {
            numThreads = CLI.getIntOption(opts, "t", 1);
        }
        alignerOptions += "-t " + numThreads;

    }


    public String getDefaultDbNameForReferenceFile(File referenceFile) {
        return FilenameUtils.getBaseName(referenceFile.getName());
    }

    public File prepareReads(File compactReadsFile) throws IOException {
        return compactReadsFile;
    }

    public File[] indexReference(File referenceFileOrDbBasename) throws IOException, InterruptedException {
        throw new UnsupportedOperationException("Indexing is not supported with GSNap at this time. ");
    }

    public void setConfiguration(final Configuration configuration) {
        pathToExecutables = configuration.getString(GobyConfiguration.EXECUTABLE_PATH_GSNAP, "");
    }

    private static final String GSNAP_EXE = SystemUtils.IS_OS_WINDOWS ? "gsnap.exe" : "gsnap";

    public File[] align(File referenceFile, File readsFile, String outputBasename) throws InterruptedException, IOException {
        final String alignCommand = FilenameUtils.concat(pathToExecutables, GSNAP_EXE);
        final CommandLine alignCommandLine = createCommandLine(alignCommand);

        if (colorSpace) {
            alignCommandLine.addArgument("--dibase");
        }

        alignCommandLine.addArguments(alignerOptions);
        alignCommandLine.addArgument("-A");
        alignCommandLine.addArgument("compact-alignment");
        alignCommandLine.addArgument("-b");
        alignCommandLine.addArgument(outputBasename);
        alignCommandLine.addArgument("-D");
        alignCommandLine.addArgument(databaseDirectory);
        alignCommandLine.addArgument("-d");
        alignCommandLine.addArgument(getDefaultDbNameForReferenceFile(referenceFile));
        alignCommandLine.addArgument(readsFile.getPath());

        LOG.info("About to execute " + alignCommandLine);
        final StopWatch timer = new StopWatch();
        timer.start();
        final DefaultExecutor alignExecutor = new DefaultExecutor();
        OutputStream saiOutputStream = null;
        OutputStream alignLogStream = null;
        try {
            saiOutputStream = new LoggingOutputStream(getClass(), Level.TRACE, "");
            alignLogStream = new LoggingOutputStream(getClass(), Level.INFO, "");
            alignExecutor.setStreamHandler(new PumpStreamHandler(saiOutputStream, alignLogStream));

            final int exitValue = alignExecutor.execute(alignCommandLine);
            LOG.info("Exit value = " + exitValue);
        } finally {
            IOUtils.closeQuietly(alignLogStream);
            IOUtils.closeQuietly(saiOutputStream);
        }
        LOG.info("Command executed in: " + timer.toString());

        return buildResults(outputBasename);
    }


}
