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

package edu.cornell.med.icb.goby.aligners;

import edu.cornell.med.icb.goby.config.GobyConfiguration;
import edu.cornell.med.icb.goby.modes.CompactToFastaMode;
import edu.cornell.med.icb.goby.util.LoggingOutputStream;
import org.apache.commons.configuration.Configuration;
import org.apache.commons.exec.CommandLine;
import org.apache.commons.exec.DefaultExecutor;
import org.apache.commons.exec.PumpStreamHandler;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang.SystemUtils;
import org.apache.commons.lang.time.StopWatch;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.log4j.Level;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

/**
 * This aligner uses Lastag.  Lastag is the name we have given our modified
 * version of the lastal alignment algorithm. Below is an attempt to gather
 * documentation for Lastag in one place, from the Last Manual and accompanying
 * documentation (manual.txt, tag-seeds.txt, ChangeLog.txt, ...).  Since Last
 * is a highly customizable tool which also involves various heuristic filters
 * for efficiency, these descriptions should be interpreted only as guides for
 * the functioning of the tool.
 * <p/>
 * Some sections from the LAST Manual are repeated here verbatim.  Further
 * documentation is included in: SVN/icb/3rdparty/last/doc/
 * <p/>
 * Alignments are computed in two steps.  First, LASTDB makes a database of the
 * reference sequence.  Then, for each sequence from a set of query sequences,
 * LASTAL searches the reference for high-scoring matching locations.  LASTAL
 * uses a two-step approach: first find initial matches, then extend alignments
 * from these matches.  In this case, the "initial matches" are: all exact
 * matches of any part of a tag to the genome, of depth >= [L], where the match
 * occurs at most [M] times in the genome.
 * <p/>
 * <h3> Steps in lastal </h3>
 * <p/>
 * <pre>
 * 1) Find initial matches:
 *      keep those with multiplicity <= [M] and depth >= [L].
 * 2) Extend gapless alignments from the initial matches:
 *      keep those with score >= [D]. (lastal recommends D = 3/5 * E)
 * 3) Extend gapped alignments from the gapless alignments:
 *      keep those with score >= [E].
 * 4) Non-redundantize the gapped alignments:
 *      remove those that share an endpoint with a higher-scoring alignment.
 * 5) Calculate probabilities (OFF by default).
 * 6) Redo the gapped extensions using centroid alignment (OFF by default).
 * </pre>
 * <p/>
 * <h3> Options for controlling Lastag </h3>
 * <p/>
 * <pre>
 * maxGapsAllowed - place limit on occurance of gaps in alignments
 * gapOpeningCost - specify cost to open a gap
 * m              - place limit on number of matches at highest score
 * matchQuality   - used to select [D] and [E] thresholds for efficiency
 * <p/>
 * XXX            - most flags that are defined for the Lastag command
 * <p/>
 * defaults are provided for many parameters but command-line values will override
 * </pre>
 * <p/>
 * <p/>
 * <h3> Merging identical tag sequences (recommended in tag-seeds.txt) </h3>
 * <p/>
 * If there are many identical sequences, we can speed up the mapping by
 * merging them.
 * [see ???]
 * <p/>
 * <h3> Discarding sub-optimal mappings (recommended in tag-seeds.txt) </h3>
 * <p/>
 * LAST will often align a tag to more than one genome location.  We may wish
 * to keep only the highest-scoring alignment(s) for each tag.
 * [see LastToCompactMode]
 * <p/>
 * <h3> Discarding ambiguous tags (recommended in tag-seeds.txt) </h3>
 * <p/>
 * There may be some tags that map to more than one location (with equal
 * scores).  We may wish to discard such multi-mapping tags.
 * [see LastToCompactMode]
 * <p/>
 * <h3> Score function </h3>
 * <p/>
 * An alignment is scored as follows:
 * <pre>
 *  score = R*num_matches + Q*num_mismatches + SUM_i gap_score_i
 *  gap_score = A + B*k, for a gap of size k
 * </pre>
 * where:
 * <pre>
 *  R = match score
 *  Q = mismatch score
 *  A = gap existence cost
 *  B = gap extension cost
 * </pre>
 * <p/>
 * <h3> Score function for colorspace </h3>
 * <p/>
 * An alignment is scored as follows:
 * <pre>
 *  score = R*num_matches + Q*num_mismatches + SUM_i gap_score_i
 *  gap_score = A + B*k, for a gap of size k
 * </pre>
 * where:
 * <pre>
 *  R = match score
 *  Q = mismatch score
 *  A = gap existence cost
 *  B = gap extension cost
 * </pre>
 * <p/>
 * <h3> lastdb: usage: lastdb [options] output-name fasta-sequence-file(s) </h3>
 * <p/>
 * <pre>
 * Main Options (default settings):
 * -p: interpret the sequences as proteins
 * -c: read the sequences case-sensitively
 * -m: periodic spaced-seed pattern (1)
 * <p/>
 * Advanced Options (default settings):
 * -w: index step (1)
 * -s: volume size (1342177280)
 * -a: user-defined alphabet
 * -b: bucket depth
 * -v: be verbose: write messages about what lastdb is doing
 * </pre>
 * <p/>
 * <h3> lastag: usage: lastal [options] lastdb-name fasta-sequence-file(s) </h3>
 * <p/>
 * <pre>
 * Main options (default settings):
 * -h: show all options and their default settings
 * -o: output file
 * -u: mask lowercase letters: 0=off, 1=softer, 2=soft, 3=hard (0)
 * -s: strand: 0=reverse, 1=forward, 2=both (2 for DNA, 1 for protein)
 * -f: output format: 0=tabular, 1=maf (1)
 * <p/>
 * Score parameters (default settings):
 * -r: match score   (1 for DNA, blosum62 for protein)
 * -q: mismatch cost (1 for DNA, blosum62 for protein)
 * -p: file for residue pair scores
 * -a: gap existence cost (7 for DNA, 11 for protein)
 * -b: gap extension cost (1 for DNA,  2 for protein)
 * -c: unaligned residue pair cost (100000)
 * -x: maximum score dropoff for gapped extensions (max[y, a+b*20])
 * -y: maximum score dropoff for gapless extensions (max-match-score * 10)
 * -d: minimum score for gapless alignments (e*3/5)
 * -e: minimum score for gapped alignments (40 for DNA, 100 for protein)
 * <p/>
 * Miscellaneous options (default settings):
 * -m: maximum multiplicity for initial matches (10)
 * -l: minimum depth for initial matches (1)
 * -k: step-size along the query sequence (1)
 * -i: query batch size (16 MiB when counting matches, 128 MiB otherwise)
 * -w: supress repeats within this distance inside large exact matches (1000)
 * -t: 'temperature' for calculating probabilities (1/lambda)
 * -g: 'gamma' parameter for gamma-centroid alignment (1)
 * -v: be verbose: write messages about what lastal is doing
 * -j: output type: 0=match counts, 1=gapless, 2=redundant gapped, 3=gapped,
 *                  4=probabilities, 5=centroid (3)
 * -z: write counts using a more compact single line format
 * </pre>
 * <p/>
 * <h3> Miscellaneous </h3>
 * <p/>
 * Depth = the number of matched, non-skipped nucleotides
 * <p/>
 * The Last Manual describes options for using quality scores (manual.txt and
 * tag-seeds.txt), although these options restrict the flexibility one has to
 * select match & mismatch scores
 *
 * @author Fabien Campagne
 *         Date: Jul 9, 2009
 *         Time: 11:39:29 AM
 */
public class LastagAligner extends LastAligner {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(LastagAligner.class);

    private static final String LASTAG_EXEC = SystemUtils.IS_OS_WINDOWS ? "lastag.exe" : "lastag";

    public LastagAligner() {
        super();
    }

    @Override
    public void setConfiguration(final Configuration configuration) {
        pathToExecutables = configuration.getString(GobyConfiguration.EXECUTABLE_PATH_LASTAG, "");
    }

    @Override
    public CompactToFastaMode getReadsCompactToFastaConverter() {
        final CompactToFastaMode toFastaConverter = super.getReadsCompactToFastaConverter();
        // for lastag, force fasta format:
        toFastaConverter.setOutputFormat(CompactToFastaMode.OutputFormat.FASTA);
        return toFastaConverter;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public File[] align(final File referenceFile, final File readsFile, final String outputBasename)
            throws IOException {
        assert pathToExecutables != null : "path to executables must be defined.";
        validateScoreOptions();

        final File nativeReads = prepareReads(readsFile);
        if (databaseName == null) {
            databaseName = getDefaultDbNameForReferenceFile(referenceFile);
        }
        indexReference(referenceFile);
        LOG.info("Searching..");
        assert outputBasename != null : "basename must not be null";
        forceMakeParentDir(outputBasename);
        // if (colorSpace) { prepareMatrixFile(); }

        // full path to lastag executable
        final String command = FilenameUtils.concat(pathToExecutables, LASTAG_EXEC);
        final CommandLine commandLine = createCommandLine(command);
        commandLine.addArgument("-a" + gapOpeningCost);
        commandLine.addArgument("-b" + gapExtentionCost);
        commandLine.addArgument("-d" + dThresholdOption(minReadLength));
        commandLine.addArgument("-e" + eThresholdOption(minReadLength));
        commandLine.addArgument("-m" + getSeedMaxMultiplicity());
        commandLine.addArgument("-l" + minSeedDepth);
        commandLine.addArgument("-s" + DEFAULT_STRAND_DIRECTION);
        commandLine.addArgument("-j" + DEFAULT_OUTPUT_TYPE);

        // this option is a lastag extension which produces a "counts" file as part of the
        // alignment process. Also, the counts are written in a different "compressed" format.
        commandLine.addArgument("-z");

        commandLine.addArgument("-i" + DEFAULT_QUERY_BATCH_SIZE);
        commandLine.addArgument("-v");
        commandLine.addArgument("-f" + DEFAULT_OUTPUT_FORMAT);
        commandLine.addArguments(alignerOptions, false);       // don't quote these options
        commandLine.addArguments(scoreOptions(), false);       // don't quote these options
        commandLine.addArgument("-o");
        commandLine.addArgument(outputBasename);
        commandLine.addArgument(getDatabasePath(databaseDirectory, databaseName));
        commandLine.addArgument(nativeReads.toString());

        LOG.info("About to execute " + commandLine);
        final StopWatch timer = new StopWatch();
        timer.start();
        final DefaultExecutor executor = new DefaultExecutor();
        OutputStream logStream = null;
        try {
            logStream = new LoggingOutputStream(getClass(), Level.INFO, "");
            executor.setStreamHandler(new PumpStreamHandler(logStream));

            final int exitValue = executor.execute(commandLine);
            LOG.info("Exit value = " + exitValue);
        } finally {
            IOUtils.closeQuietly(logStream);
        }
        LOG.info("Command executed in: " + timer.toString());

        // convert native alignment into compact reads
        final File[] buildResults = processAlignment(referenceFile, readsFile, outputBasename);
        if (!keepTemporaryFiles) {
            FileUtils.deleteQuietly(new File(outputBasename + ".maf"));
            FileUtils.deleteQuietly(new File(outputBasename + ".counts"));

        }
        // also delete the colorSpaceMatrix file, if created
        FileUtils.deleteQuietly(new File(matrixFilename));
        return buildResults;
    }
}
