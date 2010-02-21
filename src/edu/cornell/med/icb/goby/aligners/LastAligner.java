/*
 * Copyright (C) 2010 Institute for Computational Biomedicine,
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

import edu.cornell.med.icb.goby.config.GobyConfiguration;
import edu.cornell.med.icb.goby.modes.AbstractAlignmentToCompactMode;
import edu.cornell.med.icb.goby.modes.CompactToFastaMode;
import edu.cornell.med.icb.goby.modes.LastToCompactMode;
import edu.cornell.med.icb.goby.modes.FastaToCompactMode;
import edu.cornell.med.icb.goby.reads.ColorSpaceConverter;
import edu.cornell.med.icb.goby.util.LoggingOutputStream;
import edu.mssm.crover.cli.CLI;
import org.apache.commons.configuration.Configuration;
import org.apache.commons.exec.CommandLine;
import org.apache.commons.exec.DefaultExecutor;
import org.apache.commons.exec.PumpStreamHandler;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.io.LineIterator;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang.SystemUtils;
import org.apache.commons.lang.time.StopWatch;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.log4j.Level;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.util.LinkedList;
import java.util.List;

/**
 * This aligner uses <a href="http://last.cbrc.jp/">Last</a>.
 * <p/>
 * Alignments are computed in two steps.  First, LASTDB makes a database of the
 * reference sequence.  Then, for each sequence from a set of query sequences,
 * LASTAL searches the reference for high-scoring matching locations.  LASTAL
 * uses a two-step approach: first find initial matches, then extend alignments
 * from these matches.  In this case, the "initial matches" are: all exact
 * matches of any part of a tag to the genome, of depth >= [L], where the match
 * occurs at most [M] times in the genome.
 *
 * <h3> Steps in lastal </h3>
 *
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
 *
 * <h3> Score function </h3>
 *
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
 *
 * <h3> Score function for colorspace </h3>
 *
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
 *
 * <h3> lastdb: usage: lastdb [options] output-name fasta-sequence-file(s) </h3>
 *
 * <pre>
 * Main Options (default settings):
 * -p: interpret the sequences as proteins
 * -c: read the sequences case-sensitively
 * -m: periodic spaced-seed pattern (1)
 *
 * Advanced Options (default settings):
 * -w: index step (1)
 * -s: volume size (1342177280)
 * -a: user-defined alphabet
 * -b: bucket depth
 * -v: be verbose: write messages about what lastdb is doing
 * </pre>
 *
 * <h3> lastag: usage: lastal [options] lastdb-name fasta-sequence-file(s) </h3>
 *
 * <pre>
 * Main options (default settings):
 * -h: show all options and their default settings
 * -o: output file
 * -u: mask lowercase letters: 0=off, 1=softer, 2=soft, 3=hard (0)
 * -s: strand: 0=reverse, 1=forward, 2=both (2 for DNA, 1 for protein)
 * -f: output format: 0=tabular, 1=maf (1)
 *
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
 *
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
 * </pre>
 *
 * <h3> Miscellaneous </h3>
 *
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
// TODO replace getSeedMaxMultiplicity() with setSeedMaxMultiplicity()
public class LastAligner extends AbstractAligner {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(LastAligner.class);

    protected enum MatchQuality {
        PERFECT_MATCHES_ONLY,
        ONE_MISMATCH,
        TWO_MISMATCHES,
        IDENTITY_95,
        BEST_MATCH
    }

    protected static final String LASTDB_EXEC = SystemUtils.IS_OS_WINDOWS ? "lastdb.exe" : "lastdb";
    private static final String LASTAL_EXEC = SystemUtils.IS_OS_WINDOWS ? "lastal.exe" : "lastal";

    // parameters that control Lastdb and/or Lastal indirectly
    protected static final String LASTDB_DEFAULT_MEMORY     = "2G"; // -s:
    protected static final int DEFAULT_MAX_GAPS_ALLOWED     = 0;
    protected static final int MIN_READ_LENGTH_APPROX       = 20;
    protected static final String DEFAULT_MATCH_QUALITY     = "IDENTITY_95";

    // fixed values for Lastag command line - setAlignerOptions does not change these
    protected static final int DEFAULT_STRAND_DIRECTION     = 2; // -s: strand: 0 = reverse, 1 = forward, 2 = both (2 for DNA, 1 for protein)
    protected static final int DEFAULT_OUTPUT_TYPE          = 3; // -j: output type: 0 = match counts, 1 = gapless, 2 = redundant gapped, 3 = gapped, 4 = probabilities, 5 = centroid
    protected static final int DEFAULT_OUTPUT_FORMAT        = 1; // -f: output format: 0 = tabular, 1 = maf
    protected static final String DEFAULT_QUERY_BATCH_SIZE  = "1G"; // -i: query batch size

    // suggested defaults for Lastal command line - setAlignerOptions may change these
    protected static final int DEFAULT_GAP_OPENING_COST     = 2; // -a
    protected static final int DEFAULT_GAP_EXTENSION_COST   = 1; // -b
    protected static final int DEFAULT_MIN_SEED_DEPTH       = 15; // -l
    protected static final int DEFAULT_MISMATCH_COST        = 1;
    protected static final int DEFAULT_MATCH_TO_MISMATCH_RATIO_NT = 1;
    protected static final int DEFAULT_MATCH_TO_MISMATCH_RATIO_CS = 2;

    // parameters to resolve Lastal command-line options
    protected int maxGapsAllowed = DEFAULT_MAX_GAPS_ALLOWED;
    protected int gapOpeningCost = DEFAULT_GAP_OPENING_COST;
    protected int gapExtentionCost = DEFAULT_GAP_EXTENSION_COST;
    protected int minSeedDepth = DEFAULT_MIN_SEED_DEPTH;
    protected int matchScore = -1;   // -1 indicates that value has not been set by user
    protected int misMatchCost = -1; // -1 indicates that value has not been set by user

    protected MatchQuality matchQuality = MatchQuality.valueOf(DEFAULT_MATCH_QUALITY);

    // Can not yet initialize matrixFilename, since relevant members not set (e.g. work directory)
    protected String matrixFilename = "";

    public LastAligner() {
        super();
        extensions = new String[]{"prj"};
    }

    public void setConfiguration(final Configuration configuration) {
        pathToExecutables = configuration.getString(GobyConfiguration.EXECUTABLE_PATH_LAST, "");
    }

    public String getDefaultDbNameForReferenceFile(final File referenceFile) {
        return FilenameUtils.getBaseName(referenceFile.getName()).replaceAll(" ", "_");
    }

    /**
     * Return the reference file converter used by aligner algorithm.
     */
    @Override
    public CompactToFastaMode getReferenceCompactToFastaConverter() {
        final CompactToFastaMode processor = new CompactToFastaMode();
        processor.setHashOutputFilename(true);
        processor.setIndexToHeader(false);
        // Reference conversion (always from nt-space) *MAY* be needed by alignment algorithm to match colorspace reads platforms
        processor.setOutputColorMode(colorSpace); // lastal requires this conversion when dealing with colorSpace
        processor.setOutputFakeNtMode(false); // lastal uses real colorSpace
        processor.setOutputFakeQualityMode(false); // lastal currently does not use (fake) qualities
        // Filter using the default alphabet for specified output mode
        processor.setAlphabet((colorSpace) ? "0123ACTG" : "ACTG"); // lastal expects colorspace with nt-adaptors, or just nt
        return processor;
    }


    /**
     * Return the reads file converter used by aligner algorithm.
     */
    @Override
    public CompactToFastaMode getReadsCompactToFastaConverter() {
        final CompactToFastaMode processor = new CompactToFastaMode();
        processor.setIndexToHeader(true);
        processor.setHashOutputFilename(true);
        // Since colorSpace is determined by reads platform, colorspace conversion is never needed for reads
        processor.setOutputFakeNtMode(false); // lastal uses real colorSpace
        processor.setOutputFakeQualityMode(false); // lastal currently does not use (fake) qualities
        // Filter using the default alphabet for specified output mode
        processor.setAlphabet((colorSpace) ? "0123ACTG" : "ACTG"); // lastal expects colorspace with nt-adaptors, or just nt
        processor.setTrimAdaptorLength((colorSpace) ? 2 : 0); // prevent aligner from finding false matches with adaptor at the start of the sequence
        processor.setOutputFormat(CompactToFastaMode.OutputFormat.FASTQ);   // We use Fastq format with lastal.
        processor.setQualityEncoding(FastaToCompactMode.QualityEncoding.SANGER);   // We use Fastq format with Sanger encoding for LASTAL.
        return processor;
    }

    /**
     * Returns LastToCompact processor, initialized with correct input file.
     */
    @Override
    public AbstractAlignmentToCompactMode getNativeAlignmentToCompactMode(final String outputBasename) {
        // import counts & entries
        final LastToCompactMode processor = new LastToCompactMode();
        processor.setInputFile(outputBasename + ".maf");
        processor.setOnlyMafFile(false);
        processor.setOnlyCountsFile(false);
        return processor;
    }

    /**
     * Parse and validate aligner specific options.
     * Input options are comma separated, with the syntax key1=value1,key2=value2.
     * "alignerOptions" format should match aligner's command-line specification e.g. "-key1 value1 -key2 value2"
     * This method is declared abstract so that aligners have control over the parsing of their arguments.
     */
    @Override
    public void setAlignerOptions(final String options) {
        maxGapsAllowed = DEFAULT_MAX_GAPS_ALLOWED;
        gapOpeningCost = DEFAULT_GAP_OPENING_COST;
        matchQuality   = MatchQuality.valueOf(DEFAULT_MATCH_QUALITY);
        alignerOptions = "-Q 1 ";   // indicate FASTQ input in SANGER format.
        if (options == null) {
            return;
        }
        final String[] opts = options.split("[,=]");

        // options specific to LastalAligner
        maxGapsAllowed = CLI.getIntOption(opts, "maxGapsAllowed", DEFAULT_MAX_GAPS_ALLOWED);
        matchQuality = MatchQuality.valueOf(CLI.getOption(opts, "matchQuality", DEFAULT_MATCH_QUALITY));
        gapExtentionCost = CLI.getIntOption(opts, "b", DEFAULT_GAP_EXTENSION_COST);
        minSeedDepth = CLI.getIntOption(opts, "l", DEFAULT_MIN_SEED_DEPTH);


        // always favor --gapOpeningCost over redundant option -a
        if (CLI.isKeywordGiven(opts, "a")) {
            gapOpeningCost = CLI.getIntOption(opts, "a", DEFAULT_GAP_OPENING_COST);
        }
        if (CLI.isKeywordGiven(opts, "gapOpeningCost")) {
            gapOpeningCost = CLI.getIntOption(opts, "gapOpeningCost", DEFAULT_GAP_OPENING_COST);
        }

        // defaults for match score and mismatch cost are determined w.r.t.
        // colorspace but command-line values will override these defaults
        if (CLI.isKeywordGiven(opts, "r")) {
            matchScore = CLI.getIntOption(opts, "r", -999);
            assert (matchScore != -999) : "Parsing error";
        }
        if (CLI.isKeywordGiven(opts, "q")) {
            misMatchCost = CLI.getIntOption(opts, "q", -999);
            assert (misMatchCost != -999) : "Parsing error";
        }

        // options specific to "bwa aln" tool
        // FLOAT -t -g
        alignerOptions += (CLI.isKeywordGiven(opts, "t")) ? " -t " + CLI.getFloatOption(opts, "t", -999) : "";
        alignerOptions += (CLI.isKeywordGiven(opts, "g")) ? " -g " + CLI.getFloatOption(opts, "g", -999) : "";
        // INT -c -x -y -a -m -l -k -w
        alignerOptions += (CLI.isKeywordGiven(opts, "c")) ? " -c " + CLI.getIntOption(opts, "c", -999) : "";
        alignerOptions += (CLI.isKeywordGiven(opts, "x")) ? " -x " + CLI.getIntOption(opts, "x", -999) : "";
        alignerOptions += (CLI.isKeywordGiven(opts, "y")) ? " -y " + CLI.getIntOption(opts, "y", -999) : "";
        // alignerOptions += (CLI.isKeywordGiven(opts, "m")) ? " -m " + CLI.getIntOption(opts, "m", -999) : ""; // setAmbiguityThreshold()
        alignerOptions += (CLI.isKeywordGiven(opts, "k")) ? " -k " + CLI.getIntOption(opts, "k", -999) : "";
        alignerOptions += (CLI.isKeywordGiven(opts, "w")) ? " -w " + CLI.getIntOption(opts, "w", -999) : "";
        // BOOLEAN -v
        alignerOptions += (CLI.isKeywordGiven(opts, "v")) ? " -v " : "";
        // EXCLUDED -f -u -s -p -i -j -z -d -e
        // see : DEFAULT_STRAND_DIRECTION, DEFAULT_OUTPUT_TYPE,
        // DEFAULT_QUERY_BATCH_SIZE, DEFAULT_OUTPUT_FORMAT, ...
        assert (!CLI.isKeywordGiven(opts, "f")) : "Currently not supporting  -f  option via alignerOptions";
        assert (!CLI.isKeywordGiven(opts, "u")) : "Currently not supporting  -u  option via alignerOptions";
        assert (!CLI.isKeywordGiven(opts, "s")) : "Currently not supporting  -s  option via alignerOptions";
        assert (!CLI.isKeywordGiven(opts, "p")) : "Currently not supporting  -p  option via alignerOptions";
        assert (!CLI.isKeywordGiven(opts, "i")) : "Currently not supporting  -i  option via alignerOptions";
        assert (!CLI.isKeywordGiven(opts, "j")) : "Currently not supporting  -j  option via alignerOptions";
        assert (!CLI.isKeywordGiven(opts, "z")) : "Currently not supporting  -z  option via alignerOptions";
        assert (!CLI.isKeywordGiven(opts, "d")) : "Currently not supporting  -d  option via alignerOptions";
        assert (!CLI.isKeywordGiven(opts, "e")) : "Currently not supporting  -e  option via alignerOptions";
        //
        assert (!alignerOptions.contains("-999")) : "Parsing error.";
        validateScoreOptions();
    }

    /**
     * Specifies matchQuality, which can be one of:
     * <p/>
     * <LI> PERFECT_MATCHES_ONLY,
     * <LI> ONE_MISMATCH,
     * <LI> TWO_MISMATCHES,
     * <LI> IDENTITY_95,
     * <LI> BEST_MATCH
     *
     * @param matchQuality
     */
    public void setMatchQuality(final String matchQuality) {
        this.matchQuality = MatchQuality.valueOf(matchQuality);
    }

    /**
     * The match/mismatch substitutionMatrix for residue pair scores is defined
     * in ColorSpaceConverter.getColorSpaceSubstitutionMatrix() as 2 for a
     * matching pairs and -1 for all mismatches
     * <p/>
     * This relative value of matches / mismatches is appropriate when scoring
     * alignments for colorspace reads. E.g.:
     * <p/>
     * 1) Since individual colorspace errors are likely read errors (a property
     * of 2-base encoding sequencers), the smaller relative penalty for
     * mismatches (-1 versus 2) makes sense.  A more complicated score function
     * (one with more parameters) would likely not penalize these at all.
     * <p/>
     * 2) Biological sequence differences between the query and the reference,
     * are most likely SNPS, which in colorspace affect the two adjacent colors
     * (2 mismatches).  Similarly, correctly sequenced and aligned bases most
     * often occur consecutively, and therefore guarantee k-l matching colors
     * for every k consecutive matching nucleotides.  Thus the score function
     * will often be penalized twice the mismatch score, while matches will
     * only accumulate a single match reward.  The 2:-1 match/mismatch residue
     * pair scoring is an attempt to correct for this.
     */
    protected void validateScoreOptions() {
        if (misMatchCost < 0) {
            misMatchCost = DEFAULT_MISMATCH_COST;
        }
        if (matchScore < 0) {
            final int ratio = (colorSpace) ? DEFAULT_MATCH_TO_MISMATCH_RATIO_CS : DEFAULT_MATCH_TO_MISMATCH_RATIO_NT;
            matchScore = ratio * misMatchCost;
        }
        if (mParameter < 1) {
            mParameter = 1; // must be at least 1
        }
    }

    protected String scoreOptions() {
        validateScoreOptions();
        return "-r " + matchScore + " -q " + misMatchCost;
    }

    /**
     * Currently disabled so that user can specify
     *  -r <match-score> -q <mismatch-score> on the command line
     * parsed through setAlignerOptions. See also scoreOptions()
     */
    private String substitutionMatrixOption() {
        // prepare colorSpaceMatrix file, if needed
        if (colorSpace) {
            return "-p " + matrixFilename;
        } else {
            return " ";
        }
    }

    public void prepareMatrixFile() throws IOException {
        matrixFilename = FilenameUtils.concat(workDirectory, "color-space-matrix.mat");
        forceMakeParentDir(matrixFilename);
        FileWriter writer = null;
        try {
            writer = new FileWriter(matrixFilename);
            writer.write(ColorSpaceConverter.getColorSpaceSubstitutionMatrix());
        } finally {
            IOUtils.closeQuietly(writer);
        }
    }

    /**
     * Initialize lastal's seedMaxMultiplicty in terms of our mParameter. Lastal will only
     * extend seeds that are unambigous, according to this value.
     * <p/>
     * Typically: seedMaxMultiplicity >= mParameter
     * <p/>
     * Rationale: For each read, lastal extracts contiguous subsequences (seeds) that match at
     * most this number of times to the reference genome.  The initial seed alignments are
     * used to evaluate full read-length alignments at each location.  After scoring the
     * full-length alignments, the best ones are returned.  The number returned could be
     * less than seedMaxMultiplicity.
     * <p/>
     * Since this algorithm uses heuristics to seed and filter alignments, it is difficult to
     * characterize the results.
     * <p/>
     * Our strategy is to tune the lastag algorithm so that it will generate lots of results,
     * and then by filtering them, we can *try to* guarantee certain properties of the output.
     * <p/>
     * To maximize the number of alignments that match the genome at most mParameter times,
     * we select seedMaxMultiplicity greater than mParameter.
     */
    protected int getSeedMaxMultiplicity() {
        validateScoreOptions();
        return 2 * mParameter;
    }

    /**
     * Return minimum score for gapped alignments, derived from number of
     * specified matches & mismatches, and the minimum read length. This
     * parameter is used to filter seeds.  Once the number of matches /
     * mismatches is determined
     *
     * The matchQuality parameter is used to select allowable mismatches.
     * The minimumReadLength - numMisMatch is used as a lower bound on the
     * number of matches.
     *
     * The threshold is the minimum of the score function, computed under the
     * assumption that
     * 1) all mismatches occur in snps
     * 2) all mismatches occur in maxGapsAllowed gaps
     *
     * This assumes that the reads will have roughly the same length. If this
     * is not the case (i.e. for Helicos data, which has a triangular
     * distribution over read lengths), the resulting thresholds will be very
     * loose, allowing most seeds to qualify for subsequent full-length
     * alignment.  This will take a very long time.  To retain efficiency, one
     * should partition the set of reads by length.
     * <p/>
     * Parameters and score function described in the Last Manual
     * <pre>
     *  score = R*num_matches + Q*num_mismatches + SUM_i gap_score_i
     *  gap_score = A + B*k, for gap of size k
     * </pre>
     * where:
     * <pre>
     *  R = match score
     *  Q = mismatch score
     *  A = gap existence cost
     *  B = gap extension cost
     * </pre>
     * @param minReadLength
     * @return
     */
    protected int eThresholdOption(int minReadLength) {
        validateScoreOptions();
        if (minReadLength == 0) {
            minReadLength = MIN_READ_LENGTH_APPROX;
        }
        if (colorSpace) {
            minReadLength = minReadLength - 1;
        }
        final int numMatch;
        final int numMisMatch;
        switch (matchQuality) {
            case PERFECT_MATCHES_ONLY:
                numMatch = minReadLength;
                numMisMatch = 0;
                break;
            case ONE_MISMATCH:
                numMatch = minReadLength - 1;
                numMisMatch = 1;
                break;
            case TWO_MISMATCHES:
                numMatch = minReadLength - 2;
                numMisMatch = 2;
                break;
            case IDENTITY_95:
                numMisMatch = (int) Math.ceil(minReadLength * 0.05);
                numMatch = minReadLength - numMisMatch;
                break;
            case BEST_MATCH:
                numMatch = minReadLength / 3;
                numMisMatch = 0;
                // ensures threshold is low enough to search all potential candidates
                break;
            default:
                numMisMatch = (int) Math.floor(minReadLength * 0.10);
                numMatch = minReadLength - numMisMatch;
        }
        final int eThreshAssumingSNPs = numMatch * matchScore - numMisMatch * misMatchCost;
        final int eThreshAssumingGAPs = numMatch * matchScore
                - numMisMatch * gapExtentionCost - maxGapsAllowed * gapOpeningCost;
        return Math.max(0, Math.min(eThreshAssumingSNPs, eThreshAssumingGAPs));
    }

    /**
     * Value used by last v52. 3/5 of the gapped score should be achieved by the gapless extension.
     */
    protected int dThresholdOption(final int minReadLength) {
        validateScoreOptions();
        if (maxGapsAllowed == 0) {
            return eThresholdOption(minReadLength);
        } else {
            return eThresholdOption(minReadLength) * 3 / 5;
        }
    }

    /**
     * If referenceFileOrDbBasename is a reference compact (compact-reads or fasta or fa)
     * this will create a database based on database directory and database name.
     * If referenceFileOrDbBasename is a db basename, it must already exist.
     *
     * @param referenceFileOrDbBasename The compact file with the reference sequence to index OR
     * the database basename (prefix)
     * @return
     * @throws IOException
     */
    public File[] indexReference(final File referenceFileOrDbBasename) throws IOException {
        if (isDatabaseBasename(referenceFileOrDbBasename.toString())) {
            return matchingExtension(referenceFileOrDbBasename.toString(), extensions);
        }

        final String databasePathPrefix = getDatabasePath(databaseDirectory, databaseName);
        if (isDatabaseBasename(databasePathPrefix)) {
            return matchingExtension(databasePathPrefix, extensions);
        }
        final File fastaReferenceFile = prepareReference(referenceFileOrDbBasename);
        forceMakeParentDir(databasePathPrefix);

        // full path to lastdb executable
        final String command = FilenameUtils.concat(pathToExecutables, LASTDB_EXEC);
        final CommandLine commandLine = createCommandLine(command);
        commandLine.addArgument("-v");
        commandLine.addArgument("-s" + LASTDB_DEFAULT_MEMORY);
        commandLine.addArgument("-a");
        commandLine.addArgument(getAlphabet());
        commandLine.addArgument(databasePathPrefix);
        commandLine.addArgument(fastaReferenceFile.toString());

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
        listFiles(databaseDirectory);
        return matchingExtension(databasePathPrefix, extensions);
    }

    /**
     * Align.
     *
     * @param referenceFile  Compact or native database basename for reference sequences.
     * @param readsFile      Compact or native format (i.e., fasta, fastq) aligner read format.
     * @param outputBasename Basename where to write the compact alignment.
     * @return
     * @throws IOException
     */
    public File[] align(final File referenceFile, final File readsFile, final String outputBasename) throws IOException {
        assert pathToExecutables != null : "path to executables must be defined.";
        validateScoreOptions();

        final File nativeReads = prepareReads(readsFile);
        if (databaseName == null) {
            databaseName = getDefaultDbNameForReferenceFile(referenceFile);
        }
        final File[] indexedReference = indexReference(referenceFile);
        LOG.info("Searching..");
        forceMakeParentDir(outputBasename);
        // if (colorSpace) { prepareMatrixFile(); }
        align(outputBasename, nativeReads);
        count(outputBasename, nativeReads);

        // convert native alignment into compact reads
        final File[] buildResults = processAlignment(referenceFile, readsFile, outputBasename);

        FileUtils.deleteQuietly(new File(outputBasename + ".maf"));
        FileUtils.deleteQuietly(new File(outputBasename + ".counts"));

        // also delete the colorSpaceMatrix file, if created
        FileUtils.deleteQuietly(new File(matrixFilename));

        return buildResults;
    }

    private void align(final String outputBasename, final File nativeReads) throws IOException {
        // full path to lastal executable
        final String command = FilenameUtils.concat(pathToExecutables, LASTAL_EXEC);
        final CommandLine commandLine = createCommandLine(command);
        commandLine.addArgument("-a" + gapOpeningCost);
        commandLine.addArgument("-b" + gapExtentionCost);
        commandLine.addArgument("-d" + dThresholdOption(minReadLength));
        commandLine.addArgument("-e" + eThresholdOption(minReadLength));
        commandLine.addArgument("-m" + getSeedMaxMultiplicity());
        commandLine.addArgument("-l" + minSeedDepth);
        commandLine.addArgument("-s" + DEFAULT_STRAND_DIRECTION);
        commandLine.addArgument("-j" + DEFAULT_OUTPUT_TYPE);
        commandLine.addArgument("-i" + DEFAULT_QUERY_BATCH_SIZE);
        commandLine.addArgument("-v");
        commandLine.addArgument("-f" + DEFAULT_OUTPUT_FORMAT);
        commandLine.addArguments(alignerOptions, false);       // don't quote these options
        commandLine.addArguments(scoreOptions(), false);       // don't quote these options
        commandLine.addArgument("-o");
        commandLine.addArgument(outputBasename + ".maf");
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
    }

    private void count(final String outputBasename, final File nativeReads) throws IOException {
        // full path to lastal executable
        final String command = FilenameUtils.concat(pathToExecutables, LASTAL_EXEC);
        final CommandLine commandLine = createCommandLine(command);
        commandLine.addArgument("-a" + gapOpeningCost);
        commandLine.addArgument("-b" + gapExtentionCost);
        commandLine.addArgument("-d" + dThresholdOption(minReadLength));
        commandLine.addArgument("-e" + eThresholdOption(minReadLength));
        commandLine.addArgument("-m" + getSeedMaxMultiplicity());
        commandLine.addArgument("-l" + minSeedDepth);
        commandLine.addArgument("-s" + DEFAULT_STRAND_DIRECTION);
        commandLine.addArgument("-j0");
        commandLine.addArgument("-i" + DEFAULT_QUERY_BATCH_SIZE);
        commandLine.addArgument("-v");
        commandLine.addArgument("-f" + DEFAULT_OUTPUT_FORMAT);
        commandLine.addArguments(alignerOptions, false);       // don't quote these options
        commandLine.addArguments(scoreOptions(), false);       // don't quote these options
        commandLine.addArgument("-o");
        final String countsFilename = outputBasename + ".counts";
        commandLine.addArgument(countsFilename);
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

        compressCounts(countsFilename);
    }

    private void compressCounts(final String countsFilename) throws IOException {
        final List<String> linesOut = new LinkedList<String>();
        final LineIterator lineIterator = FileUtils.lineIterator(new File(countsFilename));
        int queryName = 0;
        int length = 0;
        int count = 0;
        try {
            while (lineIterator.hasNext()) {
                final String line = lineIterator.nextLine();
                if (line.length() > 0 && line.charAt(0) == '#') {
                    if (line.equals("# length\tcount")) {
                        linesOut.add("query-name\tdepth\tquery-length\tnumber-of-matches");
                    } else {
                        linesOut.add(line);
                    }
                } else {
                    if (StringUtils.isBlank(line)) {
                        // write the previous count (if we have one)
                        if (count != 0) {
                            linesOut.add(buildCountLine(queryName, length, count));
                        }
                        // reset counts
                        queryName = 0;
                        length = 0;
                        count = 0;
                    } else {
                        final String[] tokens = StringUtils.split(line);
                        if (tokens.length == 1) {
                            queryName = Integer.parseInt(tokens[0]);
                        } else if (tokens.length == 2) {
                            length = Integer.parseInt(tokens[0]);
                            count = Integer.parseInt(tokens[1]);
                        } else {
                            LOG.warn("Ignoring line: " + line);
                        }
                    }
                }
            }
        } finally {
            LineIterator.closeQuietly(lineIterator);
        }

        // write the last count (if we have one)
        if (count != 0) {
            linesOut.add(buildCountLine(queryName, length, count));
        }

        // and write the "compressed" counts to the file
        FileUtils.writeLines(new File(countsFilename), linesOut);
    }

    private String buildCountLine(final int queryName, final int length, final int count) {
        final StringBuilder countLine = new StringBuilder(80);
        countLine.append(queryName);
        countLine.append("\t");
        countLine.append(length);
        countLine.append("\t");
        countLine.append(length);
        countLine.append("\t");
        countLine.append(count);
        return countLine.toString();
    }
}
