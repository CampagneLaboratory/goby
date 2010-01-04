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

package edu.cornell.med.icb.goby.aligners;

import edu.cornell.med.icb.goby.modes.AbstractAlignmentToCompactMode;
import edu.cornell.med.icb.goby.modes.CompactToFastaMode;
import edu.cornell.med.icb.goby.modes.LastToCompactMode;
import edu.cornell.med.icb.reads.ColorSpaceConverter;
import edu.cornell.med.icb.util.ExecuteProgram;
import edu.cornell.med.icb.util.GobyPropertyKeys;
import edu.cornell.med.icb.util.GroovyProperties;
import edu.mssm.crover.cli.CLI;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang.SystemUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.log4j.Level;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

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
 *
 *
 * <h3> Options for controlling Lastag </h3>
 *
 * <pre>
 * maxGapsAllowed - place limit on occurance of gaps in alignments
 * gapOpeningCost - specify cost to open a gap
 * m              - place limit on number of matches at highest score
 * matchQuality   - used to select [D] and [E] thresholds for efficiency
 *
 * XXX            - most flags that are defined for the Lastag command
 *
 * defaults are provided for many parameters but command-line values will override
 * </pre>
 *
 *
 * <h3> Merging identical tag sequences (recommended in tag-seeds.txt) </h3>
 *
 * If there are many identical sequences, we can speed up the mapping by
 * merging them.
 * [see ???]
 *
 * <h3> Discarding sub-optimal mappings (recommended in tag-seeds.txt) </h3>
 *
 * LAST will often align a tag to more than one genome location.  We may wish
 * to keep only the highest-scoring alignment(s) for each tag.
 * [see LastToCompactMode]
 *
 * <h3> Discarding ambiguous tags (recommended in tag-seeds.txt) </h3>
 *
 * There may be some tags that map to more than one location (with equal
 * scores).  We may wish to discard such multi-mapping tags.
 * [see LastToCompactMode]
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
 * -z: write counts using a more compact single line format
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
public class LastagAligner extends AbstractAligner {
    private static final Log LOG = LogFactory.getLog(LastagAligner.class);

    enum MatchQuality {
        PERFECT_MATCHES_ONLY,
        ONE_MISMATCH,
        TWO_MISMATCHES,
        IDENTITY_95,
        BEST_MATCH
    }

    private static final String LASTDB_EXEC = SystemUtils.IS_OS_WINDOWS ? "lastdb.exe" : "lastdb";
    private static final String LASTAG_EXEC = SystemUtils.IS_OS_WINDOWS ? "lastag.exe" : "lastag";

    // parameters that control Lastdb and/or Lastag indirectly
    private static final String LASTDB_DEFAULT_MEMORY     = "2G"; // -s:
    private static final int DEFAULT_MAX_GAPS_ALLOWED     = 0;
    private static final int MIN_READ_LENGTH_APPROX       = 20;
    private static final String DEFAULT_MATCH_QUALITY     = "IDENTITY_95";

    // fixed values for Lastag command line - setAlignerOptions does not change these
    private static final int DEFAULT_STRAND_DIRECTION     = 2; // -s: strand: 0 = reverse, 1 = forward, 2 = both (2 for DNA, 1 for protein)
    private static final int DEFAULT_OUTPUT_TYPE          = 3; // -j: output type: 0 = match counts, 1 = gapless, 2 = redundant gapped, 3 = gapped, 4 = probabilities, 5 = centroid
    private static final int DEFAULT_OUTPUT_FORMAT        = 1; // -f: output format: 0 = tabular, 1 = maf
    private static final String DEFAULT_QUERY_BATCH_SIZE  = "1G"; // -i: query batch size

    // suggested defaults for Lastag command line - setAlignerOptions may change these
    private static final int DEFAULT_GAP_OPENING_COST     = 2; // -a
    private static final int DEFAULT_GAP_EXTENSION_COST   = 1; // -b
    private static final int DEFAULT_MIN_SEED_DEPTH       = 15; // -l
    private static final int DEFAULT_MISMATCH_COST        = 1;
    private static final int DEFAULT_MATCH_TO_MISMATCH_RATIO_NT = 1;
    private static final int DEFAULT_MATCH_TO_MISMATCH_RATIO_CS = 2;

    // parameters to resolve Lastag command-line options
    private int maxGapsAllowed = DEFAULT_MAX_GAPS_ALLOWED;
    private int gapOpeningCost = DEFAULT_GAP_OPENING_COST;
    private int gapExtentionCost = DEFAULT_GAP_EXTENSION_COST;
    private int minSeedDepth = DEFAULT_MIN_SEED_DEPTH;
    private int matchScore = -1; // -1 indicates that value has not been set by user
    private int misMatchCost = -1; // -1 indicates that value has not been set by user

    private MatchQuality matchQuality = MatchQuality.valueOf(DEFAULT_MATCH_QUALITY);

    // Can not yet initialize matrixFilename, since relevant members not set (e.g. work directory)
    private String matrixFilename = "";


    public LastagAligner() {
        super();
        extensions = new String[]{"prj"};
    }

    public void setProperties(final GroovyProperties properties) {
        pathToExecutables = properties.get(GobyPropertyKeys.EXECUTABLE_PATH_LASTAG, ".");
    }

    public String getDefaultDbNameForReferenceFile(final File referenceFile) {
        return FilenameUtils.getBaseName(referenceFile.getName()).replaceAll(" ", "_");
    }

    /**
     * Return the reference file converter used by aligner algorithm
     */
    public CompactToFastaMode getReferenceCompactToFastaConverter() {
        CompactToFastaMode processor = new CompactToFastaMode();
        processor.setHashOutputFilename(true);
        processor.setIndexToHeader(false);
        // Reference conversion (always from nt-space) *MAY* be needed by alignment algorithm to match colorspace reads platforms
        processor.setOutputColorMode(colorSpace); // lastag requires this conversion when dealing with colorSpace
        processor.setOutputFakeNtMode(false); // lastag uses real colorSpace
        processor.setOutputFakeQualityMode(false); // lastag currently does not use (fake) qualities
        // Filter using the default alphabet for specified output mode
        processor.setAlphabet((colorSpace) ? "0123ACTG" : "ACTG"); // lastag expects colorspace with nt-adaptors, or just nt
        return processor;
    }


    /**
     * Return the reads file converter used by aligner algorithm
     */
    public CompactToFastaMode getReadsCompactToFastaConverter() {
        CompactToFastaMode processor = new CompactToFastaMode();
        processor.setIndexToHeader(true);
        processor.setHashOutputFilename(true);
        // Since colorSpace is determined by reads platform, colorspace conversion is never needed for reads
        processor.setOutputFakeNtMode(false); // lastag uses real colorSpace
        processor.setOutputFakeQualityMode(false); // lastag currently does not use (fake) qualities
        // Filter using the default alphabet for specified output mode
        processor.setAlphabet((colorSpace) ? "0123ACTG" : "ACTG"); // lastag expects colorspace with nt-adaptors, or just nt
        processor.setTrimAdaptorLength((colorSpace) ? 2 : 0); // prevent aligner from finding false matches with adaptor at the start of the sequence
        return processor;
    }

    /**
     * Returns LastToCompact processor, initialized with correct input file
      */
    // FROM LastagAligner.align()
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
    public void setAlignerOptions(String options) {
        maxGapsAllowed = DEFAULT_MAX_GAPS_ALLOWED;
        gapOpeningCost = DEFAULT_GAP_OPENING_COST;
        matchQuality   = MatchQuality.valueOf(DEFAULT_MATCH_QUALITY);
        alignerOptions = "";
        if (options == null) {
            return;
        }
        final String[] opts = options.split("[,=]");

        // options specific to LastagAligner
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
            assert (matchScore!=-999) : "Parsing error";
        }
        if (CLI.isKeywordGiven(opts, "q")) {
            misMatchCost = CLI.getIntOption(opts, "q", -999);
            assert (misMatchCost!=-999) : "Parsing error";
        }

        // options specific to "bwa aln" tool
        // FLOAT -t -g
        alignerOptions += (CLI.isKeywordGiven(opts, "t")) ? " -t " + CLI.getFloatOption(opts, "t", -999) : "";
        alignerOptions += (CLI.isKeywordGiven(opts, "g")) ? " -g " + CLI.getFloatOption(opts, "g", -999) : "";
        // INT -c -x -y -a -m -l -k -w
        alignerOptions += (CLI.isKeywordGiven(opts, "c")) ? " -c " + CLI.getIntOption  (opts, "c", -999) : "";
        alignerOptions += (CLI.isKeywordGiven(opts, "x")) ? " -x " + CLI.getIntOption  (opts, "x", -999) : "";
        alignerOptions += (CLI.isKeywordGiven(opts, "y")) ? " -y " + CLI.getIntOption  (opts, "y", -999) : "";
        alignerOptions += (CLI.isKeywordGiven(opts, "m")) ? " -m " + CLI.getIntOption  (opts, "m", -999) : ""; // getSeedMaxMultiplicity()
        alignerOptions += (CLI.isKeywordGiven(opts, "k")) ? " -k " + CLI.getIntOption  (opts, "k", -999) : "";
        alignerOptions += (CLI.isKeywordGiven(opts, "w")) ? " -w " + CLI.getIntOption  (opts, "w", -999) : "";
        // BOOLEAN -v
        alignerOptions += (CLI.isKeywordGiven(opts, "v")) ? " -v "                                       : "";
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
    private void validateScoreOptions() {
        if (misMatchCost < 0) {
            misMatchCost = DEFAULT_MISMATCH_COST;
        }
        if (matchScore < 0) {
            final int ratio = (colorSpace) ? DEFAULT_MATCH_TO_MISMATCH_RATIO_CS : DEFAULT_MATCH_TO_MISMATCH_RATIO_NT;
            matchScore = ratio * misMatchCost;
        }
        if (mParameter<1) {
            mParameter = 1; // must be at least 1
        }
    }

    private String scoreOptions() {
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
        } catch (IOException e) {
            throw new RuntimeException(e);
        } finally {
            IOUtils.closeQuietly(writer);
        }
    }

    /**
     * Initialize lastag's seedMaxMultiplicty in terms of our mParameter. Lastag will only extend seeds that are
     * unambigous, according to this value.
     * <p/>
     * Typically: seedMaxMultiplicity >= mParameter
     * <p/>
     * Rational: For each read, lastag extracts contiguous subsequences (seeds) that match at most this number of times
     * to the reference genome.  The initial seed alignments are used to evaluate full read-length alignments at each
     * location.  After scoring the full-length alignments, the best ones are returned.  The number returned could be
     * less than seedMaxMultiplicity.
     * <p/>
     * Since this algorithm uses heuristics to seed and filter alignments, it is difficult to characterize the results.
     * <p/>
     * Our strategy is to tune the lastag algorithm so that it will generate lots of results, and then by filtering
     * them, we can *try to* guarantee certain properties of the output.
     * <p/>
     * To maximize the number of alignments that match the genome at most mParameter times, we select seedMaxMultiplicity
     * greater than mParameter.
     */
    private int getSeedMaxMultiplicity() {
        validateScoreOptions();
        return 2*mParameter;
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
        int numMatch;
        int numMisMatch;
        switch (matchQuality) {
            case PERFECT_MATCHES_ONLY:
                numMatch = minReadLength;
                numMisMatch = 0;
                break;
            case ONE_MISMATCH:
                numMatch = minReadLength-1;
                numMisMatch = 1;
                break;
            case TWO_MISMATCHES:
                numMatch = minReadLength-2;
                numMisMatch = 2;
                break;
            case IDENTITY_95:
                numMisMatch = (int) Math.ceil(minReadLength * 0.05);
                numMatch = minReadLength - numMisMatch;
                break;
            case BEST_MATCH:
                numMatch = (int) (minReadLength / 3);
                numMisMatch = 0;
                // ensures threshold is low enough to search all potential candidates
                break;
            default:
                numMisMatch = (int) Math.floor(minReadLength * 0.10);
                numMatch = minReadLength - numMisMatch;
        }
        final int eThreshAssumingSNPs = numMatch*matchScore - numMisMatch*misMatchCost;
        final int eThreshAssumingGAPs = numMatch*matchScore - numMisMatch*gapExtentionCost - maxGapsAllowed*gapOpeningCost;
        return Math.max(0, Math.min(eThreshAssumingSNPs, eThreshAssumingGAPs));
    }

    /**
     *  value used by last v52. 3/5 of the gapped score should be achieved by the gapless extension.
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
     *                                  the database basename (prefix)
     * @return
     * @throws IOException
     * @throws InterruptedException
     */
    public File[] indexReference(final File referenceFileOrDbBasename)
            throws IOException, InterruptedException {

        if (isDatabaseBasename(referenceFileOrDbBasename.toString())) {
            return matchingExtension(referenceFileOrDbBasename.toString(), extensions);
        }
        final String databasePathPrefix = getDatabasePath(databaseDirectory, databaseName);
        if (isDatabaseBasename(databasePathPrefix)) {
            return matchingExtension(databasePathPrefix, extensions);
        }
        final File fastaReferenceFile = prepareReference(referenceFileOrDbBasename);
        forceMakeParentDir(databasePathPrefix);
        final ExecuteProgram executor = new ExecuteProgram();
        final String command = String.format("%s/%s -v -s%s -a %s %s %s ",
                pathToExecutables,
                LASTDB_EXEC,
                LASTDB_DEFAULT_MEMORY, // -s
                getAlphabet(),
                databasePathPrefix,
                fastaReferenceFile.toString().replaceAll(" ", "\\ "));
        executor.executeToLog(command, LastagAligner.class, Level.INFO, "LastAG[__OUTPUT_TAG__]: ");
        listFiles("Database files", databaseDirectory);
        return matchingExtension(databasePathPrefix, extensions);
    }

    /**
     * Align.
     *
     * @param referenceFile  Compact or native database basename for reference sequences.
     * @param readsFile      Compact or native format (i.e., fasta, fastq) aligner read format.
     * @param outputBasename Basename where to write the compact alignment.
     * @return
     * @throws InterruptedException
     * @throws IOException
     */
    public File[] align(final File referenceFile, final File readsFile, final String outputBasename) throws InterruptedException, IOException {
        assert pathToExecutables != null : "path to executables must be defined.";
        validateScoreOptions();

        try {
            final File nativeReads = prepareReads(readsFile);
            if (databaseName == null) {
                databaseName = getDefaultDbNameForReferenceFile(referenceFile);
            }
            final File[] indexedReference = indexReference(referenceFile);
            LOG.info("Searching..");
            forceMakeParentDir(outputBasename);
            // if (colorSpace) { prepareMatrixFile(); }

            final ExecuteProgram executor = new ExecuteProgram();

            final String command = String.format("%s/%s -a%d -b%d -d%d -e%d " +
                    " -m%d -l%d -s%d -j%d -z -i%s -v -f%d %s %s -o %s %s %s",
                    pathToExecutables,
                    LASTAG_EXEC,
                    gapOpeningCost,                      // -a
                    gapExtentionCost,                    // -b
                    dThresholdOption(minReadLength),     // -d
                    eThresholdOption(minReadLength),     // -e
                    getSeedMaxMultiplicity(),            // -m
                    minSeedDepth,                 // -l
                    DEFAULT_STRAND_DIRECTION,     // -s
                    DEFAULT_OUTPUT_TYPE,          // -j
                    DEFAULT_QUERY_BATCH_SIZE,     // -i
                    DEFAULT_OUTPUT_FORMAT,        // -f
                    alignerOptions,
                    scoreOptions(),                      // previously substitutionMatrixOption()
                    outputBasename,                      // -o
                    getDatabasePath(databaseDirectory, databaseName),
                    nativeReads);

            executor.executeToLog(command, LastagAligner.class, Level.INFO,
                    "LastAG[__OUTPUT_TAG__]: ");

            // convert native alignment into compact reads
            final File[] buildResults = processAlignment(referenceFile, readsFile, outputBasename);
            FileUtils.deleteQuietly(new File(outputBasename + ".maf"));
            FileUtils.deleteQuietly(new File(outputBasename + ".counts"));

            // also delete the colorSpaceMatrix file, if created
            FileUtils.deleteQuietly(new File(matrixFilename));

            return buildResults;

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
