package edu.cornell.med.icb.goby.aligners;

import edu.cornell.med.icb.goby.modes.AbstractAlignmentToCompactMode;
import edu.cornell.med.icb.goby.modes.CompactToFastaMode;
import edu.cornell.med.icb.goby.modes.SAMToCompactMode;
import edu.cornell.med.icb.util.ExecuteProgram;
import edu.cornell.med.icb.util.GobyPropertyKeys;
import edu.cornell.med.icb.util.GroovyProperties;
import edu.mssm.crover.cli.CLI;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang.RandomStringUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.log4j.Level;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;


/**
 * Facade to the BWA aligner. For details about BWA, see http://maq.sourceforge.net/bwa-man.shtml
 *
 * <h3> Description </h3>
 * BWA  is  a  fast  light-weighted  tool  that  aligns  relatively  short
 * sequences  (queries)  to a sequence database (targe), such as the human
 * reference genome. It implements two different algorithms, both based on
 * Burrows-Wheeler  Transform  (BWT).  The first algorithm is designed for
 * short queries up to ~200bp with low error rate (<3%).  It  does  gapped
 * global alignment w.r.t. queries, supports paired-end reads, and is one of
 * the fastest short read alignment algorithms to date while also  vis- iting
 * suboptimal hits. The second algorithm, dBWT-SW, is designed for long reads
 * with more errors. It performs heuristic  Smith-Waterman-like alignment  to
 * find high-scoring local hits (and thus chimera). On low- error short
 * queries, dBWT-SW is slower and less accurate than the first algorithm, but
 * on long queries, it is better.
 * <p/>
 * For  both  algorithms,  the  database  file in the FASTA format must be
 * first indexed with the `index' command, which  typically  takes  a  few
 * hours.  The first algorithm is implemented via the `aln' command, which
 * finds the suffix array (SA) coordinates of good hits of each individual
 * read,  and  the `samse/sampe' command, which converts SA coordinates to
 * chromosomal coordinate and pairs reads (for `sampe'). The second  algo-
 * rithm is invoked by the `dbtwsw' command. It works for single-end reads
 * only.
 *  <p/>
 * <h3> Alignment Accuracy </h3>
 * When  seeding is disabled, BWA guarantees to find an alignment contain- ing
 * maximum maxDiff differences including maxGapO gap  opens  which  do not
 * occur  within nIndelEnd bp towards either end of the query. Longer gaps may
 * be found if maxGapE is positive, but it is not  guaranteed  to find  all
 * hits. When seeding is enabled, BWA further requires that the first seedLen
 * subsequence contains no  more  than  maxSeedDiff  differ- ences.
 *  <p/>
 * When gapped alignment is disabled, BWA is expected to generate the same
 * alignment as Eland, the Illumina alignment  program.  However,  as  BWA
 * change  `N'  in  the  database  sequence to random nucleotides, hits to
 * these random sequences will also be counted. As a consequence, BWA  may mark
 * a  unique  hit  as a repeat, if the random sequences happen to be identical
 * to the sequences which should be unqiue in the database. This random
 * behaviour will be avoided in future releases.
 *  <p/>
 * By default, if the best hit is no so repetitive (controlled by -R), BWA also
 * finds all hits contains one more mismatch;  otherwise,  BWA  finds all
 * equally best hits only. Base quality is NOT considered in evaluat- ing hits.
 * In pairing, BWA searches, among the found hits under the con- straint  of
 * the  maxOcc  option,  for pairs within maxInsSize and with proper
 * orientation.
 *
 * <h3> Memory Requirement </h3>
 * With bwtsw algorithm, 2.5GB memory is required for  indexing  the  com-
 * plete  human  genome sequences. For short reads, the `aln' command uses
 * ~2.3GB memory and the `sampe' command uses ~3.5GB.
 *
 * <h3> Speed </h3>
 * Indexing the human genome sequences takes 3 hours with bwtsw algorithm.
 * Indexing  smaller  genomes  with IS or divsufsort algorithms is several
 * times faster, but requires more memory.
 * <p/>
 * Speed of alignment is largely determined by the error rate of the query
 * sequences (r). Firstly, BWA runs much faster for near perfect hits than for
 * hits with many differences, and it stops searching for a  hit  with l+2
 * differences if a l-difference hit is found. This means BWA will be very slow
 * if r is high because in this case BWA has to visit hits  with many
 * differences and looking for these hits is expensive. Secondly, the alignment
 * algorithm behind makes the speed sensitive to  [k  log(N)/m], where  k is
 * the maximum allowed differences, N the size of database and m the length of
 * a query. In practice, we choose k w.r.t. r  and  there- fore  r is the
 * leading factor. I would not recommend to use BWA on data with r>0.02.
 * <p/>
 * Pairing is slower for shorter reads. This  is  mainly  because  shorter
 * reads  have more spurious hits and converting SA coordinates to chromo- somal
 * coordinates are very costly.
 *  <p/>
 * In a practical experiment, BWA is able to map 2 million 32bp reads to a
 * bacterial  genome  in  several minutes, map the same amount of reads to
 * human X chromosome in 8-15 minutes and to the  human  genome  in  15-25
 * minutes.  This  result  implies that the speed of BWA is insensitive to the
 * size of database and therefore BWA is more efficient when the data- base  is
 * sufficiently large. On smaller genomes, hash based algorithms are usually
 * much faster.
 *
 * <pre>
 * index  bwa index [-p prefix] [-a algoType] [-c] <in.db.fasta>
 *
 *        Index database sequences in the FASTA format.
 *
 *        OPTIONS:
 *
 *        -c        Build color-space index. The input fast should  be  in
 *                  nucleotide space.
 *
 *        -p STR    Prefix of the output database [same as db filename]
 *
 *        -a STR    Algorithm   for   constructing  BWT  index.  Available
 *                  options are:
 *
 *                  is     IS linear-time algorithm for constructing  suf-
 *                         fix  array. It requires 5.37N memory where N is
 *                         the size of  the  database.  IS  is  moderately
 *                         fast,  but  does  not work with database larger
 *                         than 2GB. IS is the default  algorithm  due  to
 *                         its  simplicity. The current codes for IS algo-
 *                         rithm are reimplemented by Yuta Mori.
 *
 *                  bwtsw  Algorithm implemented in  BWT-SW.  This  method
 *                         works  with the whole human genome, but it does
 *                         not work with database smaller than 10MB and it
 *                         is usually slower than IS.
 *
 *
 *
 * aln    bwa aln [-n maxDiff] [-o maxGapO] [-e maxGapE] [-d nDelTail] [-i
 *        nIndelEnd] [-k maxSeedDiff] [-l seedLen] [-t nThrds] [-cRN]  [-M
 *        misMsc]  [-O  gapOsc]  [-E gapEsc] <in.db.fasta> <in.query.fq> >
 *        <out.sai>
 *
 *        Find the SA coordinates of the input reads.
 *
 *        Maximum  maxSeedDiff differences  are  allowed  in  the first seedLen
 *        subsequence and maximum maxDiff differences are allowed in the whole
 *        sequence.
 *
 *        OPTIONS:
 *
 *        -n NUM    Maximum edit distance if the  value  is  INT,  or  the
 *                  fraction  of  missing alignments given 2% uniform base
 *                  error rate if FLOAT. In the latter case,  the  maximum
 *                  edit  distance  is  automatically chosen for different
 *                  read lengths. [0.04]
 *
 *                  int bwa_cal_maxdiff(int length, double nt_err_rate, double fn_rate)
 *                  {
 *                      double elambda = exp(-length * nt_err_rate);
 *                      double sum, y = 1.0;
 *                      int k, x = 1;
 *                      for (k = 1, sum = elambda; k < 1000; ++k) {
 *                          y *= length * nt_err_rate;
 *                          x *= k;
 *                          sum += elambda * y / x;
 *                          if (1.0 - sum < fn_rate) return k;
 *                      }
 *                      return 2;
 *                  }
 *
 *
 *        -o INT    Maximum number of gap opens [1]
 *
 *        -e INT    Maximum number of gap extensions, -1 for  k-difference
 *                  mode (disallowing long gaps) [-1]
 *
 *        -d INT    Disallow  a  long  deletion  within INT bp towards the
 *                  3'-end [16]
 *
 *        -i INT    Disallow an indel within INT bp towards the ends [5]
 *
 *        -l INT    Take the first INT subsequence  as  seed.  If  INT  is
 *                  larger  than  the query sequence, seeding will be dis-
 *                  abled. For long reads, this option is typically ranged
 *                  from 25 to 35 for `-k 2'. [inf]
 *
 *        -k INT    Maximum edit distance in the seed [2]
 *
 *        -t INT    Number of threads (multi-threading mode) [1]
 *
 *        -M INT    Mismatch  penalty.  BWA will not search for suboptimal
 *                  hits with a score lower than (bestScore-misMsc). [3]
 *
 *        -O INT    Gap open penalty [11]
 *
 *        -E INT    Gap extension penalty [4]
 *
 *        -R INT    Proceed with suboptimal alignments  if  there  are  no
 *                  more  than  INT  equally  best  hits. This option only
 *                  affects paired-end mapping. Increasing this  threshold
 *                  helps  to  improve the pairing accuracy at the cost of
 *                  speed, especially for short reads (~32bp).
 *
 *        -c        Reverse query but not complement it, which is required
 *                  for alignment in the color space.
 *
 *        -N        Disable  iterative  search. All hits with no more than
 *                  maxDiff differences will be found. This mode  is  much
 *                  slower than the default.
 * </pre>
 *
 *
 * @author Fabien Campagne
 *         Date: Jul 9, 2009
 *         Time: 5:21:08 PM
 */
public class BWAAligner extends AbstractAligner {
    private static final Log LOG = LogFactory.getLog(BWAAligner.class);

    private static final String BWA_INDEX_EXEC = "bwa index";
    private static final String BWA_ALIGN_EXEC = "bwa aln";
    private static final String BWA_SAMSE_EXEC = "bwa samse";
    private static final int DEFAULT_SEED_LENGTH = 35;

    private String databasePrefix;
    private String saiBinaryFilename = "";
    private String samBinaryFilename = "";

    private int seedLength = DEFAULT_SEED_LENGTH;

    public BWAAligner() {
        super();
        extensions = new String[]{"rsa", "rpac", "rbwt", "pac", "bwt", "ann", "amb", "sa"};
    }

    public void setProperties(final GroovyProperties properties) {
        pathToExecutables = properties.get(GobyPropertyKeys.EXECUTABLE_PATH_BWA, ".");
    }

    public String getDefaultDbNameForReferenceFile(final File referenceFile) {
        return FilenameUtils.getBaseName(referenceFile.getName()) + ".fasta";
    }

    /**
     * Return the reference file converter used by aligner algorithm
     */
    public CompactToFastaMode getReferenceCompactToFastaConverter() {
        CompactToFastaMode processor = new CompactToFastaMode();
        processor.setHashOutputFilename(true);
        processor.setIndexToHeader(false);
        // Reference conversion (always from nt-space) *MAY* be needed by alignment algorithm to match colorspace reads platforms
        processor.setOutputColorMode(false); // BWA handles the conversion
        processor.setOutputFakeNtMode(false); // BWA handles the conversion
        processor.setOutputFakeQualityMode(false); // BWA handles the conversion
        // Filter using the default alphabet for specified output mode
        processor.setAlphabet("ACTG"); // BWA handles the conversion
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
        processor.setOutputFakeNtMode(colorSpace); // BWA uses fakeNt to represent colorSpace data
        processor.setOutputFakeQualityMode(colorSpace); // BWA uses fakeNt to represent colorSpace data
        // Filter using the default alphabet for specified output mode
        processor.setAlphabet("ACTG"); // BWA expects fake-nt or nt, either way, "ACTG" is the alphabet
        processor.setTrimAdaptorLength((colorSpace) ? 2 : 0); // BWA's solid2fastq remove 2 bases - prevent aligner from finding false matches with adaptor at the start of the sequence
        return processor;
    }


    /**
     * Returns SAMToCompact processor, initialized with correct input file
      */
    public AbstractAlignmentToCompactMode getNativeAlignmentToCompactMode(final String outputBasename) {
        // can not initialize, unless correct input file exists
        assert (samBinaryFilename!=null) : "Can not initialize, unless SAM Binary input file exists";
        assert (new File(samBinaryFilename).exists()) : "Can not initialize, unless correct SAM Binary file exists";
        final SAMToCompactMode processor = new SAMToCompactMode();
        processor.setInputFile(samBinaryFilename);
        processor.setNumberOfReads(numberOfReads);
        return processor;
    }

    /**
     * Parse and validate aligner specific options.
     * Input options are comma separated, with the syntax key1=value1,key2=value2.
     * "alignerOptions" format should match aligner's command-line specification e.g. "-key1 value1 -key2 value2"
     * This method is declared abstract so that aligners have control over the parsing of their arguments.
     */
    public void setAlignerOptions(String options) {
        seedLength = DEFAULT_SEED_LENGTH;
        alignerOptions = "";
        if (options == null) {
            return;
        }
        final String[] opts = options.split("[,=]");

        // seedLength parameter
        if (CLI.isKeywordGiven(opts, "l")) {
            seedLength = CLI.getIntOption(opts, "l", DEFAULT_SEED_LENGTH);
        }

        alignerOptions += (CLI.isKeywordGiven(opts, "n")) ? " -n " + CLI.getIntOption(opts, "n", -999) : "";
        alignerOptions += (CLI.isKeywordGiven(opts, "o")) ? " -o " + CLI.getIntOption(opts, "o", -999) : "";
        alignerOptions += (CLI.isKeywordGiven(opts, "e")) ? " -e " + CLI.getIntOption(opts, "e", -999) : "";
        alignerOptions += (CLI.isKeywordGiven(opts, "d")) ? " -d " + CLI.getIntOption(opts, "d", -999) : "";
        alignerOptions += (CLI.isKeywordGiven(opts, "i")) ? " -i " + CLI.getIntOption(opts, "i", -999) : "";
        alignerOptions += (CLI.isKeywordGiven(opts, "k")) ? " -k " + CLI.getIntOption(opts, "k", -999) : "";
        //alignerOptions += (CLI.isKeywordGiven(opts, "l")) ? " -l " + CLI.getIntOption(opts, "l", -999) : ""; // seedLength
        alignerOptions += (CLI.isKeywordGiven(opts, "t")) ? " -t " + CLI.getIntOption(opts, "t", -999) : "";
        alignerOptions += (CLI.isKeywordGiven(opts, "M")) ? " -M " + CLI.getIntOption(opts, "M", -999) : "";
        alignerOptions += (CLI.isKeywordGiven(opts, "O")) ? " -O " + CLI.getIntOption(opts, "O", -999) : "";
        alignerOptions += (CLI.isKeywordGiven(opts, "E")) ? " -E " + CLI.getIntOption(opts, "E", -999) : "";
        alignerOptions += (CLI.isKeywordGiven(opts, "c")) ? " -c "                                     : "";
        alignerOptions += (CLI.isKeywordGiven(opts, "N")) ? " -N "                                     : "";
        //alignerOptions += (CLI.isKeywordGiven(opts, "R")) ? " -R " : ""; // R param is only applicable for paired-end read alignment
        assert (!alignerOptions.contains("-999")) : "Parsing error.";
    }


    public File[] indexReference(final File referenceFileOrDbBasename)
            throws IOException, InterruptedException {
        databasePrefix = referenceFileOrDbBasename.toString();
        if (isDatabaseBasename(databasePrefix)) {
            return matchingExtension(referenceFileOrDbBasename.toString(), extensions);
        }
        databasePrefix = getDatabasePath(databaseDirectory, databaseName);
        if (isDatabaseBasename(databasePrefix)) {
            return matchingExtension(databasePrefix, extensions);
        }
        final File fastaReferenceFile = prepareReference(referenceFileOrDbBasename);

        forceMakeParentDir(databasePrefix);
        final ExecuteProgram executor = new ExecuteProgram();
        final String command = String.format("%s/%s  -a %s %s -p %s %s ",
                pathToExecutables, BWA_INDEX_EXEC,
                databaseIndexingStyle(fastaReferenceFile),
                colorSpace ? "-c" : "",
                databasePrefix,
                fastaReferenceFile.toString().replaceAll(" ", "\\ "));
        executor.executeToLog(command, BWAAligner.class, Level.INFO, "BWA[__OUTPUT_TAG__]: ");
        listFiles("Database files", databaseDirectory);
        return matchingExtension(databasePrefix, extensions);
    }

    private String databaseIndexingStyle(final File fastaReferenceFile) {
        // use bwtsw for genome fasta files larger than 5M, if otherwise.
        return fastaReferenceFile.length() < 5 * 1024 * 1024 ? "is" : "bwtsw";

    }

    /**
     * Return string option of form "-l <seedLength>".
     *
     * If a positive seedLength has been initialized, use this value. Otherwise,
     * use -l 35 to restrict the seed to the first 35 bp of each read. This
     * considerably speeds up searches with ~100 bp reads.
     *
     * Only use this option if the seedLength is shorter than the minimum read length,
     * otherwise bwa may fail (despite what the documentation says).
     */
    private String getLOption() {
        return (minReadLength>=seedLength) ? ("-l " + seedLength) : "";
    }


    public File[] align(final File referenceFile, final File readsFile, final String outputBasename) throws InterruptedException, IOException {
        assert pathToExecutables != null : "path to executables must be defined.";

        BufferedOutputStream saiOutputStream = null;
        BufferedOutputStream samOutputStream = null;
        try {
            final File nativeReads = prepareReads(readsFile);
            if (nativeReads.toString().contains(" ")) {
                throw new RuntimeException("BWA does not support spaces in filename: " + nativeReads);

            }
            if (databaseName == null) {
                databaseName = getDefaultDbNameForReferenceFile(referenceFile);
            }
            final File[] indexedReference = indexReference(referenceFile);
            saiBinaryFilename = FilenameUtils.concat(workDirectory, File.createTempFile(RandomStringUtils.randomAlphabetic(10), ".sai").getName());
            samBinaryFilename = FilenameUtils.concat(workDirectory, File.createTempFile(RandomStringUtils.randomAlphabetic(10), ".sam").getName());
            forceMakeParentDir(saiBinaryFilename);
            forceMakeParentDir(samBinaryFilename);
            LOG.info("Searching..");
            final ExecuteProgram executor = new ExecuteProgram();

            String command = String.format("%s/%s %s %s %s %s %s",
                    pathToExecutables, BWA_ALIGN_EXEC,
                    colorSpace ? "-c" : "",
                    getLOption(),
                    alignerOptions,
                    databasePrefix, nativeReads);
            saiOutputStream =  new BufferedOutputStream(new FileOutputStream(saiBinaryFilename));
            executor.executeStderrToLog(command, saiOutputStream, BWAAligner.class, Level.INFO,
                    "BWA[__OUTPUT_TAG__]: ");
            saiOutputStream.close();

            // convert sai to SAM format:
            command = String.format("%s/%s %s %s  %s",
                    pathToExecutables, BWA_SAMSE_EXEC,
                    databasePrefix, saiBinaryFilename, nativeReads);
            samOutputStream = new BufferedOutputStream(new FileOutputStream(samBinaryFilename));
            executor.executeStderrToLog(command, samOutputStream, BWAAligner.class, Level.INFO,
                    "BWA[__OUTPUT_TAG__]: ");
            samOutputStream.close();

            // convert native alignment into compact reads
            final File[] buildResults = processAlignment(referenceFile, readsFile, outputBasename);
            FileUtils.deleteQuietly(new File(saiBinaryFilename));
            FileUtils.deleteQuietly(new File(samBinaryFilename));

            return buildResults;

        } catch (IOException e) {
            throw new RuntimeException(e);
        } finally {
            IOUtils.closeQuietly(saiOutputStream);
            IOUtils.closeQuietly(samOutputStream);
        }
    }
}
