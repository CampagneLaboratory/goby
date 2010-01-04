package edu.cornell.med.icb.goby.aligners;

import edu.cornell.med.icb.goby.util.GroovyProperties;

import java.io.File;
import java.io.IOException;

/**
 * Abstracts the details of calling an aligner.
 *
 * @author Fabien Campagne
 *         Date: Jul 9, 2009
 *         Time: 11:38:59 AM
 */
public interface Aligner {
    /**
     * Convert the compact reads file into the native read input format of the aligner. If the reads are already provided
     * in native format, no conversion is performed (e.g., Fasta for an aligner that reads Fasta).
     *
     * @param compactReadsFile Reads.
     * @return The file that contains the reads converted to native aligner format (e.g., fasta, fastq or other)
     */
    File prepareReads(File compactReadsFile) throws IOException;

    /**
     * Index the reference compact file and returns the files that make up the index.
     *
     * @param referenceFileOrDbBasename The compact file with the reference sequence to index OR
     * the database basename (prefix)
     * @return
     */
    File[] indexReference(File referenceFileOrDbBasename) throws IOException, InterruptedException;

    /**
     * Convert the native alignment format into compact alignment format.
     *
     * @param outputBasename Basename where to write the compact alignment.
     * @return The set of files that make up the result alignment.
     * @throws InterruptedException Thrown when the task is interrupted for any reason.
     */
    File[] processAlignment(File referenceFile, File readsFile, String outputBasename) throws InterruptedException, IOException;

    /**
     * Align reads to a reference. Convert reads and reference to aligner native format if necessary,
     * index the reference if needed, perform the alignment, convert the resuls to compact alignment format.
     *
     * @param referenceFile  Compact or native database basename for reference sequences.
     * @param readsFile      Compact or native format (i.e., fasta, fastq) aligner read format.
     * @param outputBasename Basename where to write the compact alignment.
     * @return The set of files that make up the result alignment.
     * @throws InterruptedException Thrown when the task is interrupted for any reason.
     */
    File[] align(File referenceFile, File readsFile, String outputBasename) throws InterruptedException, IOException;

    /**
     * Indicate whether the reads are in color-space
     *
     * @param colorSpace True for color-space reads, false otherwise.
     */
    void setColorSpace(boolean colorSpace);

    /**
     * Configure the path where databases will be written.
     *
     * @param path
     */
    void setDatabaseDirectory(String path);

    /**
     * Configure the path to aligner executables.
     *
     * @param path
     */
    void setPathToExecutables(String path);

    /**
     * Configure the work directory.
     *
     * @param path
     */
    void setWorkDirectory(String path);


    /**
     * Set parameters to filter reads for quality when importing from native format into compact format. See LastToCompactMode for
     * details.
     *
     * @param qualityFilterParameters i.e., threshold=0.05 to keep reads at most 5% bases different from reference.
     */
    void setQualityFilterParameters(String qualityFilterParameters);

    /**
     * Set parameters to filter reads for quality when importing from native format into compact format. See LastToCompactMode for
     * details.
     *
     * @param mParameter i.e. mParameter=2 to keep reads with at most 2 optimal-scoring matches across target
     */
    void setAmbiguityThreshold(int mParameter);



    /**
     * Instruct the aligner to filter reads with the read set encoded in this file.
     *
     * @param readIndexFilterFile
     */
    void setReadIndexFilter(File readIndexFilterFile);

    /**
     * Set aligner specific options. Options are comma separated, with the syntax key1=value1,key2=value2.
     * Keys are aligner dependent.
     *
     * @param options
     */
    void setAlignerOptions(String options);

    /**
     * Instruct the aligner to filter reference sequences with the read set encoded in this file.
     *
     * @param referenceIndexFilterFile
     */
    void setReferenceIndexFilter(File referenceIndexFilterFile);

    /**
     * Specify the properties to use.
     * @param properties
     */
    void setProperties(final GroovyProperties properties);

    void setDatabaseName(final String databaseName);

    String getDefaultDbNameForReferenceFile(final File referenceFile);

    String getDatabaseName();

    boolean isDatabaseBasename(final String databasePathPrefix);

    String getAlphabet();
}
