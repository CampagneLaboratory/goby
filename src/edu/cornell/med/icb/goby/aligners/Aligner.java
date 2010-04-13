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

import org.apache.commons.configuration.Configuration;

import java.io.File;
import java.io.IOException;

/**
 * Abstracts the details of calling an aligner.
 *
 * @author Fabien Campagne
 *         Date: Jul 9, 2009
 *         Time: 11:38:59 AM
 */
// TODO - replace setAmbiguityThreshold and setQualityFilterParameters with setFilterOptions()
public interface Aligner {
    /**
     * Convert the compact reads file into the native read input format of the aligner. If the
     * reads are already provided in native format, no conversion is performed (e.g., Fasta for
     * an aligner that reads Fasta).
     *
     * @param compactReadsFile Reads.
     * @return The file that contains the reads converted to native aligner format
     * (e.g., fasta, fastq or other)
     * @throws IOException
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
     * Specify the configuration to use for the aligner.
     * @param configuration The specific configuration/properties for the aligner
     */
    void setConfiguration(final Configuration configuration);

    void setDatabaseName(final String databaseName);

    String getDefaultDbNameForReferenceFile(final File referenceFile);

    String getDatabaseName();

    /**
     * Returns true if the reference is a database basename.
     * Searches for databasePathPrefix + "." + [extensions] to make sure all those
     * files exist.
     *
     * @param databasePathPrefix path prefix for database files
     * @return true if databasePathPrefix specifies an existing database
     */
    boolean isDatabaseBasename(final String databasePathPrefix);

    String getAlphabet();

    void setKeepTemporaryFiles(boolean keepTemporaryFiles);
}
