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

package edu.cornell.med.icb.goby.util;

import java.io.File;

/**
 * This class contains values and utilities related to file types used by Goby.
 *
 * @author Kevin Dorff
 */
public final class FileExtensionHelper {
    /**
     * Probable file extensions for FASTA/FASTQ files.
     */
    public static final String[] FASTX_FILE_EXTS = {
            ".fa.gz", ".fna.gz", ".fasta.gz", ".fq.gz", ".fnq.gz", ".fastq.gz",
            ".csfasta", ".csfasta.gz", ".csfastq", ".csfastq.gz",
            ".csfa", ".csfa.gz", ".csfq", ".csfq.gz",
            ".fa", ".fna", ".fasta", ".fq", ".fnq", ".fastq",
            ".txt", ".txt.gz", ".gz"
    };

    /**
     * File extenstion for sequence data in "compact reads" format.
     */
    public static final String[] COMPACT_READS_FILE_EXTS = {
            ".compact-reads"
    };

    /**
     * File extensions for alignment data in "compact reads" format.
     */
    public static final String[] COMPACT_ALIGNMENT_FILE_EXTS = {
            ".entries", ".header", ".tmh", ".stats", ".counts", ".index"
    };

   


    /**
     * Possible types for compact files.
     */
    public enum CompactFileType {
        /**
         * The file contains alignment data.
         */
        alignment,
        /**
         * The file contains reads or possibly reference data.
         */
        reads,
        /**
         * The file contents cannot be determined or is not a supported type.
         */
        unknown
    }

    /**
     * Determine the type of data within a file in compact reads format.
     *
     * @param file The file to check
     * @return The type of data this file is likely to contain
     */
    public static CompactFileType determineCompactFileType(final File file) {
        return determineCompactFileType(file.toString());
    }

    /**
     * Determine the type of data within a file in compact reads format.
     *
     * @param file The name of file to check
     * @return The type of data this file is likely to contain
     */
    public static CompactFileType determineCompactFileType(final String file) {
        if (hasCompactAlignmentExtension(file) || new File(file+".entries").canRead()) {
            // The file represents an alignment
            return CompactFileType.alignment;
        } else if (hasCompactReadsExtension(file)|| new File(file+".compact-reads").canRead()) {
            // The file represents a series of reads
            return CompactFileType.reads;
        } else {
            // the file extension does not match a known Goby compact type
            return CompactFileType.unknown;
        }
    }

    /**
     * Does the given file have an extension that matches a known Goby compact reads format?
     *
     * @param file The file to check
     * @return true if the file is likely to be a compact reads file
     */
    public static boolean hasCompactReadsExtension(final File file) {
        return hasCompactReadsExtension(file.toString());
    }

    /**
     * Does the given filename have an extension that matches a known Goby compact reads format?
     *
     * @param file The name of file to check
     * @return true if the filename is likely to be a compact reads file
     */
    public static boolean hasCompactReadsExtension(final String file) {
        for (final String ext : COMPACT_READS_FILE_EXTS) {
            if (file.endsWith(ext)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Does the given file have an extension that matches a known Goby compact alignment format?
     *
     * @param file The file to check
     * @return true if the file is likely to be a compact alignment file
     */
    public static boolean hasCompactAlignmentExtension(final File file) {
        return hasCompactAlignmentExtension(file.toString());
    }

    /**
     * Does the given filename have an extension that matches a known Goby compact alignment format?
     *
     * @param file The name of file to check
     * @return true if the filename is likely to be a compact alignment file
     */
    public static boolean hasCompactAlignmentExtension(final String file) {
        for (final String ext : COMPACT_ALIGNMENT_FILE_EXTS) {
            if (file.endsWith(ext)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Private constructor for utility class.
     */
    private FileExtensionHelper() {
        super();
    }
}
