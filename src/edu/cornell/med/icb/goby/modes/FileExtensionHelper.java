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

package edu.cornell.med.icb.goby.modes;

import java.io.File;

/**
 * Some basic stats on MAQ Map, FASTA, and FASTQ files.
 *
 * @author Kevin Dorff
 */
public final class FileExtensionHelper {
    static final String[] FASTX_FILE_EXTS = {
            ".fa.gz", ".fna.gz", ".fasta.gz", ".fq.gz", ".fnq.gz", ".fastq.gz",
            ".csfasta", ".csfasta.gz", ".csfastq", ".csfastq.gz",
            ".csfa", ".csfa.gz", ".csfq", ".csfq.gz",
            ".fa", ".fna", ".fasta", ".fq", ".fnq", ".fastq",
            ".txt", ".txt.gz"};

    static final String[] COMPACT_READS_FILE_EXTS = {".compact-reads"};

    static final String[] COMPACT_ALIGNMENT_FILE_EXTS = {
            ".entries", ".header", ".tmh", ".stats"
    };

    /**
     * Possible types for compact files.
     */
    public enum CompactFileType {
        alignment,
        reads,
        unknown
    }

    /**
     * Determine the type of data within a file in compact reads format.
     * @param file The file to check
     * @return The type of data this file is likely to contain
     */
    public static CompactFileType determineCompactFileType(final File file) {
        return determineCompactFileType(file.toString());
    }

    /**
     * Determine the type of data within a file in compact reads format.
     * @param file The name of file to check
     * @return The type of data this file is likely to contain
     */
    public static CompactFileType determineCompactFileType(final String file) {
        // does the file represent an alignment?
        for (final String ext : COMPACT_ALIGNMENT_FILE_EXTS) {
            if (file.endsWith(ext)) {
                return CompactFileType.alignment;
            }
        }

        // does the file represent a series of reads?
        for (final String ext : COMPACT_READS_FILE_EXTS) {
            if (file.endsWith(ext)) {
                return CompactFileType.reads;
            }
        }

        return CompactFileType.unknown;
    }

    private FileExtensionHelper() {
        super();
    }
}
