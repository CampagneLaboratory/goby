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

package edu.cornell.med.icb.goby.alignments;

import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;

import java.io.IOException;

/**
 * Represents a slice of an alignment file. The slice is defined for this specific file by a start and
 * end offset, in bytes into the file.
 *
 * @author Fabien Campagne
 *         Date: 10/19/11
 *         Time: 1:08 PM
 */
public class FileSlice {
    public long startOffset;
    public long endOffset;
    public String basename;

    public FileSlice(long startOffset, long endOffset) {
        this.startOffset = startOffset;
        this.endOffset = endOffset;

    }

    public FileSlice(long startOffset, long endOffset, String basename) {
        this.startOffset = startOffset;
        this.endOffset = endOffset;
        this.basename = basename;
    }

    public static FileSlice COMPLETE_FILE(String basename) {
        return new FileSlice(0, Long.MAX_VALUE, basename);
    }

    public static FileSlice getSlice(AlignmentReaderFactory factory, String basename, GenomicRange range) throws IOException {
        if (range == null) {
            return FileSlice.COMPLETE_FILE(basename);
        }

        AlignmentReader reader = null;
        try {
            reader = factory.createReader(basename);
            reader.readHeader();

            DoubleIndexedIdentifier ids = new DoubleIndexedIdentifier(reader.getTargetIdentifiers());
            range.resolveChromosomeIndices(ids);
            reader.readIndex();

            long startOffset = reader.getStartByteOffset(range.startReferenceIndex, range.startPosition);
            long endOffset = reader.getEndByteOffset(range.endReferenceIndex, range.endPosition);
            if (startOffset==Long.MIN_VALUE) {
                startOffset=0;
            }
            return new FileSlice(startOffset, endOffset, basename);
        } finally {
            if (reader != null) {
                reader.close();
            }
        }

    }
}
