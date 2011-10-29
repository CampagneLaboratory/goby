/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 *  This file is part of the Goby IO API.
 *
 *     The Goby IO API is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     The Goby IO API is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with the Goby IO API.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.cornell.med.icb.goby.counts;

import edu.cornell.med.icb.goby.modes.CompactAlignmentToCountsMode;
import org.apache.commons.io.FileUtils;
import org.bdval.io.compound.CompoundDataOutput;
import org.bdval.io.compound.CompoundFileWriter;

import java.io.ByteArrayOutputStream;
import java.io.Closeable;
import java.io.File;
import java.io.IOException;

/**
 * Writes archives of counts information  for several sequences. CountsArchiveWriter leverages
 * the {@link org.bdval.io.compound.CompoundFileWriter}.
 *
 * @author Fabien Campagne
 *         Date: May 13, 2009
 *         Time: 7:06:57 PM
 */
public class CountsArchiveWriter implements Closeable {
    private final CompoundFileWriter compoundWriter;
    private CountsWriterI currentCountsWriterI;
    private String currentId;

    private ByteArrayOutputStream stream;
    private long totalBitsWritten;
    private int totalTransitions;
    private long totalBasesSeen;
    private long totalSitesSeen;
    private CountIndexBuilder indexBuilder = new CountIndexBuilder();
    private boolean verbose = false;

    /**
     * Initialize with a basename. Count information will be written to a file basename+".counts"
     *
     * @param basename Basename for the alignment/counts.
     * @throws IOException If an error occurs preparing the packaged counts.
     */
    public CountsArchiveWriter(final String basename) throws IOException {
        this(basename, CompactAlignmentToCountsMode.COUNT_ARCHIVE_MODIFIER_DEFAULT);
    }

    /**
     * Initialize with a basename. Count information will be written to a file basename+"."+countArchiveModifier
     *
     * @param basename             Basename for the alignment/counts.
     * @param countArchiveModifier the extension to write the archive.
     * @throws IOException If an error occurs preparing the packaged counts.
     */
    public CountsArchiveWriter(final String basename, final String countArchiveModifier) throws IOException {
        final String physicalFilename = basename + "." + countArchiveModifier;
        final File tobedeleted = new File(physicalFilename);
        // delete before writing since the compound file writer supports appending.
        FileUtils.deleteQuietly(tobedeleted);
        compoundWriter = new CompoundFileWriter(physicalFilename);

    }

    /**
     * Obtain a countsWriter to write counts for a sequence. Identifier provides a means to keep
     * track of multiple sequences.
     *
     * @param referenceIndex index of the sequence for which counts need to be recorded.
     * @param identifier     Identifier of the sequence for which counts need to be recorded.
     * @return A ready CountWriter implementation.
     * @throws IOException If an error occurs.
     */
    public CountsWriterI newCountWriter(final int referenceIndex, final String identifier) throws IOException {
        stream = new ByteArrayOutputStream(100000);
        currentCountsWriterI = new CountsWriter(stream);
        currentId = Integer.toString(referenceIndex) + "," + identifier;
        return currentCountsWriterI;
    }

    /**
     * Obtain a countsWriter to write counts for a sequence. CountInfoIndex provides a means to
     * keep track of multiple sequences.
     *
     * @param countInfoIndex Identifier of the sequence for which counts need to be recorded.
     * @return A ready CountWriter implementation.
     * @throws IOException If an error occurs.
     */
    public CountsWriterI newCountWriter(final int countInfoIndex) throws IOException {
        return newCountWriter(countInfoIndex, String.valueOf(countInfoIndex));
    }

    /**
     * Return a count Writer to the counts archive. Count writers must be returned to the
     * archive after they have been populated with count information.
     *
     * @param writerI The countWriter being returned.
     * @throws IOException If an error occurs packaging the count information in the archive.
     */
    public void returnWriter(final CountsWriterI writerI) throws IOException {
        assert writerI == currentCountsWriterI : "You must return the current counts writer.";
        writerI.close();
        totalBitsWritten += writerI.getNumberOfBitsWritten();
        totalTransitions += writerI.getNumberOfTransitions();
        totalBasesSeen += writerI.getNumberOfBasesSeen();
        totalSitesSeen += writerI.getNumberOfSitesSeen();
        final byte[] bytes = stream.toByteArray();

        final org.bdval.io.compound.CompoundDataOutput part = compoundWriter.addFile(currentId);
        part.write(bytes);
        part.close();
        compoundWriter.finishAddFile();

        final org.bdval.io.compound.CompoundDataOutput indexPart = compoundWriter.addFile("#index:" + currentId);
        indexBuilder.buildIndex(bytes, indexPart);
        indexPart.close();
        compoundWriter.finishAddFile();
        currentId = null;
        currentCountsWriterI = null;
    }

    /**
     * Closes this archive.
     *
     * @throws IOException
     */
    @Override
    public void close() throws IOException {
        // a file whose name starts with # is special and not interpreted as a reference sequence index,id
        final CompoundDataOutput statsFile = compoundWriter.addFile("#stats");
        statsFile.writeUTF("totalBasesSeen");
        statsFile.writeLong(totalBasesSeen);
        statsFile.writeUTF("totalSitesSeen");
        statsFile.writeLong(totalSitesSeen);
        statsFile.writeUTF("END");
        statsFile.close();
        compoundWriter.finishAddFile();
        if (verbose) {
            System.out.println("Global statististics:%n");
            System.out.printf("Bits written: %d %n", totalBitsWritten);
            System.out.printf("Number of transitions: %d %n", totalTransitions);
            System.out.printf("Bits per transitions: %2.2g %n", ((double) (totalBitsWritten) / (double) totalTransitions));
        }
        compoundWriter.close();
    }

    /**
     * Set the verbose flag. When the flag is true, prints some statistics to the console upon close(). Otherwise, no
     * output to the console.
     *
     * @param verbose Verbose flag.
     */
    public void setVerbose(final boolean verbose) {
        this.verbose = verbose;
    }
}
