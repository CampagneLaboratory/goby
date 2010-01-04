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

package edu.cornell.med.icb.goby.counts;

import org.apache.commons.io.FileUtils;
import org.bdval.io.compound.CompoundFileWriter;

import java.io.ByteArrayOutputStream;
import java.io.Closeable;
import java.io.DataOutput;
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
    private CountsWriter currentCountsWriter;
    private String currentId;

    private ByteArrayOutputStream stream;

    /**
     * Initialize with a basename. Count information will be written to a file basename+".counts"
     *
     * @param basename Basename for the alignment/counts.
     * @throws IOException If an error occurs preparing the packaged counts.
     */
    public CountsArchiveWriter(final String basename) throws IOException {
      this(basename, "counts");
    }
    /**
     * Initialize with a basename. Count information will be written to a file basename+"."+countArchiveModifier
     *
     * @param basename Basename for the alignment/counts.
     * @param countArchiveModifier the extension to write the archive.
     * @throws IOException If an error occurs preparing the packaged counts.
     */
    public CountsArchiveWriter(String basename, String countArchiveModifier) throws IOException {
        final String physicalFilename = basename + "."+countArchiveModifier;
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
     * @param identifier Identifier of the sequence for which counts need to be recorded.
     * @return A ready CountWriter implementation.
     * @throws IOException If an error occurs.
     */
    public CountsWriter newCountWriter(final int referenceIndex, final String identifier) throws IOException {
        stream = new ByteArrayOutputStream(100000);
        currentCountsWriter = new CountsWriter(stream);
        currentId = Integer.toString(referenceIndex)+","+identifier;
        return currentCountsWriter;
    }

    /**
     * Obtain a countsWriter to write counts for a sequence. CountInfoIndex provides a means to
     * keep track of multiple sequences.
     *
     * @param countInfoIndex Identifier of the sequence for which counts need to be recorded.
     * @return A ready CountWriter implementation.
     * @throws IOException If an error occurs.
     */
    public CountsWriter newCountWriter(final int countInfoIndex) throws IOException {
        return newCountWriter(countInfoIndex,String.valueOf(countInfoIndex));
    }

    /**
     * Return a count Writer to the counts archive. Count writers must be returned to the
     * archive after they have been populated with count information.
     *
     * @param writer The countWriter being returned.
     * @throws IOException If an error occurs packaging the count information in the archive.
     */
    public void returnWriter(final CountsWriter writer) throws IOException {
        assert writer == currentCountsWriter : "You must return the current counts writter.";
        writer.close();
        final byte[] bytes = stream.toByteArray();

        final DataOutput part = compoundWriter.addFile(currentId);
        part.write(bytes);

        compoundWriter.finishAddFile();
        currentId = null;
        currentCountsWriter = null;
    }

    /**
     * Closes this archive.
     *
     * @throws IOException
     */
    public void close() throws IOException {
        compoundWriter.close();
    }
}
