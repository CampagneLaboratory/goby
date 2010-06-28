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

package edu.cornell.med.icb.goby.reads;

import com.google.protobuf.ByteString;
import com.google.protobuf.CodedInputStream;
import edu.cornell.med.icb.goby.exception.GobyRuntimeException;
import edu.cornell.med.icb.goby.util.FileExtensionHelper;
import it.unimi.dsi.fastutil.io.FastBufferedInputStream;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang.StringUtils;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.zip.GZIPInputStream;

/**
 * Reads sequences in the compact format from a stream produced with MessageChunkWriter.
 *
 * @author Fabien Campagne
 *         Date: Apr 24, 2009
 *         Time: 6:44:29 PM
 */
public class ReadsReader implements Iterator<Reads.ReadEntry>, Iterable<Reads.ReadEntry>,
        Closeable {
    private final MessageChunksReader reader;
    private Reads.ReadCollection collection;

    /**
     * Initialize the reader.
     *
     * @param path Path to the input file
     * @throws IOException If an error occurs reading the input
     */
    public ReadsReader(final String path) throws IOException {
        this(FileUtils.openInputStream(new File(path)));
    }

    /**
     * Initialize the reader.
     *
     * @param file The input file
     * @throws IOException If an error occurs reading the input
     */
    public ReadsReader(final File file) throws IOException {
        this(FileUtils.openInputStream(file));
    }

    /**
     * Initialize the reader.
     *
     * @param stream Stream over the input
     */
    public ReadsReader(final InputStream stream) {
        super();
        reader = new MessageChunksReader(stream);
    }

    /**
     * Initialize the reader to read a segment of the input. Sequences represented by a
     * collection which starts between the input position start and end will be returned
     * upon subsequent calls to {@link #hasNext()} and {@link #next()}.
     *
     * @param start Start offset in the input file
     * @param end   End offset in the input file
     * @param path  Path to the input file
     * @throws IOException If an error occurs reading the input
     */
    public ReadsReader(final long start, final long end, final String path) throws IOException {
        this(start, end, new FastBufferedInputStream(FileUtils.openInputStream(new File(path))));
    }

    /**
     * Initialize the reader to read a segment of the input. Sequences represented by a
     * collection which starts between the input position start and end will be returned
     * upon subsequent calls to {@link #hasNext()} and {@link #next()}.
     *
     * @param start  Start offset in the input file
     * @param end    End offset in the input file
     * @param stream Stream over the input file
     * @throws IOException If an error occurs reading the input.
     */
    public ReadsReader(final long start, final long end, final FastBufferedInputStream stream)
            throws IOException {
        super();
        reader = new FastBufferedMessageChunksReader(start, end, stream);
    }

    /**
     * Returns true if the input has more sequences.
     *
     * @return true if the input has more sequences, false otherwise.
     */
    public boolean hasNext() {
        final boolean hasNext =
                reader.hasNext(collection, collection != null ? collection.getReadsCount() : 0);
        final GZIPInputStream uncompressStream = reader.getUncompressStream();
        try {
            if (uncompressStream != null) {
                final CodedInputStream codedInput = CodedInputStream.newInstance(uncompressStream);
                codedInput.setSizeLimit(Integer.MAX_VALUE);
                collection = Reads.ReadCollection.parseFrom(codedInput);
                if (collection.getReadsCount() == 0) {
                    return false;
                }
            }
        } catch (IOException e) {
            throw new GobyRuntimeException(e);
        } finally {
            IOUtils.closeQuietly(uncompressStream);
        }
        return hasNext;
    }

    /**
     * Returns the next read entry from the input stream.
     * TODO: The current implementation will throw an exception if this is called before hasNext
     *
     * @return the next read entry from the input stream.
     */
    public final Reads.ReadEntry next() {
        if (!reader.hasNext(collection, collection.getReadsCount())) {
            throw new NoSuchElementException();
        }
        return collection.getReads(reader.incrementEntryIndex());
    }

    /**
     * This operation is not supported.
     */
    public void remove() {
        throw new UnsupportedOperationException("Cannot remove from a reader.");
    }
    /**
     * Decode the sequence in this entry to the sequence MutableString.
     *
     * @param entry    The entry which provides the sequence in encoded format.
     * @param sequence Where to write the decoded sequence.
     */
    public static void decodeSequence(final Reads.ReadEntry entry, final MutableString sequence) {
        decodeSequence(entry, sequence, false);
    }
    /**
     * Decode the sequence in this entry to the sequence MutableString.
     *
     * @param entry    The entry which provides the sequence in encoded format.
     * @param sequence Where to write the decoded sequence.
     * @param decodePair True: decodes the pair sequence. False: decodes the primary sequence.
     */
    public static void decodeSequence(final Reads.ReadEntry entry, final MutableString sequence, boolean decodePair) {
        final ByteString seq = decodePair ? entry.getSequencePair() : entry.getSequence();
        final int length =decodePair ? entry.getReadLengthPair(): entry.getReadLength();
        sequence.setLength(length);
        for (int i = 0; i < length; ++i) {
            sequence.setCharAt(i, (char) seq.byteAt(i));
        }
    }

    /**
     * Decode the quality scores in this entry to qualityScores MutableString.
     *
     * @param entry The entry which provides the sequence in encoded format.
     * @return the quality scores byte array (or an empty byte array if no scores)
     */
    public static byte[] decodeQualityScores(final Reads.ReadEntry entry) {
        final ByteString scores = entry.getQualityScores();
        if (scores == null) {
            return null;
        }

        final byte[] result = new byte[scores.size()];
        final int length = entry.getReadLength();
        for (int i = 0; i < length; ++i) {
            result[i] = scores.byteAt(i);
        }
        return result;
    }

    /**
     * Make the reader "iterable" for java "for each" loops and such.
     *
     * @return this object
     */
    public Iterator<Reads.ReadEntry> iterator() {
        return this;
    }

    /**
     * {@inheritDoc}
     */
    public void close() throws IOException {
        reader.close();
    }

    /**
     * Return the basename corresponding to the input reads filename.  Note
     * that if the filename does have the extension known to be a compact read
     * the returned value is the original filename
     *
     * @param filename The name of the file to get the basename for
     * @return basename for the alignment file
     */
    public static String getBasename(final String filename) {
        for (final String ext : FileExtensionHelper.COMPACT_READS_FILE_EXTS) {
            if (StringUtils.endsWith(filename, ext)) {
                return StringUtils.removeEnd(filename, ext);
            }
        }

        // perhaps the input was a basename already.
        return filename;
    }

    /**
     * Return the basenames corresponding to the input filenames. Less basename than filenames
     * may be returned (if several filenames reduce to the same baseline after removing
     * the extension).
     *
     * @param filenames The names of the files to get the basnames for
     * @return An array of basenames
     */
    public static String[] getBasenames(final String... filenames) {
        final ObjectSet<String> result = new ObjectArraySet<String>();
        if (filenames != null) {
            for (final String filename : filenames) {
                result.add(getBasename(filename));
            }
        }
        return result.toArray(new String[result.size()]);
    }


}
