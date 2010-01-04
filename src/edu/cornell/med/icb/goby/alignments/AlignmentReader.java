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

package edu.cornell.med.icb.goby.alignments;

import com.google.protobuf.CodedInputStream;
import edu.cornell.med.icb.goby.reads.FastBufferedMessageChunksReader;
import edu.cornell.med.icb.goby.reads.MessageChunksReader;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.io.FastBufferedInputStream;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.Closeable;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.util.Iterator;
import java.util.Properties;
import java.util.zip.GZIPInputStream;

/**
 * Reads alignments written with AlignmentWriter.
 *
 * @author Fabien Campagne
 *         Date: Apr 30, 2009
 *         Time: 6:36:04 PM
 */
public class AlignmentReader extends AbstractAlignmentReader implements Closeable,
        Iterator<Alignments.AlignmentEntry>, Iterable<Alignments.AlignmentEntry> {
    private FileInputStream headerStream;
    private static final Log LOG = LogFactory.getLog(AlignmentReader.class);
    private int numberOfAlignedReads;
    private final MessageChunksReader alignmentEntryReader;
    private Alignments.AlignmentCollection collection;
    private Properties stats;

    public AlignmentReader(final String basename) throws FileNotFoundException {
        final FileInputStream stream = new FileInputStream(basename + ".entries");
        alignmentEntryReader = new MessageChunksReader(stream);
        headerStream = new FileInputStream(basename + ".header");
        stats = new Properties();
        try {
            final String filename = basename + ".stats";
            final File statsFile = new File(filename);
            if (statsFile.exists()) {
                stats.load(new FileReader(statsFile));
            }
        } catch (IOException e) {
            LOG.warn("cannot load properties for basename: " + basename);
        }
    }

    public AlignmentReader(final InputStream entriesStream) {
        alignmentEntryReader = new MessageChunksReader(entriesStream);
    }

    /**
     * Initialize the reader to read a segment of the input. Sequences represented by a collection which
     * starts between the input position start and end will be returned upon subsequent calls to hasNext(), next().
     *
     * @param start  Start offset in the input file.
     * @param end    End offset in the input file.
     * @param stream Stream over the input file
     * @throws IOException If an error occurs reading the input.
     */
    public AlignmentReader(final long start, final long end, final FastBufferedInputStream stream) throws IOException {
        alignmentEntryReader = new FastBufferedMessageChunksReader(start, end, stream);

    }

    /**
     * Returns true if the input has more sequences.
     *
     * @return true if the input has more sequences, false otherwise.
     */
    @Override
    public boolean hasNextAligmentEntry() {
        final boolean hasNext = alignmentEntryReader.hasNext(collection, numberOfEntries());
        if (!hasNext) {
            collection = null;
        }
        final GZIPInputStream uncompressStream = alignmentEntryReader.getUncompressStream();
        try {
            if (uncompressStream != null) {
                collection = Alignments.AlignmentCollection.parseFrom(uncompressStream);
                if (collection.getAlignmentEntriesCount() == 0) {
                    return false;
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return hasNext;

    }

    private int numberOfEntries() {
        return collection != null ? collection.getAlignmentEntriesCount() : 0;
    }

    /**
     * Returns the next read entry from the input stream.
     *
     * @return the next read entry from the input stream.
     */
    @Override
    public final Alignments.AlignmentEntry nextAlignmentEntry() {

        if (!alignmentEntryReader.hasNext(collection, numberOfEntries())) {
            throw new IllegalStateException();
        }
        return collection.getAlignmentEntries(alignmentEntryReader.incrementEntryIndex());
    }

    public boolean hasNext() {
        return hasNextAligmentEntry();
    }

    public Alignments.AlignmentEntry next() {
        return nextAlignmentEntry();
    }

    /**
     * This operation is not supported.
     */
    public void remove() {
        throw new UnsupportedOperationException("Cannot remove from a reader.");
    }

    /**
     * Read the header of this alignment.
     *
     * @throws IOException If an error occurs.
     */
    @Override
    public void readHeader() throws IOException {
        // accept very large header messages, since these may contain query identifiers:
        final CodedInputStream codedInput = CodedInputStream.newInstance(headerStream);
        codedInput.setSizeLimit(Integer.MAX_VALUE);
        final Alignments.AlignmentHeader header = Alignments.AlignmentHeader.parseFrom(codedInput);


        queryIdentifiers = parseIdentifiers(header.getQueryNameMapping());
        targetIdentifiers = parseIdentifiers(header.getTargetNameMapping());
        if (header.getQueryLengthCount() > 0) {
            queryLengths = new IntArrayList(header.getQueryLengthList()).toIntArray();
        }
        numberOfQueries = header.getNumberOfQueries();
        numberOfTargets = header.getNumberOfTargets();
        setHeaderLoaded(true);
        numberOfAlignedReads=header.getNumberOfAlignedReads();
    }

    private IndexedIdentifier parseIdentifiers(final Alignments.IdentifierMapping nameMapping) {
        final IndexedIdentifier result = new IndexedIdentifier();

        for (final Alignments.IdentifierInfo info : nameMapping.getMappingsList()) {
            result.put(new MutableString(info.getName()), info.getIndex());
        }
        return result;
    }

    /**
     * {@inheritDoc}
     */
    public void close() {
        alignmentEntryReader.close();
    }

    /**
     * Return the basename corresponding to the input filename.
     */
    public static String getBasename(final String inputFile) {
        for (final String ext : new String[]{".entries", ".tmh", ".header", ".counts", ".stats"}) {
            if (inputFile.endsWith(ext)) {
                return inputFile.replace(ext, "");
            }
        }
        // perhaps the input was a basename alread..
        return inputFile;
    }

    /**
     * Return the basenames corresponding to the input filenames. Less basename than filenames may be returned (if
     * several filenames reduce to the same baseline after removing the extension).
     */
    public static String[] getBasenames(final String... inputFiles) {
        final ObjectSet<String> result = new ObjectArraySet<String>();

        for (final String filename : inputFiles) {
            result.add(getBasename(filename));
        }
        return result.toArray(new String[result.size()]);
    }

    public Iterator<Alignments.AlignmentEntry> iterator() {
        return this;
    }

    public Properties getStatistics() {
        return stats;
    }

    public int getNumberOfAlignedReads() {
        return numberOfAlignedReads;
    }
}
