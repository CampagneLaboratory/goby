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

import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.bdval.io.compound.CompoundDataInput;
import org.bdval.io.compound.CompoundDirectoryEntry;
import org.bdval.io.compound.CompoundFileReader;

import java.io.ByteArrayInputStream;
import java.io.Closeable;
import java.io.IOException;
import java.io.InputStream;
import java.util.Collection;

/**
 * Read an archive of count information. Multiple sequences are typically stored in an archive.
 * This reader provides support to obtain count information for each sequence.
 *
 * @author Fabien Campagne
 *         Date: May 14, 2009
 *         Time: 12:22:30 PM
 */
public class CountsArchiveReader implements Closeable {
    private final CompoundFileReader compoundReader;
    private final Int2ObjectMap<String> indexToIdentifierMap;
    private Object2IntMap<String> identifierToIndexMap;

    /**
     * Initialize the MultiCountReader. Will look for count archive information in basename".count"
     *
     * @param basename Basename associated with the count information archive. Opens the defaults
     *                 ".counts" archive.
     * @throws IOException if the file cannot be accessed
     */
    public CountsArchiveReader(final String basename) throws IOException {
        this(basename, "counts");
    }

    /**
     * Initialize the MultiCountReader. Will look for count archive information in basename".count"
     *
     * @param basename Basename associated with the count information archive.
     * @param alternativeCountArchiveExtension
     *                 Which extension to use to open the count archive.
     * @throws IOException if the file cannot be accessed
     */
    public CountsArchiveReader(String basename, String alternativeCountArchiveExtension) throws IOException {
        compoundReader = new CompoundFileReader(basename + "." + alternativeCountArchiveExtension);
        indexToIdentifierMap = new Int2ObjectOpenHashMap<String>();
        identifierToIndexMap = new Object2IntOpenHashMap<String>();
        scanDirectory();
    }


    /**
     * Obtain the count reader over the count information identified by id.
     *
     * @param identifier The indentifier of the count information to retrieve.
     * @return A countReader over count information associated with the id.
     * @throws IOException If an error occurs reading count information.
     */
    public CountsReader getCountReader(final String identifier) throws IOException {

        final CompoundDataInput input = compoundReader.readFile(makeFileIdentifier(identifier));
        final byte[] bytes = new byte[(int) input.length()];
        input.readFully(bytes);
        final InputStream stream = new ByteArrayInputStream(bytes);
        return new CountsReader(stream);

    }

    /**
     * Obtain the count reader over the count information identified by index.
     * See getNumberOfIndices to determine the number of indices stored. Indices start at zero.
     *
     * @param countInfoIndex The index of the count information to retrieve.
     * @return A countReader over count information associated with the index.
     * @throws IOException If an error occurs reading count information.
     */
    public CountsReader getCountReader(final int countInfoIndex) throws IOException {
        return getCountReader(makeFileIdentifier(countInfoIndex));
    }

    private String makeFileIdentifier(int countInfoIndex) {
        return String.valueOf(countInfoIndex) + "," + indexToIdentifierMap.get(countInfoIndex);
    }

    public String getIdentifier(int index) {
        return indexToIdentifierMap.get(index);
    }

    private String makeFileIdentifier(String countId) {
        return countId.indexOf(',') == -1 ?
                String.valueOf(identifierToIndexMap.get(countId)) + "," + countId
                : countId;

    }

    /**
     * Obtain the number of count information available in this archive.
     *
     * @return an integer.
     */

    public int getNumberOfIndices() {
        return compoundReader.getDirectory().size();

    }

    /**
     * Obtain identifiers associated with each countReader defined in the archive.
     *
     * @return set of count information identifiers, or null if identifiers were not written in the archive.
     */
    public ObjectSet<String> getIdentifiers() {
        final ObjectSet<String> result = new ObjectOpenHashSet<String>();
        final Collection<CompoundDirectoryEntry> directory = compoundReader.getDirectory();
        for (final CompoundDirectoryEntry entry : directory) {
            final String name = entry.getName();
            final String tokens[] = name.split(",");
            if (tokens.length < 2) return null;
            result.add(tokens[1]);
        }
        return result;
    }

    private void scanDirectory() {

        final Collection<CompoundDirectoryEntry> directory = compoundReader.getDirectory();
        for (final CompoundDirectoryEntry entry : directory) {
            final String name = entry.getName();
            final String tokens[] = name.split(",");
            assert tokens.length == 2 : "archive count filenames must be of the form int,String";
            final int index = Integer.parseInt(tokens[0]);
            final String id = tokens[1];
            indexToIdentifierMap.put(index, id);
            identifierToIndexMap.put(id, index);
        }

    }

    /**
     * Obtain indices associated with each countReader defined in the archive. These indices are the index
     * of each reference sequence for which counts where written to this archive.
     *
     * @return set of reference sequence indices.
     */
    public IntSet getIndices() {
        final IntSet result = new IntArraySet();
        final Collection<CompoundDirectoryEntry> directory = compoundReader.getDirectory();
        for (final CompoundDirectoryEntry entry : directory) {
            final String name = entry.getName();
            final String tokens[] = name.split(",");
            if (tokens.length < 1) return null;


            result.add(Integer.parseInt(tokens[0]));
        }
        return result;
    }

    /**
     * Closes this archive.
     *
     * @throws IOException
     */
    public void close() throws IOException {
        compoundReader.close();
    }
}
