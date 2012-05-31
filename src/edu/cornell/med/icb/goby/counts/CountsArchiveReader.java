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

import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.io.FastByteArrayInputStream;
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.apache.log4j.Logger;
import org.bdval.io.compound.CompoundDataInput;
import org.bdval.io.compound.CompoundDirectoryEntry;
import org.bdval.io.compound.CompoundFileReader;

import java.io.*;
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
    protected final CompoundFileReader compoundReader;
    private final Int2ObjectMap<String> indexToIdentifierMap;
    private Object2IntMap<String> identifierToIndexMap;

    /**
     * The total number of bases seen in the counts data stored in this archive.
     * This is defined as the sum of count*length over all transitions stored. A normalization
     * factor for count data can be defined as   getNumberOfBasesSeen()/  getNumberOfSitesSeen() : this represents
     * the average coverage per site observed.
     *
     * @return the total number of bases seen.
     */


    public long getTotalBasesSeen() {

        return totalBasesSeen;
    }

    /**
     * The total number of sites observed at which count!=0.
     *
     * @return number of sites seen in the archive.
     */
    public long getTotalSitesSeen() {
        return totalSitesSeen;
    }

    /**
     * Indicate whether the archive had statistics.  Count archives generated with Goby 1.9.7+ include normalization statistics calculated
     * at the time the counts were written.
     *
     * @return whether the archive had statistics.
     */
    public boolean isStatsParsed() {
        return totalBasesSeen != 0 && totalSitesSeen != 0;
    }

    private long totalBasesSeen;
    private long totalSitesSeen;

    /**
     * Initialize the MultiCountReader. Will look for count archive information in basename".count"
     *
     * @param basename Basename associated with the count information archive. Opens the defaults
     *                 ".counts" archive.
     * @throws IOException if the file cannot be accessed
     */
    public CountsArchiveReader(final String basename) throws IOException {
        this(basename, COUNT_ARCHIVE_MODIFIER_DEFAULT);
    }

    final static String COUNT_ARCHIVE_MODIFIER_DEFAULT = "counts";

    /**
     * Initialize the MultiCountReader. Will look for count archive information in basename".count"
     *
     * @param basename Basename associated with the count information archive.
     * @param alternativeCountArchiveExtension
     *                 Which extension to use to open the count archive.
     * @throws IOException if the file cannot be accessed
     */
    public CountsArchiveReader(final String basename, final String alternativeCountArchiveExtension) throws IOException {
        compoundReader = new CompoundFileReader(basename + "." + alternativeCountArchiveExtension);
        indexToIdentifierMap = new Int2ObjectOpenHashMap<String>();
        identifierToIndexMap = new Object2IntOpenHashMap<String>();
        scanDirectory();
    }


    /**
     * Obtain the count reader over the count information identified by id.
     *
     * @param identifier The identifier of the count information to retrieve.
     * @return A countReader over count information associated with the id.
     * @throws IOException If an error occurs reading count information.
     */
    public CountsReader getCountReader(final String identifier) throws IOException {

        String name = makeFileIdentifier(identifier);
        final CompoundDataInput input = compoundReader.readFile(name);
        final byte[] bytes = new byte[(int) input.length()];
        input.readFully(bytes);
        // warning: the countStream implementation has to support RepositionableStream
        final InputStream countStream = new FastByteArrayInputStream(bytes);
        final String indexName = "#index:" + name;
        if (compoundReader.containsFile(indexName)) {

            final CompoundDataInput indexInput = compoundReader.readFile(indexName);
            // this archive contained an index:
            return new CountsReader(countStream, indexInput);
        } else {
            return new CountsReader(countStream);
        }

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

    private String makeFileIdentifier(final int countInfoIndex) {
        return countInfoIndex + "," + indexToIdentifierMap.get(countInfoIndex);
    }

    public String getIdentifier(final int index) {
        return indexToIdentifierMap.get(index);
    }

    protected String makeFileIdentifier(final String countId) {
        return countId.indexOf(',') == -1
                ? identifierToIndexMap.get(countId) + "," + countId : countId;

    }

    /**
     * Obtain the number of count information available in this archive.
     *
     * @return an integer.
     */

    public int getNumberOfIndices() {
        return identifierToIndexMap.size();

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

            if (!name.startsWith("#")) {
                final String[] tokens = name.split(",");
                if (tokens.length < 2) {
                    return null;
                }
                result.add(tokens[1]);
            }
        }
        return result;
    }

    private void scanDirectory() {
        final Collection<CompoundDirectoryEntry> directory = compoundReader.getDirectory();
        for (final CompoundDirectoryEntry entry : directory) {
            final String name = entry.getName();
            if (!name.startsWith("#")) {
                final String[] tokens = name.split(",");
                assert tokens.length == 2 : "archive count filenames must be of the form int,String";
                final int index = Integer.parseInt(tokens[0]);
                final String id = tokens[1];
                indexToIdentifierMap.put(index, id);
                identifierToIndexMap.put(id, index);
            } else {
                parseSpecialFile(entry);
            }

        }

    }

    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(CountsArchiveReader.class);

    private void parseSpecialFile(CompoundDirectoryEntry entry) {
        if ("#stats".equals(entry.getName())) {
            try {

                final CompoundDataInput input = compoundReader.readFile("#stats");
                boolean done = false;
                while (!done) {
                    String key = input.readUTF();
                    if ("totalBasesSeen".equals(key)) {
                        totalBasesSeen = input.readLong();
                    } else {
                        if ("totalSitesSeen".equals(key)) {
                            totalSitesSeen = input.readLong();
                        }
                    }
                    if ("END".equals(key)) {
                        done = true;
                    }
                }
            } catch (IOException e) {
                LOG.error("could not access special file #stats in counts archive");
            }

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

            if (!name.startsWith("#")) {
                final String[] tokens = name.split(",");

                if (tokens.length < 1) {
                    return null;
                }

                result.add(Integer.parseInt(tokens[0]));
            }
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
