/*
 * Copyright (C) 2009-2012 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.alignments.perms;

import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
import it.unimi.dsi.fastutil.objects.Object2ByteMap;
import it.unimi.dsi.fastutil.objects.Object2ByteOpenHashMap;
import it.unimi.dsi.fastutil.objects.Object2IntAVLTreeMap;
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

/**
 * A class to convert read names to query indices. Provides logic to keep a read name to index association
 * for as long as needed, but no longer (to reduce memory consumption). Used when converting SAM/BAM to Goby
 * alignments.
 *
 * @author Fabien Campagne
 *         Date: 3/5/12
 *         Time: 5:10 PM
 */
public class ReadNameToIndex {
    /**
     * Used to log informational and debug messages.
     */
    private static final Log LOG = LogFactory.getLog(ReadNameToIndex.class);

    private final String basename;
    private int globalQueryMaxOccurences = 2;


    public ReadNameToIndex(String filename) {
        this.basename = AlignmentReaderImpl.getBasename(filename);
        reset();

    }

    private void reset() {
        namesToIndex.clear();
        namesToIndex.defaultReturnValue(-1);
        timesRequested.clear();
        timesRequested.defaultReturnValue((byte)0);
    }

    private void pushToPreStorage(MutableString readName, int queryIndex) {
        //TODO write to some file.
 //       System.out.printf("%d\t%s%n",  queryIndex,readName);
    }

    private int smallIndexCounter;
    private final Object2ByteMap<MutableString> timesRequested = new Object2ByteOpenHashMap<MutableString>();
    private final Object2IntMap<MutableString> namesToIndex = new Object2IntAVLTreeMap<MutableString>();

    public int getQueryIndex(final String readName, final int maxObservations) {
        final MutableString readNameMutable = new MutableString(readName).compact();
        final int timesRequestedInt = timesRequested.getByte(readNameMutable);
        final int queryIndex;
        if (timesRequestedInt == 0) {
            queryIndex = smallIndexCounter++;
            if (maxObservations > 1) {
                namesToIndex.put(readNameMutable, queryIndex);
            }
        } else {
            queryIndex = namesToIndex.getInt(readNameMutable);
        }
        final int timesSeen = timesRequestedInt + 1;
        // decide if we have reached max observations for this query index:
        if (timesSeen >= maxObservations) {
            // if yes, remove the index from the map, it will not be asked again.
            namesToIndex.remove(readNameMutable);
            timesRequested.remove(readNameMutable);
            pushToPreStorage(readNameMutable, queryIndex);
        } else {
            // if not, keep it in the map until requested that many times.
            if (maxObservations > 1) {
                timesRequested.put(readNameMutable, (byte)timesSeen);
            }
        }
        return queryIndex;

    }


    public void setPruneLimit(byte limit) {
        globalQueryMaxOccurences = limit;
    }


    public void close() {
        // TODO close output writer.
    }
}
