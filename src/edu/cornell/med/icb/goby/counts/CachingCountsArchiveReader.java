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

package edu.cornell.med.icb.goby.counts;


import it.unimi.dsi.fastutil.io.FastByteArrayInputStream;
import org.bdval.io.compound.CompoundDataInput;

import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStream;

/**
 * A count archive reader that caches the last chromosome accessed.
 *
 * @author Fabien Campagne
 *         Date: 6/12/11
 *         Time: 2:09 PM
 */
public class CachingCountsArchiveReader extends CountsArchiveReader {

    private byte[] bytes;
    private String previousId;
    private byte[] indexBytes;
    private boolean hasIndex;

    public CachingCountsArchiveReader(final String basename) throws IOException {
        super(basename);
    }

    int previousIndex = -1;

    /**
     * Determine if the previous counts reader  returned supports position().
     * @return
     */
    public boolean hasIndex() {
        return hasIndex;
    }

    /**
     * Obtain the count reader over the count information identified by id.
     *
     * @param identifier The indentifier of the count information to retrieve.
     * @return A countReader over count information associated with the id.
     * @throws IOException If an error occurs reading count information.
     */
    @Override
    public CountsReader getCountReader(final String identifier) throws IOException {
        if (identifier == null || previousId == null || !identifier.equals(previousId)) {
            String name = makeFileIdentifier(identifier);
            final CompoundDataInput input = compoundReader.readFile(name);
            bytes = new byte[(int) input.length()];
            input.readFully(bytes);
            String indexName = "#index:" + name;
            if (compoundReader.containsFile(indexName)) {

                final CompoundDataInput indexInput = compoundReader.readFile(indexName);
                indexBytes = new byte[(int) indexInput.length()];
                indexInput.readFully(indexBytes);
                hasIndex = true;
            }

            previousId = identifier;
        }
        final InputStream stream = new FastByteArrayInputStream(bytes);
        if (hasIndex) {


            final DataInputStream indexDataInput = new DataInputStream(new FastByteArrayInputStream(indexBytes));
            return new CountsReader(stream, indexDataInput);

        } else {

            return new CountsReader(stream);

        }

    }
}
