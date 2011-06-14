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

import org.bdval.io.compound.CompoundDataInput;

import java.io.ByteArrayInputStream;
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

    public CachingCountsArchiveReader(final String basename) throws IOException {
        super(basename);
    }

    int previousIndex = -1;


    /**
     * Obtain the count reader over the count information identified by id.
     *
     * @param identifier The indentifier of the count information to retrieve.
     * @return A countReader over count information associated with the id.
     * @throws IOException If an error occurs reading count information.
     */
    @Override
    public CountsReader getCountReader(final String identifier) throws IOException {
        if (identifier==null || previousId==null || !identifier.equals(previousId)) {
            final CompoundDataInput input = compoundReader.readFile(makeFileIdentifier(identifier));
            bytes = new byte[(int) input.length()];
            input.readFully(bytes);
            previousId=identifier;
        }

        final InputStream stream = new ByteArrayInputStream(bytes);
        return new CountsReader(stream);

    }
}
