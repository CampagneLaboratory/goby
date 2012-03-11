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

import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: 3/11/12
 *         Time: 3:54 PM
 */
public interface PermutationReaderInterface {
    /**
     * Return the query index associated with a small index, or -1 if the association was not defined.
     *
     * @param smallIndex
     * @return
     * @throws java.io.IOException
     */
    int getQueryIndex(int smallIndex) throws IOException;
}
