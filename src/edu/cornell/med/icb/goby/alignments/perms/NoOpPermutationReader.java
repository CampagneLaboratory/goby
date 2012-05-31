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
 * Use this class when there is no permutation.
 * @author Fabien Campagne
 *         Date: 3/11/12
 *         Time: 3:53 PM
 */
public class NoOpPermutationReader implements PermutationReaderInterface {

    @Override
    /**
     * This class always return the value provided as argument.
     * @return i
     */
    public final int getQueryIndex(final int i) throws IOException {
        return i;
    }
}
