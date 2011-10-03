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

package edu.cornell.med.icb.goby.util;

import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.counts.CountsReader;
import edu.cornell.med.icb.goby.reads.ReadsReader;

import java.io.IOException;

/**
 * @author campagne
 *         Date: 9/30/11
 *         Time: 1:49 PM
 */
public class IOUtil {
    /**
     * Close quietly a reads reader.
     *
     * @param reader
     */
    public static void closeQuietly(ReadsReader reader) {
        if (reader==null) {
            return;
        }

        try {
            reader.close();

        } catch (IOException e) {
            return;
        }
    }

    /**
     * Close quietly a reads reader.
     *
     * @param reader
     */
    public static void closeQuietly(AlignmentReader reader) {
     if (reader==null) {
            return;
        }
        reader.close();

    }

    /**
     * Close quietly a counts reader.
     *
     * @param reader
     */
    public static void closeQuietly(CountsReader reader) {
        if (reader==null) {
            return;
        }
        try {
            reader.close();
        } catch (IOException e) {

        }

    }
}
