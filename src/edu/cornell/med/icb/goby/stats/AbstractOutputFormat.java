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

package edu.cornell.med.icb.goby.stats;

import edu.cornell.med.icb.goby.modes.DiscoverSequenceVariantsMode;
import edu.cornell.med.icb.goby.util.OutputInfo;
import edu.cornell.med.icb.goby.util.OutputInfoFromWriter;

import java.io.PrintWriter;

/**
 * @author Fabien Campagne
 *         Date: 2/20/12
 *         Time: 12:21 PM
 */
public abstract class AbstractOutputFormat {
    public void defineColumns(final PrintWriter writer, final DiscoverSequenceVariantsMode mode) {
        defineColumns(new OutputInfoFromWriter(writer), mode);
    }

    public abstract void defineColumns(final OutputInfo outputInfo, final DiscoverSequenceVariantsMode mode);

}
