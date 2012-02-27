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

import edu.cornell.med.icb.goby.algorithmic.algorithm.dmr.ObservationWriter;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import org.apache.commons.io.output.NullWriter;

import java.io.Writer;

/**
 * This implementation does nothing.
 *
 * @author Fabien Campagne
 *         Date: 2/27/12
 *         Time: 12:38 PM
 */

public class DummyObservationWriter extends ObservationWriter {
    @Override
    public void close() {

    }

    public DummyObservationWriter() {
        super(new NullWriter());

    }

    @Override
    public void setElementIds(String[] ids) {
    }

    @Override
    public void setTypeOfPair(TypeOfPair type) {

    }

    @Override
    public void setHeaderIds(String[] headerIds) {

    }

    @Override
    public void writeHeader(String[] headerValuesA, String[] headerValuesB, String[] headerCovariatesA, String[] headerCovariatesB) {
    }

    @Override
    public void observed(IntArrayList valuesA, IntArrayList valuesB, IntArrayList covariatesA, IntArrayList covariatesB) {

    }
}
