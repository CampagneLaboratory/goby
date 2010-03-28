/*
 * Copyright (C) 2010 Institute for Computational Biomedicine,
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

import it.unimi.dsi.lang.MutableString;

/**
 * Interface for global normalization methods.
 *
 * @author Fabien Campagne
 *         Date: Mar 27, 2010
 *         Time: 3:02:34 PM
 */
public interface NormalizationMethod {
    /**
     * Return the name of the method.
     *
     * @return a string that uniquely identifies this method.
     */
    public String getIdentifier();

    /**
     * Return an abbreviation to be used as prefix to statistics normalized with this method. All upper-case three letters
     * abbreviations are strongly suggested to make the output look consistent.
     *
     * @return
     */
    public String getAbbreviation();

    void normalize(DifferentialExpressionCalculator calculator, String... groups);


    double getNormalizedExpressionValue(final DifferentialExpressionCalculator deCalculator,
                                        final String sampleId,
                                        final MutableString elementId);
}
