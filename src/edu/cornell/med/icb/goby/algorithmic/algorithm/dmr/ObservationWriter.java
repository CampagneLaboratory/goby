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

package edu.cornell.med.icb.goby.algorithmic.algorithm.dmr;

import edu.cornell.med.icb.goby.algorithmic.data.GroupComparison;
import it.unimi.dsi.fastutil.ints.IntArrayList;

import java.io.PrintWriter;
import java.io.Writer;

/**
 * Writes observations to disk in the format recognized by empirical-p mode.
 *
 * @author Fabien Campagne
 *         Date: 2/27/12
 *         Time: 12:14 PM
 */
public class ObservationWriter {
    private TypeOfPair typeOfPair;
    private boolean headerWritten = false;
    private String[] headerIds;

    public void close() {
        outputWriter.close();
    }


    public void setNullComparison(String nameGroup) {
        comparison=nameGroup+"/"+nameGroup;
    }

    public static enum TypeOfPair {
        UNDEFINED, WITHIN_GROUP_PAIR, BETWEEN_GROUP_PAIR
    }

    PrintWriter outputWriter;
    String[] elementIds;

    public ObservationWriter(Writer writer) {
        this.outputWriter = new PrintWriter(writer);
    }

    public void setElementIds(final String[] ids) {
        elementIds = ids;
    }

    public void setTypeOfPair(final TypeOfPair type) {
        this.typeOfPair = type;
    }

    public void setHeaderIds(String[] headerIds) {
        this.headerIds = headerIds;
    }

    public void writeHeader(String[] headerValuesA, String[] headerValuesB, String[] headerCovariatesA, String[] headerCovariatesB) {
        if (headerWritten) {
            return;
        } else {
            assert headerIds != null : " setHeaderIds must be called before calling writeHeader";

            outputWriter.write("PAIR_TYPE");
            outputWriter.write("\tCOMPARISON");
            writeHeaderColumns(headerIds);
            outputWriter.write("\tDELIMITER1");
            writeHeaderColumns(headerValuesA);
            outputWriter.write("\tDELIMITER2");
            writeHeaderColumns(headerValuesB);
            outputWriter.write("\tDELIMITER3");
            writeHeaderColumns(headerCovariatesA);
            outputWriter.write("\tDELIMITER4");
            writeHeaderColumns(headerCovariatesB);
            outputWriter.write('\n');
            headerWritten = true;
        }
    }

    private void writeHeaderColumns(String[] array) {
        for (String id : array) {
            outputWriter.write('\t');
            outputWriter.write(id);
        }
    }

    public void setComparison(GroupComparison comparison) {
        this.comparison = comparison.nameGroup1+"/"+comparison.nameGroup2;
    }

    String comparison;

    public void observed(final IntArrayList valuesA, final IntArrayList valuesB,
                         final IntArrayList covariatesA, final IntArrayList covariatesB) {

        outputWriter.write(typeOfPair.toString());
        outputWriter.write('\t');
        outputWriter.write(comparison);

        for (final String elementId : elementIds) {
            outputWriter.write("\t");
            outputWriter.write(elementId);

        }
        write("VALUES_A", valuesA);
        write("VALUES_B", valuesB);
        write("COVARIATES_A", covariatesA);
        write("COVARIATES_B", covariatesB);
        outputWriter.write("\n");

    }

    private void write(final String id, final IntArrayList vals) {
        outputWriter.write('\t');
        outputWriter.write(id);
        for (final int value : vals) {
            outputWriter.write('\t');
            outputWriter.write(Integer.toString(value));
        }
    }

}
