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

package edu.cornell.med.icb.goby.modes;

import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import static org.junit.Assert.assertEquals;

/**
 * Test TabToColumnInfoMode.
 */
public class TestTabToColumnInfoMode {

    @Test
    public void testTabToColumnInfo() throws IOException {
        final TabToColumnInfoMode columnsInfo = new TabToColumnInfoMode();
        final String inputFilename = "test-data/tsv/col-info-test.tsv";
        columnsInfo.addInputFile(new File(inputFilename));
        columnsInfo.setNoOutput(true);
        columnsInfo.execute();
        final Map<String, TabToColumnInfoMode.ColumnTypes> columnDetailsMap = columnsInfo.getFilenameToDetailsMap().get(inputFilename);
        assertEquals("column type is wrong", "INT", columnDetailsMap.get("int").toString());
        assertEquals("column type is wrong", "DOUBLE", columnDetailsMap.get("double").toString());
        assertEquals("column type is wrong", "DOUBLE", columnDetailsMap.get("double2").toString());
        assertEquals("column type is wrong", "STRING", columnDetailsMap.get("string1").toString());
        assertEquals("column type is wrong", "STRING", columnDetailsMap.get("string2").toString());
        assertEquals("column type is wrong", "STRING", columnDetailsMap.get("string3").toString());
    }

    @Test
    public void testTabToColumnInfoGzip() throws IOException {
        final TabToColumnInfoMode columnsInfo = new TabToColumnInfoMode();
        final String inputFilename = "test-data/tsv/col-info-test.tsv.gz";
        columnsInfo.addInputFile(new File(inputFilename));
        columnsInfo.setNoOutput(true);
        columnsInfo.execute();
        final Map<String, TabToColumnInfoMode.ColumnTypes> columnDetailsMap = columnsInfo.getFilenameToDetailsMap().get(inputFilename);
        assertEquals("column type is wrong", "INT", columnDetailsMap.get("int").toString());
        assertEquals("column type is wrong", "DOUBLE", columnDetailsMap.get("double").toString());
        assertEquals("column type is wrong", "DOUBLE", columnDetailsMap.get("double2").toString());
        assertEquals("column type is wrong", "STRING", columnDetailsMap.get("string1").toString());
        assertEquals("column type is wrong", "STRING", columnDetailsMap.get("string2").toString());
        assertEquals("column type is wrong", "STRING", columnDetailsMap.get("string3").toString());
    }

    @Test
    public void testObserveValue() {
        assertEquals("observed type is wrong", "INT", TabToColumnInfoMode.typeFromValue("4").toString());
        assertEquals("observed type is wrong", "DOUBLE", TabToColumnInfoMode.typeFromValue("3.2").toString());
        assertEquals("observed type is wrong", "INT", TabToColumnInfoMode.typeFromValue("5").toString());
        assertEquals("observed type is wrong", "STRING", TabToColumnInfoMode.typeFromValue("3.2a").toString());
        assertEquals("observed type is wrong", "INT", TabToColumnInfoMode.typeFromValue("4").toString());
        assertEquals("observed type is wrong", "DOUBLE", TabToColumnInfoMode.typeFromValue("3.2").toString());
        assertEquals("observed type is wrong", "INT", TabToColumnInfoMode.typeFromValue("93").toString());
        assertEquals("observed type is wrong", "DOUBLE", TabToColumnInfoMode.typeFromValue("-2e2").toString());
        assertEquals("observed type is wrong", "DOUBLE", TabToColumnInfoMode.typeFromValue("-72.3e-02").toString());
        assertEquals("observed type is wrong", "DOUBLE", TabToColumnInfoMode.typeFromValue("+72.3e-2").toString());
        assertEquals("observed type is wrong", "DOUBLE", TabToColumnInfoMode.typeFromValue("-3.2").toString());
        assertEquals("observed type is wrong", "DOUBLE", TabToColumnInfoMode.typeFromValue("3.2").toString());
        assertEquals("observed type is wrong", "INT", TabToColumnInfoMode.typeFromValue("4").toString());
        assertEquals("observed type is wrong", "INT", TabToColumnInfoMode.typeFromValue("-32").toString());
        // Values are trimmed.
        assertEquals("observed type is wrong", "INT", TabToColumnInfoMode.typeFromValue("  -32  ").toString());
        assertEquals("observed type is wrong", "DOUBLE", TabToColumnInfoMode.typeFromValue("NaN").toString());
        assertEquals("observed type is wrong", "DOUBLE", TabToColumnInfoMode.typeFromValue("-Inf").toString());
        assertEquals("observed type is wrong", "STRING", TabToColumnInfoMode.typeFromValue("woot").toString());
        assertEquals("observed type is wrong", "DOUBLE", TabToColumnInfoMode.typeFromValue("+Infinity").toString());
        // Integers (parsed from String) cannot start with + but Doubles can. Strange but true.
        assertEquals("observed type is wrong", "DOUBLE", TabToColumnInfoMode.typeFromValue("+32").toString());
        assertEquals("observed type is wrong", "DOUBLE", TabToColumnInfoMode.typeFromValue("+INFINITY").toString());
        assertEquals("observed type is wrong", "DOUBLE", TabToColumnInfoMode.typeFromValue("-INF").toString());
        assertEquals("observed type is wrong", "STRING", TabToColumnInfoMode.typeFromValue("no").toString());
        assertEquals("observed type is wrong", "INT", TabToColumnInfoMode.typeFromValue("9").toString());
        assertEquals("observed type is wrong", "DOUBLE", TabToColumnInfoMode.typeFromValue("nan").toString());
        assertEquals("observed type is wrong", "DOUBLE", TabToColumnInfoMode.typeFromValue("+inf").toString());
        assertEquals("observed type is wrong", "DOUBLE", TabToColumnInfoMode.typeFromValue("inf").toString());
        assertEquals("observed type is wrong", "UNKNOWN", TabToColumnInfoMode.typeFromValue("    ").toString());
        assertEquals("observed type is wrong", "UNKNOWN", TabToColumnInfoMode.typeFromValue("").toString());
        assertEquals("observed type is wrong", "UNKNOWN", TabToColumnInfoMode.typeFromValue(null).toString());
    }
}
