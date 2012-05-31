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

import edu.cornell.med.icb.goby.readers.vcf.ColumnType;
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
        columnsInfo.addInputFilename(inputFilename);
        columnsInfo.setCreateCache(false);
        columnsInfo.execute();
        final Map<String, ColumnType> columnDetailsMap = columnsInfo.getFilenameToDetailsMap().get(inputFilename);
        assertEquals("column type is wrong", ColumnType.Integer, columnDetailsMap.get("int"));
        assertEquals("column type is wrong", ColumnType.Float, columnDetailsMap.get("double"));
        assertEquals("column type is wrong", ColumnType.Float, columnDetailsMap.get("double2"));
        assertEquals("column type is wrong", ColumnType.String, columnDetailsMap.get("string1"));
        assertEquals("column type is wrong", ColumnType.String, columnDetailsMap.get("string2"));
        assertEquals("column type is wrong", ColumnType.String, columnDetailsMap.get("string3"));
    }

    @Test
    public void testTabToColumnInfoGzip() throws IOException {
        final TabToColumnInfoMode columnsInfo = new TabToColumnInfoMode();
        final String inputFilename = "test-data/tsv/col-info-test.tsv.gz";
        columnsInfo.addInputFile(new File(inputFilename));
        columnsInfo.setCreateCache(false);
        columnsInfo.execute();
        final Map<String, ColumnType> columnDetailsMap = columnsInfo.getDetailsAtIndex(0);
        assertEquals("column type is wrong", ColumnType.Integer, columnDetailsMap.get("int"));
        assertEquals("column type is wrong", ColumnType.Float, columnDetailsMap.get("double"));
        assertEquals("column type is wrong", ColumnType.Float, columnDetailsMap.get("double2"));
        assertEquals("column type is wrong", ColumnType.String, columnDetailsMap.get("string1"));
        assertEquals("column type is wrong", ColumnType.String, columnDetailsMap.get("string2"));
        assertEquals("column type is wrong", ColumnType.String, columnDetailsMap.get("string3"));
    }

    @Test
    public void testTabToColumnInfoMultiIndex() throws IOException {
        final TabToColumnInfoMode columnsInfo = new TabToColumnInfoMode();
        final String inputFilename1 = "test-data/tsv/col-info-test.tsv";
        final String inputFilename2 = "test-data/tsv/col-info-test.tsv.gz";
        columnsInfo.addInputFilename(inputFilename1);
        columnsInfo.addInputFile(new File(inputFilename2));
        columnsInfo.setCreateCache(false);
        columnsInfo.execute();
        final Map<String, ColumnType> columnDetailsMap1 = columnsInfo.getDetailsAtIndex(0);
        final Map<String, ColumnType> columnDetailsMap2 = columnsInfo.getDetailsAtIndex(1);
    }
    
    @Test(expected=IndexOutOfBoundsException.class)
    public void testTabToColumnInfoMultiBadIndex() throws IOException {
        final TabToColumnInfoMode columnsInfo = new TabToColumnInfoMode();
        final String inputFilename1 = "test-data/tsv/col-info-test.tsv";
        final String inputFilename2 = "test-data/tsv/col-info-test.tsv.gz";
        columnsInfo.addInputFile(new File(inputFilename1));
        columnsInfo.addInputFile(new File(inputFilename2));
        columnsInfo.setCreateCache(false);
        columnsInfo.execute();
        final Map<String, ColumnType> columnDetailsMap = columnsInfo.getDetailsAtIndex(2);
    }

    @Test
    public void testObserveValue() {
        final TabToColumnInfoMode tabToColumnInfoMode = new TabToColumnInfoMode();
        tabToColumnInfoMode.setVerbose(true);
        assertEquals("observed type is wrong", ColumnType.Integer, tabToColumnInfoMode.typeFromValue("4"));
        assertEquals("observed type is wrong", ColumnType.Float, tabToColumnInfoMode.typeFromValue("3.2"));
        assertEquals("observed type is wrong", ColumnType.Integer, tabToColumnInfoMode.typeFromValue("5"));
        assertEquals("observed type is wrong", ColumnType.String, tabToColumnInfoMode.typeFromValue("3.2a"));
        assertEquals("observed type is wrong", ColumnType.Integer, tabToColumnInfoMode.typeFromValue("4"));
        assertEquals("observed type is wrong", ColumnType.Float, tabToColumnInfoMode.typeFromValue("3.2"));
        assertEquals("observed type is wrong", ColumnType.Integer, tabToColumnInfoMode.typeFromValue("93"));
        assertEquals("observed type is wrong", ColumnType.Float, tabToColumnInfoMode.typeFromValue("-2e2"));
        assertEquals("observed type is wrong", ColumnType.Float, tabToColumnInfoMode.typeFromValue("-72.3e-02"));
        assertEquals("observed type is wrong", ColumnType.Float, tabToColumnInfoMode.typeFromValue("+72.3e-2"));
        assertEquals("observed type is wrong", ColumnType.Float, tabToColumnInfoMode.typeFromValue("-3.2"));
        assertEquals("observed type is wrong", ColumnType.Float, tabToColumnInfoMode.typeFromValue("3.2"));
        assertEquals("observed type is wrong", ColumnType.Integer, tabToColumnInfoMode.typeFromValue("4"));
        assertEquals("observed type is wrong", ColumnType.Integer, tabToColumnInfoMode.typeFromValue("-32"));
        // Values are trimmed.
        assertEquals("observed type is wrong", ColumnType.Integer, tabToColumnInfoMode.typeFromValue("  -32  "));
        assertEquals("observed type is wrong", ColumnType.Float, tabToColumnInfoMode.typeFromValue("NaN"));
        assertEquals("observed type is wrong", ColumnType.Float, tabToColumnInfoMode.typeFromValue("-Inf"));
        assertEquals("observed type is wrong", ColumnType.String, tabToColumnInfoMode.typeFromValue("woot"));
        assertEquals("observed type is wrong", ColumnType.Float, tabToColumnInfoMode.typeFromValue("+Infinity"));
        // Integers (parsed from String) cannot start with + but Doubles can. Strange but true.
        assertEquals("observed type is wrong", ColumnType.Float, tabToColumnInfoMode.typeFromValue("+32"));
        assertEquals("observed type is wrong", ColumnType.Float, tabToColumnInfoMode.typeFromValue("+INFINITY"));
        assertEquals("observed type is wrong", ColumnType.Float, tabToColumnInfoMode.typeFromValue("-INF"));
        assertEquals("observed type is wrong", ColumnType.String, tabToColumnInfoMode.typeFromValue("no"));
        assertEquals("observed type is wrong", ColumnType.Integer, tabToColumnInfoMode.typeFromValue("9"));
        assertEquals("observed type is wrong", ColumnType.Float, tabToColumnInfoMode.typeFromValue("nan"));
        assertEquals("observed type is wrong", ColumnType.Float, tabToColumnInfoMode.typeFromValue("+inf"));
        assertEquals("observed type is wrong", ColumnType.Float, tabToColumnInfoMode.typeFromValue("inf"));
        assertEquals("observed type is wrong", ColumnType.Unknown, tabToColumnInfoMode.typeFromValue("    "));
        assertEquals("observed type is wrong", ColumnType.Unknown, tabToColumnInfoMode.typeFromValue(""));
        assertEquals("observed type is wrong", ColumnType.Unknown, tabToColumnInfoMode.typeFromValue(null));
    }
}
