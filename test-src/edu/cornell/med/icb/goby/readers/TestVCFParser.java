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

package edu.cornell.med.icb.goby.readers;

import edu.cornell.med.icb.goby.readers.vcf.ColumnInfo;
import edu.cornell.med.icb.goby.readers.vcf.ColumnType;
import edu.cornell.med.icb.goby.readers.vcf.Columns;
import edu.cornell.med.icb.goby.readers.vcf.VCFParser;
import org.junit.Test;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import static org.junit.Assert.*;

/**
 * Test the VCF parser
 *
 * @author Fabien Campagne
 *         Date: Mar 26, 2011
 *         Time: 3:15:13 PM
 */
public class TestVCFParser {


    @Test
    public void testParse1() throws FileNotFoundException, VCFParser.SyntaxException {
        VCFParser parser = new VCFParser(new FileReader("test-data/vcf/example1.vcf"));
        parser.readHeader();
        Columns cols = parser.getColumns();
        ColumnInfo info = cols.find("INFO");
        assertTrue(info.hasField("DP"));
        assertTrue(info.hasField("MQ"));
        assertTrue(info.hasField("FQ"));
        ColumnInfo format = cols.find("FORMAT");
        assertTrue(format.hasField("GT"));
        assertTrue(format.hasField("GQ"));
        assertTrue(format.hasField("GL"));
        assertTrue(format.hasField("SP"));
        assertEquals("Likelihoods for RR,RA,AA genotypes (R=ref,A=alt)",
                format.getField("GL").description);
    }

    /**
     * Verify we can handle both TSV headers that contain spaces AND we correctly identify
     * the column types.
     * @throws IOException
     * @throws VCFParser.SyntaxException
     */
    @Test
    public void tsvHeaderHasSpacesAndDetectTypes() throws IOException, VCFParser.SyntaxException {
        final VCFParser parser = new VCFParser("test-data/vcf/tsv-with-header-spaces.tsv");
        parser.setCacheTsvColumnTypes(false);
        parser.readHeader();
        assertEquals("Incorrect number of columns", 48, parser.countAllFields());
        assertEquals("Incorrect number of columns", 48, parser.getNumberOfColumns());
        for (int i = 0; i < parser.getNumberOfColumns(); i++) {
            if (i <= 1) {
                assertEquals("Incorrect column type", ColumnType.String, parser.getColumnType(i));
            } else {
                assertEquals("Incorrect column type", ColumnType.Float, parser.getColumnType(i));
            }
        }
        parser.close();
    }

    @Test
    public void testParse1FilenameGzipped() throws IOException, VCFParser.SyntaxException {
        // same as Parse1 with gzipped filename
        VCFParser parser = new VCFParser("test-data/vcf/example1.vcf.gz");
        parser.readHeader();
        Columns cols = parser.getColumns();
        ColumnInfo info = cols.find("INFO");
        assertTrue(info.hasField("DP"));
        assertTrue(info.hasField("MQ"));
        assertTrue(info.hasField("FQ"));
        ColumnInfo format = cols.find("FORMAT");
        assertTrue(format.hasField("GT"));
        assertTrue(format.hasField("GQ"));
        assertTrue(format.hasField("GL"));
        assertTrue(format.hasField("SP"));
        assertEquals("Likelihoods for RR,RA,AA genotypes (R=ref,A=alt)",
                format.getField("GL").description);
    }

    @Test
    public void testParse1WithGroups() throws FileNotFoundException, VCFParser.SyntaxException {
        VCFParser parser = new VCFParser(new FileReader("test-data/vcf/example-with-groups.vcf"));
        parser.readHeader();
        Columns cols = parser.getColumns();
        ColumnInfo info = cols.find("INFO");
        assertTrue(info.hasField("DP"));
        assertTrue(info.getField("DP").group.equals("DEPTHS"));
        assertTrue(info.hasField("MQ"));
        assertTrue(info.getField("MQ").group.equals("RMS"));
    }

    @Test
    public void testParseExample2() throws FileNotFoundException, VCFParser.SyntaxException {
        VCFParser parser = new VCFParser(new FileReader("test-data/vcf/example2.tsv"));
        parser.readHeader();
        Columns cols = parser.getColumns();

        assertNotNull("Fixed CHROM column must exist", cols.find("CHROM"));
        assertNotNull("Fixed POS column must exist", cols.find("POS"));
        assertNotNull("Fixed ID column must exist", cols.find("ID"));
        assertNotNull("Fixed ALT column must exist", cols.find("ALT"));
        assertNotNull("Fixed QUAL column must exist", cols.find("QUAL"));
        assertNotNull("Fixed QUAL column must exist", cols.find("INFO"));
        assertNotNull("optional sample column 1 must exist", cols.find("results/IPBKRNW/IPBKRNW-replicate.bam"));
        assertNotNull("optional sample column w must exist", cols.find("results/IPBKRNW/IPBKRNW-sorted.bam"));

        assertTrue("file must have one line at least", parser.hasNextDataLine());
        assertEquals("0", parser.getColumnValue(0).toString());
        assertEquals("145497099", parser.getColumnValue(1).toString());
        assertEquals("1/1:25,3,0:42", parser.getColumnValue(10).toString());
    }


    @Test
    public void testParseExample3() throws FileNotFoundException, VCFParser.SyntaxException {
        VCFParser parser = new VCFParser(new FileReader("test-data/vcf/example3.tsv"));
        parser.readHeader();
        Columns cols = parser.getColumns();

        assertNotNull("Fixed CHROM column must exist", cols.find("CHROM"));
        assertNotNull("Fixed POS column must exist", cols.find("POS"));
        assertNotNull("Fixed ID column must exist", cols.find("ID"));
        assertNotNull("Fixed ALT column must exist", cols.find("ALT"));
        assertNotNull("Fixed QUAL column must exist", cols.find("QUAL"));
        assertNotNull("Fixed QUAL column must exist", cols.find("INFO"));
        assertNotNull("optional sample column 1 must exist", cols.find("results/IPBKRNW/IPBKRNW-replicate.bam"));
        assertNotNull("optional sample column w must exist", cols.find("results/IPBKRNW/IPBKRNW-sorted.bam"));

    }

    @Test
    public void testParseExampleBug() throws IOException, VCFParser.SyntaxException {
        VCFParser parser = new VCFParser("test-data/stats-writer/GT-bug-example-vcf");
        parser.readHeader();
        Columns cols = parser.getColumns();


        int formatGlobalFieldIndex = parser.getGlobalFieldIndex("FORMAT", "GT");

        while (parser.hasNextDataLine()) {

            assertEquals("GT", parser.getStringFieldValue(formatGlobalFieldIndex));
            parser.next();

        }


    }

    @Test
    public void testParseValue2() throws FileNotFoundException, VCFParser.SyntaxException {
        VCFParser parser = new VCFParser(new FileReader("test-data/vcf/example2.tsv"));
        parser.readHeader();
        assertTrue(parser.hasNextDataLine());

        assertEquals("145497099", parser.getStringColumnValue(1));
        assertEquals(".", parser.getStringColumnValue(2));
        assertEquals("A", parser.getStringColumnValue(3));
        parser.next();
        assertTrue(parser.hasNextDataLine());
        assertEquals("29389393", parser.getStringColumnValue(1));
        assertEquals("1", parser.getStringColumnValue(0));
        assertEquals("AC", parser.getStringColumnValue(3));
    }

    @Test
    public void testParseValue3() throws FileNotFoundException, VCFParser.SyntaxException {
        VCFParser parser = new VCFParser(new FileReader("test-data/vcf/example3.tsv"));
        parser.readHeader();
        assertTrue(parser.hasNextDataLine());

        assertEquals("145497099", parser.getStringColumnValue(1));
        assertEquals(".", parser.getStringColumnValue(2));
        assertEquals("A", parser.getStringColumnValue(3));
        parser.next();
        assertTrue(parser.hasNextDataLine());
        assertEquals("29389393", parser.getStringColumnValue(1));
        assertEquals("1", parser.getStringColumnValue(0));
        assertEquals("AC", parser.getStringColumnValue(3));
    }

    public void testHasNext() throws FileNotFoundException, VCFParser.SyntaxException {
        VCFParser parser = new VCFParser(new FileReader("test-data/vcf/example2.tsv"));
        parser.readHeader();
        int lineCounter = 0;
        while (parser.hasNextDataLine()) {
            lineCounter++;
            parser.next();

        }
        assertEquals(3177, lineCounter);


    }

    @Test
    public void testParseFixedColumns() throws FileNotFoundException, VCFParser.SyntaxException {
        VCFParser parser = new VCFParser(new FileReader("test-data/vcf/example1.vcf"));
        parser.readHeader();

        Columns cols = parser.getColumns();
        assertNotNull("Fixed CHROM column must exist", cols.find("CHROM"));
        assertNotNull("Fixed POS column must exist", cols.find("POS"));
        assertNotNull("Fixed ID column must exist", cols.find("ID"));
        assertNotNull("Fixed ALT column must exist", cols.find("ALT"));
        assertNotNull("Fixed QUAL column must exist", cols.find("QUAL"));
        assertNotNull("Fixed QUAL column must exist", cols.find("INFO"));
        assertNotNull("optional sample column 1 must exist", cols.find("results/IPBKRNW/IPBKRNW-replicate.bam"));
        assertNotNull("optional sample column w must exist", cols.find("results/IPBKRNW/IPBKRNW-sorted.bam"));
    }

    @Test
    public void testColumnOrder() throws FileNotFoundException, VCFParser.SyntaxException {
        VCFParser parser = new VCFParser(new FileReader("test-data/vcf/example1.vcf"));
        parser.readHeader();

        Columns cols = parser.getColumns();
        assertEquals(0, cols.find("CHROM").columnIndex);
        assertEquals(1, cols.find("POS").columnIndex);
        assertEquals(2, cols.find("ID").columnIndex);
        assertEquals(3, cols.find("REF").columnIndex);
        assertEquals(4, cols.find("ALT").columnIndex);
        assertEquals(5, cols.find("QUAL").columnIndex);
        assertEquals(6, cols.find("FILTER").columnIndex);
        assertEquals(7, cols.find("INFO").columnIndex);
        assertEquals(8, cols.find("FORMAT").columnIndex);
        assertEquals(9, cols.find("results/IPBKRNW/IPBKRNW-replicate.bam").columnIndex);
        assertEquals(10, cols.find("results/IPBKRNW/IPBKRNW-sorted.bam").columnIndex);
    }

    @Test
    public void testParseFields() throws FileNotFoundException, VCFParser.SyntaxException {
        VCFParser parser = new VCFParser(new FileReader("test-data/vcf/example-small.vcf"));
        parser.readHeader();
        assertTrue(parser.hasNextDataLine());
        for (int i = 0; i < parser.countAllFields(); i++) {
            System.out.printf("field %s gfi:%d value: %s%n", parser.getFieldName(i), i,
                    parser.getStringFieldValue(i));
        }
        assertEquals(25, parser.countAllFields());
        assertEquals("145497099", parser.getStringFieldValue(1));
        assertEquals("0", parser.getStringFieldValue(0));
        assertEquals("A", parser.getStringFieldValue(3));

        assertEquals("1/1", parser.getStringFieldValue(parser.getGlobalFieldIndex("results/IPBKRNW/IPBKRNW-replicate.bam", "GT")));
        assertEquals("015,4,0", parser.getStringFieldValue(parser.getGlobalFieldIndex("results/IPBKRNW/IPBKRNW-replicate.bam", "PL")));
        assertEquals("11", parser.getStringFieldValue(parser.getGlobalFieldIndex("results/IPBKRNW/IPBKRNW-replicate.bam", "GQ")));


        assertEquals("1/1", parser.getStringFieldValue(parser.getGlobalFieldIndex("results/IPBKRNW/IPBKRNW-sorted.bam", "GT")));
        assertEquals("25,3,0", parser.getStringFieldValue(parser.getGlobalFieldIndex("results/IPBKRNW/IPBKRNW-sorted.bam", "PL")));
        assertEquals("42", parser.getStringFieldValue(parser.getGlobalFieldIndex("results/IPBKRNW/IPBKRNW-sorted.bam", "GQ")));
        parser.next();
        assertTrue(parser.hasNextDataLine());

    }

    @Test
    public void testParseVariableOrder() throws FileNotFoundException, VCFParser.SyntaxException {
        VCFParser parser = new VCFParser(new FileReader("test-data/vcf/example-flags.vcf"));
        parser.readHeader();
        int indelFieldIndex = parser.getGlobalFieldIndex("INFO", "INDEL");
        int af1FieldIndex = parser.getGlobalFieldIndex("INFO", "AF1");

        assertTrue(parser.hasNextDataLine());

        for (int i = 0; i < parser.countAllFields(); i++) {
            System.out.printf("field %s value: %s%n", parser.getFieldName(i),
                    parser.getStringFieldValue(i));
        }

        assertEquals("INDEL", parser.getStringFieldValue(indelFieldIndex));
        assertEquals("0.9999", parser.getStringFieldValue(af1FieldIndex));
        parser.next();
        assertTrue(parser.hasNextDataLine());

        assertEquals("", parser.getStringFieldValue(indelFieldIndex));
        parser.next();
        assertTrue(parser.hasNextDataLine());
        assertEquals("INDEL", parser.getStringFieldValue(indelFieldIndex));
    }

    @Test
    public void testParseTricky2() throws IOException, VCFParser.SyntaxException {
        VCFParser parser = new VCFParser("test-data/vcf/tricky2.vcf");
        parser.readHeader();
        while (parser.hasNextDataLine()) {
            for (int i = 0; i < parser.countAllFields(); i++) {
                final String name = parser.getFieldName(i);
                final String stringFieldValue = parser.getStringFieldValue(i);
                /*  System.out.printf("field %s gfi:%d value: %s%n", name, i,
                      stringFieldValue);
                */

            }
            parser.next();
        }
    }



    /*
    @Test
    public void testParseTrickyLarge() throws IOException, VCFParser.SyntaxException {
        //VCFParser parser = new VCFParser("/home/gobyweb/GOBYWEB_RESULTS/campagne/NXMONDD/NXMONDD.vcf.gz");
        VCFParser parser = new VCFParser("test-data/vcf/tricky-large.vcf");
        parser.readHeader();
        while (parser.hasNextDataLine()) {
            for (int i = 0; i < parser.countAllFields(); i++) {
                final String name = parser.getFieldName(i);
                               final String stringFieldValue = parser.getStringFieldValue(i);
                               //System.out.printf("field %s gfi:%d value: %s%n", name, i,
                                 //      stringFieldValue);



            }
            parser.next();
        }
    }*/
    /*
 @Test
    public void testParseTrickyLarge2() throws IOException, VCFParser.SyntaxException {
        //VCFParser parser = new VCFParser("/home/gobyweb/GOBYWEB_RESULTS/campagne/NXMONDD/NXMONDD.vcf.gz");
        VCFParser parser = new VCFParser("test-data/vcf/tricky-large.vcf.gz");
        parser.readHeader();
        while (parser.hasNextDataLine()) {
            for (int i = 0; i < parser.countAllFields(); i++) {
                final String name = parser.getFieldName(i);
                               final String stringFieldValue = parser.getStringFieldValue(i);
                               //System.out.printf("field %s gfi:%d value: %s%n", name, i,
                                 //      stringFieldValue);
               
            }
            parser.next();
        }
    }@Test
    public void testParseTrickyLarge3() throws IOException, VCFParser.SyntaxException {
        //VCFParser parser = new VCFParser("/home/gobyweb/GOBYWEB_RESULTS/campagne/NXMONDD/NXMONDD.vcf.gz");
        VCFParser parser = new VCFParser("test-data/vcf/tricky-large-regizpp.vcf.gz");
        parser.readHeader();
        while (parser.hasNextDataLine()) {
            for (int i = 0; i < parser.countAllFields(); i++) {
                final String name = parser.getFieldName(i);
                               final String stringFieldValue = parser.getStringFieldValue(i);
                               //System.out.printf("field %s gfi:%d value: %s%n", name, i,
                                 //      stringFieldValue);

            }
            parser.next();
        }
    }

    @Test
    public void testParseTrickyLarge4() throws IOException, VCFParser.SyntaxException {
        //VCFParser parser = new VCFParser("/home/gobyweb/GOBYWEB_RESULTS/campagne/NXMONDD/NXMONDD.vcf.gz");
        VCFParser parser = new VCFParser("test-data/vcf/tricky-large-rebgzipp.vcf.gz");
        parser.readHeader();
        while (parser.hasNextDataLine()) {
            for (int i = 0; i < parser.countAllFields(); i++) {
                final String name = parser.getFieldName(i);
                               final String stringFieldValue = parser.getStringFieldValue(i);
                               //System.out.printf("field %s gfi:%d value: %s%n", name, i,
                                 //      stringFieldValue);

            }
            parser.next();
        }
    }       */
}

