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

package edu.cornell.med.icb.goby.stats;

import edu.cornell.med.icb.goby.readers.vcf.ColumnType;
import edu.cornell.med.icb.goby.readers.vcf.GroupAssociations;
import edu.cornell.med.icb.goby.readers.vcf.VCFParser;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.lang.MutableString;
import org.junit.Assert;
import org.junit.Test;

import java.io.StringReader;
import java.io.StringWriter;
import java.util.Arrays;

import static junit.framework.Assert.assertNotNull;
import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

/**
 * @author Fabien Campagne
 *         Date: Apr 7, 2011
 *         Time: 10:14:34 AM
 */
public class TestVCFWriter {
    @Test
    public void testCodeGenotypes() {
        VCFWriter writer = new VCFWriter(new StringWriter());

        assertExpectedGenotype(writer, new String[]{"A", "C"}, 'A', "C", "0/1");
        assertExpectedGenotype(writer, new String[]{"C", "A"}, 'A', "C", "0/1");

        assertExpectedGenotype(writer, new String[]{"C", "A", "T"}, 'A', "C,T", "0/1/2");
        assertExpectedGenotype(writer, new String[]{"C", "A", "T"}, 'C', "A,T", "0/1/2");

    }


    private void assertExpectedGenotype(VCFWriter writer, String[] alleles, char refBase, String altBases, String expectedGenotype) {
        final ObjectArrayList<String> refAlleles = new ObjectArrayList<String>();
        final ObjectArrayList<String> altAlleles = new ObjectArrayList<String>();

        refAlleles.clear();
        refAlleles.add(Character.toString(refBase));
        altAlleles.clear();
        altAlleles.addAll(Arrays.asList(altBases.split(",")));
        final MutableString calculatedGenotype = writer.codeGenotype(alleles, refAlleles, altAlleles);
        assertEquals("Coded genotype must match expected", expectedGenotype, calculatedGenotype.toString());
    }

    @Test
    public void testWriteFieldGroups() {
        final StringWriter stringWriter = new StringWriter();
        VCFWriter writer = new VCFWriter(stringWriter);

        int fieldIndex0 = writer.defineField("INFO",
                "p-value1",
                1, ColumnType.Float,
                "A P-value", /* groups: */"p-value");

        int fieldIndex1 = writer.defineField("INFO",
                "p-value2",
                1, ColumnType.Float,
                "A P-value", /* groups: */"p-value");
        String samples[] = new String[]{"SampleA", "SampleB"};
        writer.defineSamples(samples);
        int zygFieldIndex = writer.defineField("FORMAT", "Zygosity", 1, ColumnType.String, "Zygosity", /* groups: */ "zygozity", "sample-data");
        int secondFormatField = writer.defineField("FORMAT", "Another", 1, ColumnType.String, "Another", /* groups: */ "another", "sample-data");
        writer.setWriteFieldGroupAssociations(true);
        writer.writeHeader();
        writer.close();

        final String textContent = stringWriter.getBuffer().toString();
        System.out.println(textContent);
        assertTrue(textContent.contains("##FieldGroupAssociations="));
        assertTrue(textContent.contains("CHROM=cross-sample-field"));
        assertTrue(textContent.contains("POS=genomic-coordinate"));
        assertTrue(textContent.contains("INFO/p-value1=p-value"));
        assertTrue(textContent.contains("INFO/p-value1=p-value"));
        assertTrue(textContent.contains("INFO/p-value2=p-value"));
        assertTrue(textContent.contains("FORMAT/Zygosity=sample-data"));
        assertTrue(textContent.contains("FORMAT/Zygosity=zygozity"));
        assertTrue(textContent.contains("FORMAT/Another=sample-data"));
        assertTrue(textContent.contains("FORMAT/Another=another"));

        // now check that we can read group associations with VCFParser:
        VCFParser parser = new VCFParser(new StringReader(textContent));
        try {
            parser.readHeader();
            final GroupAssociations associations = parser.getGroupAssociations();

            // check that an INFO/column is associated with its groups:
            Assert.assertTrue(associations.listGroups("p-value1").contains("p-value"));
            Assert.assertTrue(associations.listGroups("INFO/p-value1").contains("p-value"));
            final ObjectArraySet<String> columnsWithGroup = associations.getColumnsWithGroup("p-value");
            assertNotNull(columnsWithGroup);
            Assert.assertTrue(columnsWithGroup.contains("INFO/p-value1"));
            Assert.assertTrue(columnsWithGroup.contains("INFO/p-value2"));
            final ObjectArraySet<String> sampleColumnsWithGroup = associations.getColumnsWithGroup("SampleA");
            Assert.assertTrue(sampleColumnsWithGroup.contains("SampleA[Zygosity]"));
            Assert.assertTrue(sampleColumnsWithGroup.contains("SampleA[Another]"));
            junit.framework.Assert.assertEquals("cross-sample-field,p-value", associations.listGroupsAsString("INFO/p-value1"));
        } catch (VCFParser.SyntaxException e) {
            e.printStackTrace();
            fail("some syntax error was reported: " + e.getMessage());
        }
    }
}
