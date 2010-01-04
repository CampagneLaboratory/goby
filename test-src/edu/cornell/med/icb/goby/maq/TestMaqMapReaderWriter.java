/*
 * Copyright (C) 2009 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 * WEILL MEDICAL COLLEGE OF CORNELL UNIVERSITY MAKES NO REPRESENTATIONS
 * ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY PURPOSE. IT IS PROVIDED
 * "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY. THE WEILL MEDICAL COLLEGE
 * OF CORNELL UNIVERSITY SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY
 * THE USERS OF THIS SOFTWARE.
 */

package edu.cornell.med.icb.goby.maq;

import edu.cornell.med.icb.iterators.TextFileLineIterator;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.zip.GZIPInputStream;

/**
 * Test reading / writing MaqMap files.
 *
 * @author Kevin Dorff
 */
public class TestMaqMapReaderWriter {
    private static final Log LOG = LogFactory.getLog(TestMaqMapReaderWriter.class);

    /**
     * Make sure the temporary output file exists. If "clean" was run it will be gone.
     * @throws IOException if the output directory cannot be created
     */
    @Before
    public void setUp() throws IOException {
        FileUtils.forceMkdir(new File("test-results/maq-temp-output"));
    }

    /**
     * Read and check MAQ MAP Reader & Writer.
     * @throws java.io.IOException error reading
     */
    @Test
    public void testMaqMapReaderWriter() throws IOException {
        final int maxReadLen = 128;

        final MaqMapReader maqMapReader = new MaqMapReader(
                "test-data/maq/sample-map.map", maxReadLen, null);
        final MaqMapWriter maqMapWriter = new MaqMapWriter(
                "test-results/maq-temp-output/sample-map-rewritten.map",
                maxReadLen, maqMapReader.getMaqMapHeader());

        final MaqMapHeader header = maqMapReader.getMaqMapHeader();
        for (final MaqMapEntry entry : maqMapReader) {
            maqMapWriter.writeMaqMapEntry(entry);
        }
        assertEquals("Incorrect number of records read from source map file",
                82313, header.getNumberOfReads());
        maqMapReader.close();
        maqMapWriter.close();

        assertSameUnGZippedContents(
                "test-data/maq/sample-map.map",
                "test-results/maq-temp-output/sample-map-rewritten.map");

        // Verify the text file output of the re-written map
        testMaqMapReaderToText();
    }

    /**
     * This isn't it's own test. It is called by testMaqMapReaderWriter() as it depends
     * on the output of that test.
     * @throws IOException error reading / writing
     */
    public void testMaqMapReaderToText() throws IOException {
        final int maxReadLen = 128;
        MaqMapReader maqMapReader = null;
        PrintWriter out = null;
        try {
            maqMapReader = new MaqMapReader(
                    "test-results/maq-temp-output/sample-map-rewritten.map", maxReadLen, null);
            out = new PrintWriter("test-results/maq-temp-output/sample-map-short.txt");
            final MaqMapHeader header = maqMapReader.getMaqMapHeader();
            out.println(header.toString());
            int numRead = 0;
            while (numRead++ < 100) {
                if (!maqMapReader.hasNext()) {
                    break;
                }
                out.println(maqMapReader.next().format(header));
            }
            assertEquals(101, numRead);

        } finally {
            IOUtils.closeQuietly(out);
            if (maqMapReader != null) {
                maqMapReader.close();
            }
        }

        assertSameTextContents(
                "test-data/maq/sample-map-short-expected.txt",
                "test-results/maq-temp-output/sample-map-short.txt");
    }

    /**
     * Test the eland reader.
     * @throws IOException error reading
     */
    @Test
    public void testElandReader() throws IOException {
        final int maxReadLen = 128;

        final MaqMapReader elandMapReader = new ElandMapReader(
                "test-data/eland/eland-std-sample.eland", maxReadLen, null);
        final MaqMapWriter maqMapWriter = new MaqMapWriter(
                "test-results/maq-temp-output/eland-std-sample.map",
                maxReadLen, elandMapReader.getMaqMapHeader());

        final MaqMapHeader header = elandMapReader.getMaqMapHeader();
        int numRead = 0;
        for (final MaqMapEntry entry : elandMapReader) {
            numRead++;
            maqMapWriter.writeMaqMapEntry(entry);
        }
        elandMapReader.close();
        maqMapWriter.close();
        assertSameUnGZippedContents(
                "test-data/eland/eland-std-sample-expected.map",
                "test-results/maq-temp-output/eland-std-sample.map");
        LOG.debug("numRead=" + numRead + " maq header =" + header.toString());
    }

    /**
     * Take two files that contain contents that are GZIP'd and verify that the
     * underlying data is the same (when it is uncompressed).
     * @param expectedName the filname that contains the expected contents
     * @param actualName the filname that contains the actual contents
     * @throws IOException error reading
     */
    private void assertSameUnGZippedContents(final String expectedName, final String actualName)
            throws IOException {
        InputStream expectedStream = null;
        InputStream actualStream = null;
        try {
            expectedStream = new GZIPInputStream(new FileInputStream(expectedName));
            actualStream = new GZIPInputStream(new FileInputStream(actualName));

            assertSameInputStreamContents(expectedStream,  actualStream);
        } finally {
            IOUtils.closeQuietly(expectedStream);
            IOUtils.closeQuietly(actualStream);
        }
    }

    /**
     * The advantage of this over FileUtils.contentEquals() is this knows it is dealing
     * with a text file and it doesn't matter if the line separates are exactly the same.
     * @param expectedName the filname that contains the expected contents
     * @param actualName the filname that contains the actual contents
     * @throws IOException error reading
     */
    private void assertSameTextContents(final String expectedName, final String actualName)
            throws IOException {
        InputStream expectedIs = null;
        InputStream actualIs = null;
        try {
            expectedIs = new FileInputStream(expectedName);
            actualIs = new FileInputStream(actualName);
            final Iterator<String> expected = new TextFileLineIterator(expectedIs).iterator();
            final Iterator<String> actual = new TextFileLineIterator(actualIs).iterator();
            int lineNo = 0;
            while (true) {
                assertEquals("hasNext return incorrect", expected.hasNext(), actual.hasNext());
                if (!expected.hasNext()) {
                    break;
                }
                assertEquals("Line number " + lineNo + " (0 based) was incorrect",
                        expected.next(), actual.next());
                lineNo++;
            }
        } finally {
            IOUtils.closeQuietly(expectedIs);
            IOUtils.closeQuietly(actualIs);
        }
    }

    /**
     * Given two input streams, verify their binary contents are identical.
     * @param expectedStream expected stream
     * @param actualStream actual stream
     * @throws IOException error reading
     */
    private void assertSameInputStreamContents(
            final InputStream expectedStream, final InputStream actualStream) throws IOException {
        final int bufferSize = 4 * 1024 * 1024;
        final byte[] expectedBuffer = new byte[bufferSize];
        final byte[] actualBuffer = new byte[bufferSize];
        int pos = 0;
        while (true) {
            final int expectedReadCount = smarterReadFully(expectedStream, expectedBuffer);
            final int actualReadCount = smarterReadFully(actualStream,  actualBuffer);
            assertEquals("Read count differs", expectedReadCount, actualReadCount);
            if (expectedReadCount == -1) {
                break;
            }
            for (int i = 0; i < expectedReadCount; i++) {
                if (expectedBuffer[i] != actualBuffer[i]) {
                    fail("Inputstreams differ at position = " + pos);
                }
                pos++;
            }
        }
    }

    /**
     * Attempt to read from in to fill all of b. The buffer b will always be filled
     * UNLESS end of file is reached. This will return the number of bytes read which
     * should be equal to b.length unless the end of file is reached.
     * This will return -1 when the end of file is hit (and no bytes were read).
     * This method is a modification of DataInputStream.readFully(...).
     * @param in the InputStream to read from.
     * @param b the buffer to fill
     * @return the number of bytes read or -1 of no more bytes are available.
     * @throws IOException error reading
     */
    private int smarterReadFully(final InputStream in, final byte[] b) throws IOException {
        final int off = 0;
        final int len = b.length;
        int n = 0;
        while (n < len) {
            final int count = in.read(b, off + n, len - n);
            if (count < 0) {
                if (n == 0) {
                    return -1;
                } else {
                    return n;
                }
            }
            n += count;
        }
        return n;
    }

    /**
     * Main method for running the code outside of the normal JUnit framework.
     * @param args command line args
     * @throws IOException error reading / writing
     */
    public static void main(final String[] args) throws IOException {
        final TestMaqMapReaderWriter processor = new TestMaqMapReaderWriter();
        processor.testMaqMapReaderWriter();
    }
}
