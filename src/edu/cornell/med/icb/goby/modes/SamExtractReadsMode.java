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

package edu.cornell.med.icb.goby.modes;

import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.reads.ReadsWriter;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.logging.ProgressLogger;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * Display the sequence variations found in alignments.
 *
 * @author Fabien Campagne
 */
public class SamExtractReadsMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "sam-extract-reads";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Extract reads from a SAM file.";

    /**
     * The input filename.
     */
    private String inputFilename;

    /**
     * The output file.
     */
    private String outputFilename;

    private final int minimumUniqueReadIndices = 1;
    private static final Logger LOG = Logger.getLogger(SamExtractReadsMode.class);


    @Override
    public String getModeName() {
        return MODE_NAME;
    }


    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }

    /**
     * Configure.
     *
     * @param args command line arguments
     * @return this object for chaining
     * @throws java.io.IOException error parsing
     * @throws com.martiansoftware.jsap.JSAPException
     *                             error parsing
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args) throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        inputFilename = jsapResult.getString("input");
        outputFilename = jsapResult.getString("output");
        return this;
    }

    /**
     * Display sequence variations.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {

        final ReadsWriter writer = new ReadsWriter(new FileOutputStream(outputFilename));
        try {
            final ProgressLogger progress = new ProgressLogger(LOG);
            final SAMFileReader parser = new SAMFileReader(new File(inputFilename));
            parser.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

            progress.start();

            final CloseableIterator<SAMRecord> recordCloseableIterator = parser.iterator();

            while (recordCloseableIterator.hasNext()) {
                final SAMRecord samRecord = recordCloseableIterator.next();
                final String readId = samRecord.getReadName();

                writer.setIdentifier(readId);
                writer.setSequence(byteToString(samRecord.getReadBases()));
                // How are quality scores encoded in a SAM file?
                writer.setQualityScores(remove33(samRecord.getReadBases()));
                writer.appendEntry();
                progress.lightUpdate();
            }
        }


        finally {

            writer.close();
        }
    }

    private byte[] remove33(final byte[] input) {
        for (int i = 0; i < input.length; i++) {
            input[i] -= 33;
        }
        return input;
    }


    private CharSequence byteToString(final byte[] input) {
        final MutableString buffer = new MutableString();
        for (int i = 0; i < buffer.length(); i++) {
            buffer.setCharAt(i, (char) input[i]);

        }
        return buffer;
    }

    /**
     * Main method.
     *
     * @param args command line args.
     * @throws com.martiansoftware.jsap.JSAPException
     *                             error parsing
     * @throws java.io.IOException error parsing or executing.
     */
    public static void main(final String[] args) throws JSAPException, IOException {
        new SamExtractReadsMode().configure(args).execute();
    }
}
