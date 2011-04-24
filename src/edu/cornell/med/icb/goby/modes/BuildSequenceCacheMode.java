/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
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
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceCache;
import org.apache.commons.io.IOUtils;
import org.apache.commons.io.FilenameUtils;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

/**
 * Build a random access sequence cache. Can be used to provide random access to the
 * sequence of a mammalian genome.
 *
 * @author Fabien Campagne
 */
public class BuildSequenceCacheMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "build-sequence-cache";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION =
            "Converts a fasta input file to a random access cache.";

    /**
     * The Fasta input file.
     */
    private String inputFile;

    /**
     * The basename to build the cache.
     */

    private String basename;

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
     * @throws IOException error parsing
     * @throws JSAPException error parsing
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args) throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        inputFile = jsapResult.getString("input");
        basename = jsapResult.getString("basename");
        if (basename==null) {
            String filename=inputFile;
            if (filename.endsWith(".gz")) {
                // remove the .gz before removing the next extension
                // this is done so that file.fasta.gz results in basename=file.
                filename=  FilenameUtils.removeExtension(filename);                
            }
            basename= FilenameUtils.removeExtension(filename);
            System.out.printf("Automatically determined basename=%s, please use --basename if you prefer " +
                    "different output filenames.%n", basename);
        }
        return this;
    }

    /**
     * Convert the Fasta input file to a random access cache.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        final RandomAccessSequenceCache cacheBuilder = new RandomAccessSequenceCache();

        InputStream input = null;
        try {
            System.out.println("Loading and compressing genome..");
            if (inputFile.endsWith(".compact-reads")) {
                input = new FileInputStream(inputFile);
                cacheBuilder.loadCompact(input);
            } else {
                if (inputFile.endsWith(".gz")) {
                    input = new GZIPInputStream(new FileInputStream(inputFile));

                } else {
                    input = new FileInputStream(inputFile);
                }
                cacheBuilder.loadFasta(new InputStreamReader(input));
            }

            System.out.println("Done loading input. Starting to write random access cacheBuilder.");
            cacheBuilder.save(basename);
            System.out.println("Compressed genome was written to basename "+basename);
        } finally {
            IOUtils.closeQuietly(input);
        }
        System.out.println("Sequence cacheBuilder written to disk with basename " + basename);
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
        new BuildSequenceCacheMode().configure(args).execute();
    }
}
