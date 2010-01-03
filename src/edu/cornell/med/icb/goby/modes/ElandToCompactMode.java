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

package edu.cornell.med.icb.goby.modes;

import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.alignments.AlignmentWriter;
import edu.cornell.med.icb.alignments.Alignments;
import edu.cornell.med.icb.goby.maq.ElandMapReader;
import edu.cornell.med.icb.goby.maq.MaqMapEntry;
import edu.cornell.med.icb.goby.maq.MaqMapHeader;
import edu.cornell.med.icb.goby.maq.MaqMapReader;
import edu.cornell.med.icb.reads.Reads;
import edu.cornell.med.icb.reads.ReadsReader;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Properties;

/**
 * Convert an Eland file to the alignment compact format.
 *
 * @author Fabien Campagne
 * @author Kevin Dorff
 */
public class ElandToCompactMode extends AbstractGobyMode {

    /**
     * The mode name.
     */
    public static final String MODE_NAME = "eland-to-compact";
    public static final String MODE_DESCRIPTION = "Convert an Eland file to the alignment compact format."; 


    /**
     * The input file.
     */
    private String inputFile;

    /**
     * The output file.
     */
    private String outputFile;

    /**
     * Max read length for reading / writing MAQ Map files.
     */
    private int maqMapMaxReadLen;
    private String queryReadIdsFilename;
    private String targetReferenceIdsFilename;

    public String getModeName() {
        return MODE_NAME;
    }

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

        inputFile = jsapResult.getString("input");
        outputFile = jsapResult.getString("output");
        queryReadIdsFilename = jsapResult.getString("query-reads-ids");
        targetReferenceIdsFilename = jsapResult.getString("target-reference-ids");

        maqMapMaxReadLen = jsapResult.getInt("maq-map-max-read-len");
        return this;
    }

    /**
     * Run the eland-to-compact mode.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        final Properties props = new Properties();

        final AlignmentWriter writer = new AlignmentWriter(outputFile);

        // first write reference ids to compact header, if these ids are provided on the command line:
        if (targetReferenceIdsFilename != null) {
            final ObjectArrayList<String> ids = processIds(targetReferenceIdsFilename);
            if (ids.size() > 0) {
                writer.setTargetIdentifiersArray(ids.toArray(new String[ids.size()]));
            }
        }
        // write query/read ids to compact header, if provided as well:
        if (queryReadIdsFilename != null) {
            final ObjectArrayList<String> ids = processIds(queryReadIdsFilename);
            if (ids.size() > 0) {
                writer.setQueryIdentifiersArray(ids.toArray(new String[ids.size()]));
            }
        }


        // now convert entries from eland format to compact format:
        final MaqMapReader elandMapReader = new ElandMapReader(
                inputFile, maqMapMaxReadLen, null);
        // turn off indexing of read names. This tool expects readnames that code for read index as integers.
        elandMapReader.setReadNameIdentifiers(null);
        elandMapReader.configure(props);


        final MaqMapHeader header = elandMapReader.getMaqMapHeader();

        int numRead = 0;
        for (final MaqMapEntry entry : elandMapReader) {
            numRead++;
            final Alignments.AlignmentEntry.Builder currentEntry = writer.getAlignmentEntry();
            final float score = -entry.getNumMisMatches();
            final int targetPosition = (int) entry.getActualPosition();
            final int queryIndex = Integer.parseInt(entry.getName().toString());
            final int targetIndex = (int) entry.getReferenceSequenceId();
            final boolean reverseStrand = entry.isMatchingReverseStrand();

            currentEntry.setScore(score);
            currentEntry.setPosition(targetPosition);
            currentEntry.setQueryIndex(queryIndex);
            currentEntry.setTargetIndex(targetIndex);
            currentEntry.setMatchingReverseStrand(reverseStrand);
            writer.appendEntry();
        }
        elandMapReader.close();
        writer.close();
        System.out.println("numRead=" + numRead + " maq header =" + header.toString());
    }

    private ObjectArrayList<String> processIds(final String idsFilename) throws FileNotFoundException {
        final ObjectArrayList<String> ids = new ObjectArrayList<String>();
        final ReadsReader idsReader = new ReadsReader(new FileInputStream(idsFilename));
        boolean atLeastOneId = false;
        while (idsReader.hasNext()) {
            final Reads.ReadEntry readEntry = idsReader.next();
            if (readEntry.hasReadIdentifier()) {
                // resize as necessary:
                ids.size(readEntry.getReadIndex());
                // set element:
                ids.set(readEntry.getReadIndex(), readEntry.getReadIdentifier());
                atLeastOneId = true;
            }
        }
        return ids;
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
        new ElandToCompactMode().configure(args).execute();
    }
}
