/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                         Weill Medical College of Cornell University
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
import edu.cornell.med.icb.goby.alignments.AlignedSequence;
import edu.cornell.med.icb.goby.alignments.AlignmentStats;
import edu.cornell.med.icb.goby.alignments.AlignmentTooManyHitsWriter;
import edu.cornell.med.icb.goby.alignments.AlignmentWriter;
import edu.cornell.med.icb.goby.alignments.filters.AlignmentQualityFilter;
import edu.cornell.med.icb.goby.alignments.filters.PercentMismatchesQualityFilter;
import edu.cornell.med.icb.goby.reads.ReadSet;
import edu.cornell.med.icb.goby.reads.Reads;
import edu.cornell.med.icb.goby.reads.ReadsReader;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

/**
 * JSAP files must include baseclass options
 * <p/>
 * "input"
 * "output"
 * "query-reads-ids"
 * "target-reference-ids"
 * "propagate-query-ids"
 * "propagate-target-ids"
 * "read-index-filter"
 * "ambiguity-threshold"
 * "quality-filter-parameters"
 * <p/>
 * <p/>
 * User: stu
 * Date: Sep 11, 2009
 * Time: 9:30:16 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class AbstractAlignmentToCompactMode extends AbstractGobyMode {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(AbstractAlignmentToCompactMode.class);

    /**
     * default ambiguity threshold.
     */
    protected static final int DEFAULT_M_PARAM = 2;

    /**
     * Input file.
     */
    protected String inputFile;

    /**
     * Output file.
     */
    protected String outputFile;

    /**
     * Query / target identifiers.
     */
    protected String queryReadIdsFilename;
    protected String targetReferenceIdsFilename;
    protected boolean propagateQueryIds;
    protected boolean propagateTargetIds;

    /**
     * Conversion parameters.
     */
    protected String qualityFilterParameters = "";
    protected AlignmentQualityFilter qualityFilter;
    protected File readIndexFilterFile;
    protected int mParameter = DEFAULT_M_PARAM;
    protected int numberOfReads;

    /**
     * Scan.
     *
     * @param readIndexFilter
     * @param writer
     * @param targetIds
     * @param tmhWriter
     * @return number of alignment entries written
     * @throws java.io.IOException error parsing
     */
    protected abstract int scan(ReadSet readIndexFilter, IndexedIdentifier targetIds, AlignmentWriter writer,
                                AlignmentTooManyHitsWriter tmhWriter) throws IOException;


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
        //
        final JSAPResult jsapResult = parseJsapArguments(args);

        inputFile = jsapResult.getString("input");
        outputFile = jsapResult.getString("output");

        queryReadIdsFilename = jsapResult.getString("query-reads-ids");
        targetReferenceIdsFilename = jsapResult.getString("target-reference-ids");
        propagateQueryIds = jsapResult.getBoolean("propagate-query-ids");
        propagateTargetIds = jsapResult.getBoolean("propagate-target-ids");

        readIndexFilterFile = jsapResult.getFile("read-index-filter");
        mParameter = jsapResult.getInt("ambiguity-threshold");
        qualityFilterParameters = jsapResult.getString("quality-filter-parameters");

        return this;
    }

    /**
     * Run the last-to-compact mode.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {

        // read target/query identifier lookup table, and initialize output alignment file with this information
        final TransferIds transferIds = new TransferIds().invoke();
        final ReadSet readIndexFilter = transferIds.getReadIndexFilter();
        final AlignmentWriter writer = transferIds.getWriter();
        final IndexedIdentifier targetIds = transferIds.getTargetIds();

        // initialize too-many-hits output file
        final AlignmentTooManyHitsWriter tmhWriter = new AlignmentTooManyHitsWriter(outputFile, mParameter);

        try {

            // initialize numberOfReads
            if (transferIds.numberOfReads != 0) {
                numberOfReads = transferIds.numberOfReads;
            }
            if (numberOfReads <= 0) {
                System.err.println("Cannot determine number of reads. Must set property or provide reads file with -q");
                return;
            }


            // initialize quality filter
            qualityFilter = new PercentMismatchesQualityFilter();
            qualityFilter.setParameters(qualityFilterParameters);

            final int numAligns = scan(readIndexFilter, targetIds, writer, tmhWriter);
            System.out.println("Number of alignments written: " + numAligns);
            if (queryIds.size() > 0 && propagateQueryIds) {
                // we collected query ids, let's write them to the header:
                writer.setQueryIdentifiers(queryIds);
            }
        } finally {
            writer.close();
            tmhWriter.close();
        }
    }

    protected void evaluateStatistics(final AlignedSequence reference, final AlignedSequence query, final AlignmentStats stats) {
        final MutableString queryAligned = query.alignment;
        final MutableString targetAligned = reference.alignment;
        //  assert reference.alignedLength == query.alignedLength :"aligned length differ for queryIndex="+query.sequenceIdentifier;
        final int length = Math.max(query.alignedLength, reference.alignedLength);
        int numIndels = 0;
        int numMismatches = 0;
        for (int i = 0; i < length; i++) {
            final char queryBase = queryAligned.charAt(i);
            final char targetBase = targetAligned.charAt(i);
            if (queryBase == '-' && targetBase != '-' || queryBase != '-' && targetBase == '-') {
                numIndels++;
            }
            if (queryBase != '-' && targetBase != '-' && queryBase != targetBase) {
                numMismatches++;
            }
        }
        stats.numberOfIndels = numIndels;
        stats.numberOfMismatches = numMismatches;
    }


    public void setInputFile(final String inputFile) {
        this.inputFile = inputFile;
    }

    public void setOutputFile(final String outputFile) {
        this.outputFile = outputFile;
    }

    public String getOutputFile() {
        return outputFile;
    }

    public void setQueryReadIdsFilename(final String queryReadIdsFilename) {
        this.queryReadIdsFilename = queryReadIdsFilename;
    }

    public void setTargetReferenceIdsFilename(final String targetReferenceIdsFilename) {
        this.targetReferenceIdsFilename = targetReferenceIdsFilename;
    }

    public void setPropagateTargetIds(final boolean propagateTargetIds) {
        this.propagateTargetIds = propagateTargetIds;
    }

    public void setPropagateQueryIds(final boolean propagateQueryIds) {
        this.propagateQueryIds = propagateQueryIds;
    }

    public void setAmbiguityThreshold(final int mParameter) {
        this.mParameter = mParameter;
    }

    public void setQualityFilterParameters(final String qualityFilterParameters) {
        this.qualityFilterParameters = qualityFilterParameters;
    }

    public void setNumberOfReads(final int numberOfReads) {
        this.numberOfReads = numberOfReads;
    }


    public class TransferIds {
        private ReadSet readIndexFilter;
        private AlignmentWriter writer;
        private IndexedIdentifier targetIds;
        public int numberOfReads;
        private int numberOfTargets;

        public ReadSet getReadIndexFilter() {
            return readIndexFilter;
        }

        public AlignmentWriter getWriter() {
            return writer;
        }

        public IndexedIdentifier getTargetIds() {
            return targetIds;
        }

        /**
         * Postcondition: output Ids.size() == maximumSequenceIndex + 1.
         */
        private ObjectArrayList<String> processIds(final String idsFilename) throws FileNotFoundException {
            final ObjectArrayList<String> ids = new ObjectArrayList<String>(500000);
            int maxSequenceIndex = -1;
            ReadsReader idsReader = null;
            try {
                idsReader = new ReadsReader(new FileInputStream(idsFilename));

                boolean atLeastOneId = false;
                ids.size(500000);
                while (idsReader.hasNext()) {
                    final Reads.ReadEntry readEntry = idsReader.next();
                    final int readIndex = readEntry.getReadIndex();
                    if (readEntry.hasReadIdentifier()) {
                        // resize as necessary:
                        if (readIndex >= ids.size()) {
                            final double newSize = ids.size() * 1.5;
                            //  System.out.println("resizing to " + newSize);
                            ids.size((int) newSize);
                        }
                        // set element:
                        ids.set(readIndex, readEntry.getReadIdentifier());
                        atLeastOneId = true;
                    }
                    maxSequenceIndex = Math.max(maxSequenceIndex, readIndex);
                }
            } finally {
                if (idsReader != null) {
                    try {
                        idsReader.close();
                    } catch (IOException e) {
                        LOG.warn("Error closing " + idsFilename, e);
                    }
                }
            }
            ids.size(maxSequenceIndex + 1);
            ids.trim();
            assert ids.size() == maxSequenceIndex + 1;
            return ids;
        }

        public AbstractAlignmentToCompactMode.TransferIds invoke() throws IOException {
            // setup multiplicity set:
            readIndexFilter = new ReadSet();
            if (readIndexFilterFile == null) {
                readIndexFilter = null;
            } else {
                readIndexFilter.load(readIndexFilterFile);
            }

            writer = new AlignmentWriter(outputFile);
            targetIds = new IndexedIdentifier();

            // first write reference ids to compact header, if these ids are provided on the command line:
            if (targetReferenceIdsFilename != null) {
                System.out.println("Scanning target file..");
                // read reference ids from file
                final ObjectArrayList<String> ids = processIds(targetReferenceIdsFilename);
                this.numberOfTargets = ids.size();
                System.out.println("Target file had " + this.numberOfTargets + " entries.");
                // write ids to header
                writer.setNumTargets(this.numberOfTargets);
                if (this.numberOfTargets > 0 && propagateTargetIds) {
                    writer.setTargetIdentifiersArray(ids.toArray(new String[ids.size()]));
                    System.out.println("Wrote " + ids.size() + " target ids to alignment header.");
                } else {
                    System.out.println("Target ids are not propagated to output header.");
                }

                for (int index = 0; index < ids.size(); index++) {
                    final String id = ids.get(index);
                    if (id != null) {
                        targetIds.put(new MutableString(id), index);
                    }
                }
            }
            // write query/read ids to compact header, if provided as well:
            if (queryReadIdsFilename != null) {
                // read read ids from file
                System.out.println("Scanning query file..");
                final ObjectArrayList<String> ids = processIds(queryReadIdsFilename);
                this.numberOfReads = ids.size();
                System.out.println("Query file had " + this.numberOfReads + " entries.");
                for (int queryIndex = 0; queryIndex < ids.size(); ++queryIndex) {
                    final String id = ids.get(queryIndex);
                    if (id != null) {

                        queryIds.registerIdentifier(new MutableString(id));
                    }
                }
                // write ids to header
                writer.setNumQueries(this.numberOfReads);
                if (this.numberOfReads > 0 && propagateQueryIds) {
                    writer.setQueryIdentifiersArray(ids.toArray(new String[ids.size()]));
                    System.out.println("Wrote " + ids.size() + " query ids to alignment header.");
                } else {
                    System.out.println("Query ids are not propagated to output header.");
                }
            }
            return this;
        }
    }

    final IndexedIdentifier queryIds = new IndexedIdentifier();
}
