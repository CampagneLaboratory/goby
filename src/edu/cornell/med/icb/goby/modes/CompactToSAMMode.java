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

package edu.cornell.med.icb.goby.modes;

import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.alignments.*;
import edu.cornell.med.icb.goby.reads.DualRandomAccessSequenceCache;
import edu.cornell.med.icb.goby.reads.QualityEncoding;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import edu.cornell.med.icb.util.VersionUtils;
import it.unimi.dsi.Util;
import it.unimi.dsi.fastutil.ints.Int2ObjectAVLTreeMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.logging.ProgressLogger;
import net.sf.samtools.*;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * Export a Goby alignment to the BAM format.
 *
 * @author Kevin Dorff
 */
public class CompactToSAMMode extends AbstractGobyMode {

    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(CompactToSAMMode.class);

    private static final DateFormat GOBY_DATE_FORMAT = new SimpleDateFormat("dd:MMM:yyyy");

    /**
     * The mode name.
     */
    private static final String MODE_NAME = "compact-to-sam";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Exports a compact alignment to the SAM/BAM format." +
            "This tool tries to export a Goby alignment comprehensively to the BAM format. Since Goby 2.0.";

    /**
     * For logging progress.
     */
    ProgressLogger progress;

    /**
     * BAM/SAM output alignment filename.
     */
    private String output;

    /**
     * Compact-alignment input basename (or one of the files of the input, such as "/path/something.entries").
     */
    private String inputBasename;

    /**
     * The input genome in either 'compact random-access-genome' format or 'fa.gz + fa.gz.fai' format.
     * The random-access-genome file can be made from the .fa.gz reference using the build-sequence-cache mode.
     * If using the random-access-genome input, specify any one of the files in the random-access-genome.
     * If using the 'fa.gz + fa.gz.fai' input, specify the 'fa.gz' file but make sure the '.fa.gz.fai' file
     * is located in the same directory.
     * Using the random-access-genome format will be considerably faster.
     */
    private String inputGenome;

    /**
     * An "unset" value for startPosition and endPosition.
     */
    private boolean hasStartOrEndPosition;

    /**
     * The start position for the reformat.
     */
    private long startPosition;

    /**
     * The end position for the reformat.
     */
    private long endPosition = Long.MAX_VALUE;

    private QualityEncoding qualityEncoding = QualityEncoding.PHRED;

    private CompactToSAMIterateAlignments alignmentIterator;

    private SAMFileWriter outputSam;

    private SAMFileHeader samHeader;

    private SAMRecordFactory samRecordFactory;

    private RandomAccessSequenceInterface genome;

    private ExportableAlignmentEntryData exportData;

    private Int2ObjectMap<Int2ObjectMap<ExportableAlignmentEntryData>> queryIndexToFragmentsMap;

    /**
     * Flag to indicate if log4j was configured.
     */
    private boolean debug;

    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }

    /**
     * Get the quality encoding scale used for the input fastq file.
     *
     * @return the quality encoding scale used for the input fastq file
     */
    public QualityEncoding getQualityEncoding() {
        return qualityEncoding;
    }

    /**
     * Set the quality encoding scale to be used for the input fastq file.
     * Acceptable values are "Illumina", "Sanger", and "Solexa".
     *
     * @param qualityEncoding the quality encoding scale to be used for the input fastq file
     */
    public void setQualityEncoding(final QualityEncoding qualityEncoding) {
        this.qualityEncoding = qualityEncoding;
    }

    public String getOutput() {
        return output;
    }

    public void setOutput(final String output) {
        this.output = output;
    }

    public String getInputBasename() {
        return inputBasename;
    }

    public void setInputBasename(final String inputBasename) {
        this.inputBasename = inputBasename;
    }

    public String getInputGenome() {
        return inputGenome;
    }

    public void setInputGenome(final String inputGenome) {
        this.inputGenome = inputGenome;
    }

    public long getEndPosition() {
        return endPosition;
    }

    public void setEndPosition(final long endPosition) {
        this.endPosition = endPosition;
    }

    public long getStartPosition() {
        return startPosition;
    }

    public void setStartPosition(final long startPosition) {
        this.startPosition = startPosition;
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

        qualityEncoding = QualityEncoding.valueOf(jsapResult.getString("quality-encoding").toUpperCase());
        // don't even dare go through the debugging code if log4j was not configured. The debug code
        // is way too slow to run unintentionally in production!

        inputBasename = jsapResult.getString("input-basename");
        inputGenome = jsapResult.getString("input-genome");
        output = jsapResult.getString("output");
        alignmentIterator = new CompactToSAMIterateAlignments();
        alignmentIterator.parseIncludeReferenceArgument(jsapResult);
        genome = new DualRandomAccessSequenceCache();
        try {
            ((DualRandomAccessSequenceCache)genome).load(inputGenome);
        } catch (ClassNotFoundException e) {
            throw new IOException("Could not load genome", e);
        }

        if (jsapResult.contains("start-position") || jsapResult.contains("end-position")) {
            hasStartOrEndPosition = true;
            startPosition = jsapResult.getLong("start-position", 0L);
            endPosition = jsapResult.getLong("end-position", Long.MAX_VALUE);
            if (startPosition < 0) {
                startPosition = 0;
            }
            if (endPosition < 0) {
                endPosition = Long.MAX_VALUE;
            }
            if (startPosition == 0 && endPosition == 0) {
                endPosition = Long.MAX_VALUE;
            }
        }

        if (startPosition > endPosition) {
            throw new JSAPException("Start position must not be greater than the end position");
        }
        return this;
    }

    @Override
    public void execute() throws IOException {
        debug = Util.log4JIsConfigured();
        queryIndexToFragmentsMap = new Int2ObjectAVLTreeMap<Int2ObjectMap<ExportableAlignmentEntryData>>();

        final AlignmentReader gobyReader = new AlignmentReaderImpl(inputBasename);
        gobyReader.readHeader();
        final DoubleIndexedIdentifier targetIdentifiers = new DoubleIndexedIdentifier(gobyReader.getTargetIdentifiers());
        gobyReader.close();

        exportData = new ExportableAlignmentEntryData(genome, qualityEncoding, targetIdentifiers);
        final String[] basenames = new String[1];
        basenames[0] = inputBasename;
        progress = new ProgressLogger(LOG);
        progress.displayFreeMemory = true;
        progress.start();
        if (alignmentIterator == null) {
            alignmentIterator = new CompactToSAMIterateAlignments();
        }

        final int seekTargetIndex = -1;
        final String seekTargetName = "";
        final int seekStartPosition = -1;

        // For debugging set to a targetName, targetIndex, and startPosition
        /*
        final int seekTargetIndex = 10;
        final String seekTargetName = "11";
        final int seekStartPosition = 67675226;
        */

        if (hasStartOrEndPosition) {
            alignmentIterator.iterate(new FileSlice(startPosition, endPosition), basenames);
        } else {
            if (seekStartPosition > 0 && seekTargetIndex > 0) {
                final GenomicRange range = new GenomicRange();
                range.startReferenceIndex = seekTargetIndex;
                range.startPosition = seekStartPosition;
                range.startChromosome = seekTargetName;
                range.endReferenceIndex = seekTargetIndex;
                range.endPosition = seekStartPosition + 1;
                range.endChromosome = seekTargetName;
                alignmentIterator.iterate(range, basenames);
            } else {
                alignmentIterator.iterate(basenames);
            }
        }
        if (outputSam != null) {
            outputSam.close();
        }
        progress.stop();
    }

    private class CompactToSAMIterateAlignments extends IterateAlignments {

        private long numWritten;
        private ReadOriginInfo readOriginInfo;
        private boolean hasReadGroups;
        private final List<ExportableAlignmentEntryData> completeSpliceFragments;
        final LinkedList<Integer> needFragmentIndexes;
        final LinkedList<Integer> foundFragmentIndexes;

        private CompactToSAMIterateAlignments() {
            completeSpliceFragments = new ObjectArrayList<ExportableAlignmentEntryData>();
            needFragmentIndexes = new LinkedList<Integer>();
            foundFragmentIndexes = new LinkedList<Integer>();
        }

        private void initializeSam(final AlignmentReader alignmentReader) {
            // Because splices cannot be written in a sorted manner, we can never consider the output to be sorted.
            final boolean outputIsSorted = false;

            // Gather the target identifiers, supply them to the SAM file
            samHeader = new SAMFileHeader();
            final SAMSequenceDictionary samTargetDictionary = new SAMSequenceDictionary();
            final IndexedIdentifier gobyTargetIdentifiers = alignmentReader.getTargetIdentifiers();
            final DoubleIndexedIdentifier gobyBackTargetIdentifiers =
                    new DoubleIndexedIdentifier(gobyTargetIdentifiers);
            for (int i = 0; i < gobyTargetIdentifiers.size(); i++) {
                final String gobyTargetName = gobyBackTargetIdentifiers.getId(i).toString();
                final int gobyTargetLength = alignmentReader.getTargetLength()[i];
                final SAMSequenceRecord samSequenceRecord = new SAMSequenceRecord(gobyTargetName, gobyTargetLength);
                samTargetDictionary.addSequence(samSequenceRecord);
            }
            exportReadGroups(alignmentReader);
            samHeader.setSequenceDictionary(samTargetDictionary);
            final SAMProgramRecord gobyVersionProgRec = new SAMProgramRecord("Goby");
            gobyVersionProgRec.setProgramVersion(VersionUtils.getImplementationVersion(GobyDriver.class));
            samHeader.addProgramRecord(gobyVersionProgRec);
            outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(samHeader, outputIsSorted, new File(output));
            samRecordFactory = new DefaultSAMRecordFactory();
        }

        private void exportReadGroups(final AlignmentReader alignmentReader) {
            readOriginInfo = alignmentReader.getReadOriginInfo();
            if (readOriginInfo.size() > 0) {
                hasReadGroups = true;
                // Goby alignment has read origin information, export as BAM read groups:
                for (final Alignments.ReadOriginInfo roi : readOriginInfo.getPbList()) {
                    final SAMReadGroupRecord readGroup = new SAMReadGroupRecord(roi.getOriginId());
                    if (roi.hasSample()) {
                        readGroup.setSample(roi.getSample());
                    }
                    if (roi.hasPlatform()) {
                        readGroup.setSample(roi.getPlatform());
                    }
                    if (roi.hasPlatformUnit()) {
                        readGroup.setSample(roi.getPlatformUnit());
                    }
                    if (roi.hasLibrary()) {
                        readGroup.setSample(roi.getLibrary());
                    }
                    if (roi.hasRunDate()) {
                        final String runDate = roi.getRunDate();
                        try {

                            readGroup.setRunDate(GOBY_DATE_FORMAT.parse(runDate));
                        } catch (ParseException e) {
                            LOG.error("Unable to parse Goby date: " + runDate + " ignoring runDate read origin.");
                        }
                    }
                    samHeader.addReadGroup(readGroup);
                }
                exportData.setReadGroupInfo(readOriginInfo);
            }
        }

        @Override
        public void processAlignmentEntry(final AlignmentReader alignmentReader,
                                          final Alignments.AlignmentEntry alignmentEntry) {

            if (outputSam == null) {
                initializeSam(alignmentReader);
            }

            exportData.buildFrom(alignmentEntry);

            completeSpliceFragments.clear();
            Int2ObjectMap<ExportableAlignmentEntryData> fragIndexToAlignmentsMap;
            boolean splicedFragment = false;
            if (alignmentEntry.hasSplicedForwardAlignmentLink() ||
                    alignmentEntry.hasSplicedBackwardAlignmentLink()) {
                splicedFragment = true;
                final int queryIndex = alignmentEntry.getQueryIndex();
                fragIndexToAlignmentsMap = queryIndexToFragmentsMap.get(queryIndex);
                if (fragIndexToAlignmentsMap == null) {
                    fragIndexToAlignmentsMap = new Int2ObjectAVLTreeMap<ExportableAlignmentEntryData>();
                    queryIndexToFragmentsMap.put(queryIndex, fragIndexToAlignmentsMap);
                }
                fragIndexToAlignmentsMap.put(exportData.getAlignmentEntry().getFragmentIndex(),
                        ExportableAlignmentEntryData.duplicateFrom(exportData));
                findCompleteSpliceFragments(fragIndexToAlignmentsMap, queryIndexToFragmentsMap, queryIndex);

            }

            if (splicedFragment) {
                if (!completeSpliceFragments.isEmpty()) {
                    outputSplicedFragments(completeSpliceFragments);
                }
            } else {
                // This is suitable for single alignment entries or paired end
                outputSingle(exportData);
            }

        }

        /**
         * Given a list of all of the fragments for a single query index, return a list of
         * fragments that represent a complete splice formation. If the alignment is paired end,
         * this will represent a single end of the pair, not both of the pairs.
         * If none of the splices are complete, this will return null.
         * If not null, this will return a NEW LIST, not the same fragIndexToAlignmentsMap list,
         * additionally, the elements in the returned list will be removed from fragIndexToAlignmentsMap.
         * NOTE: for performance completeSpliceFragments is used repeatedly.
         *
         * @param fragIndexToAlignmentsMap the list of found fragments for a desired query index.
         * @param queryIndexToFragmentsMap
         * @param queryIndex
         */
        private void findCompleteSpliceFragments(
                final Int2ObjectMap<ExportableAlignmentEntryData> fragIndexToAlignmentsMap, Int2ObjectMap<Int2ObjectMap<ExportableAlignmentEntryData>> queryIndexToFragmentsMap, int queryIndex) {
            for (final Map.Entry<Integer, ExportableAlignmentEntryData> entry :
                    fragIndexToAlignmentsMap.entrySet()) {
                // For THIS entry, see if all related forward/backward fragments exists in in  fragIndexToAlignmentsMap
                // We need to do this for every entry in our map because this map COULD contain segments
                // that span multiple alignments (such as if we are aligning ambiguity > 1 or with pairs).
                needFragmentIndexes.clear();
                foundFragmentIndexes.clear();

                final Alignments.AlignmentEntry alignmentEntry = entry.getValue().getAlignmentEntry();
                needFragmentIndexes.add(alignmentEntry.getFragmentIndex());
                foundFragmentIndexes.add(alignmentEntry.getFragmentIndex());
                if (alignmentEntry.hasSplicedForwardAlignmentLink()) {
                    walkFragments(fragIndexToAlignmentsMap,
                            alignmentEntry.getSplicedForwardAlignmentLink(), true);
                }
                if (alignmentEntry.hasSplicedBackwardAlignmentLink()) {
                    walkFragments(fragIndexToAlignmentsMap,
                            alignmentEntry.getSplicedBackwardAlignmentLink(), false);
                }
                // We've now walked the entire available distance, forward and backward. See if the fragment
                // indexes that we need could be found.
                if (needFragmentIndexes.equals(foundFragmentIndexes)) {
                    // We have a complete set of fragments
                    for (final int fragIndex : needFragmentIndexes) {
                        completeSpliceFragments.add(fragIndexToAlignmentsMap.get(fragIndex));
                        fragIndexToAlignmentsMap.remove(fragIndex);
                    }
                    // since we have a complete list, we can now forget the queryIndex and its fragments:
                    queryIndexToFragmentsMap.remove(queryIndex);
                    break;
                }
            }
        }

        /**
         * Used to walk, forward or backward, the splice fragments. Helps with determining if we have
         * all the fragments for a single spliced alignment entry.
         *
         * @param fragIndexToAlignmentsMap the map of fragment index to alignment entry
         * @param requiredRelated          the related splice fragment we are walking from
         * @param walkForward              true if we are walking forward
         */
        private void walkFragments(final Int2ObjectMap<ExportableAlignmentEntryData> fragIndexToAlignmentsMap,
                                   final Alignments.RelatedAlignmentEntry requiredRelated, final boolean walkForward) {
            final int fragmentIndex = requiredRelated.getFragmentIndex();
            if (walkForward) {
                // Build needFragmentIndexes, foundFragmentIndexes in the order of the aligned splices
                needFragmentIndexes.addLast(fragmentIndex);
            } else {
                needFragmentIndexes.addFirst(fragmentIndex);
            }

            final ExportableAlignmentEntryData entry = fragIndexToAlignmentsMap.get(fragmentIndex);
            if (entry != null) {
                // The current needed entry was found. Look forward or backward for more required fragments
                if (walkForward) {
                    foundFragmentIndexes.addLast(fragmentIndex);
                } else {
                    foundFragmentIndexes.addFirst(fragmentIndex);
                }
                final Alignments.AlignmentEntry alignmentEntry = entry.getAlignmentEntry();
                if (walkForward && alignmentEntry.hasSplicedForwardAlignmentLink()) {
                    walkFragments(fragIndexToAlignmentsMap,
                            alignmentEntry.getSplicedForwardAlignmentLink(), true);
                }
                if (!walkForward && alignmentEntry.hasSplicedBackwardAlignmentLink()) {
                    walkFragments(fragIndexToAlignmentsMap,
                            alignmentEntry.getSplicedBackwardAlignmentLink(), false);
                }
            }
        }

        private void outputSplicedFragments(final List<ExportableAlignmentEntryData> spliceFragments) {
            final ExportableAlignmentEntryData exportData =
                    ExportableAlignmentEntryData.mergeSpliceFragments(spliceFragments);
            outputSingle(exportData);
        }

        /**
         * Output a single alignment segment. If paired end, this will output of the pair halfs.
         * Spliced alignments will be written with outputMultiFragments once all the fragments of
         * that alignment have bee encountered.
         *
         * @param toExport the alignment entry to output
         */
        private void outputSingle(final ExportableAlignmentEntryData toExport) {
            if (toExport.isInvalid()) {
                LOG.warn(toExport.toString());
                return;
            }
            if (debug) {
                LOG.debug("Wrote qi=" + toExport.getQueryIndex());
            }
            final SAMRecord samRecord = samRecordFactory.createSAMRecord(samHeader);

            samRecord.setReferenceIndex(toExport.getTargetIndex());
            samRecord.setFlags(toExport.getPairFlags());
            samRecord.setReadName(toExport.getReadName());
            samRecord.setAlignmentStart(toExport.getStartPosition());
            samRecord.setMappingQuality(toExport.getMappingQuality());

            samRecord.setReadString(toExport.getReadBasesOriginal());
            samRecord.setBaseQualities(toExport.getReadQualities().toByteArray());
            samRecord.setCigarString(toExport.getCigarString());
            samRecord.setAttribute("MD", toExport.getMismatchString());
            for (String bamAttribute : toExport.getBamAttributesList()) {
                bamAttribute = bamAttribute.replaceAll("[\n\r]", "");
                final String[] tokens = bamAttribute.split(":");
                samRecord.setAttribute(tokens[0], getValue(tokens));
                if (debug) {
                    LOG.debug(String.format("Writing %s:%s", tokens[0], getValue(tokens)));
                }
            }
            if (toExport.hasMate()) {
                samRecord.setMateReferenceIndex(toExport.getMateReferenceIndex());
                samRecord.setMateAlignmentStart(toExport.getMateAlignmentStart());
                samRecord.setInferredInsertSize(toExport.getInferredInsertSize());
            }
            if (hasReadGroups) {
                samRecord.setAttribute("RG", toExport.getReadGroup());
            }
            try {
                outputSam.addAlignment(samRecord);
            } catch (RuntimeException e) {
                System.out.println("entry: \n" + toExport.getAlignmentEntry().toString());
                System.out.println(toExport.toString());
                throw e;
            }
            numWritten++;
            progress.lightUpdate();
        }
    }

    private Object getValue(final String[] tokens) {
        final String type = tokens[1];
        if ("Z".equals(type)) {
            return tokens[2];
        }
        if ("i".equals(type)) {
            return Integer.parseInt(tokens[2]);
        }
        if ("A".equals(type)) {
            return  tokens[2].charAt(0);
        }
        LOG.warn("Attribute type %c is currently not supported, storing as string type");
        return tokens[2];
    }

    private String getTag(final String bamAttribute) {
        return bamAttribute.split(":")[0];
    }

    public void setGenome(RandomAccessSequenceInterface genome) {
        this.genome=genome;
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
        final CompactToSAMMode processor = new CompactToSAMMode();
        processor.configure(args);
        processor.execute();
    }


}
