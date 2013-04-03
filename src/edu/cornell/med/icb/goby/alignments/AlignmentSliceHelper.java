package edu.cornell.med.icb.goby.alignments;

import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

/**
 * A class to help parse command line arguments that represent a genomic slice.
 * @author Fabien Campagne
 *         Date: 4/2/13
 *         Time: 11:57 AM
 */
public class AlignmentSliceHelper {
    /**
     * The string representation of the start offset. Either an long type offset, or chr,position format.
     */
    private String startOffsetArgument;
    private String endOffsetArgument;
    /**
     * Used to log informational and debug messages.
     */
    private static final Log LOG = LogFactory.getLog(AlignmentSliceHelper.class);
    GenomicRange genomicRange=new GenomicRange();

    private long startOffset;
    private long endOffset;


    /*

     */

    /**
     * Parse the string of reference sequences to process during the iteration. The JSAP
     * argument include-reference-names must be defined.
     *
     * @param jsapResult The jsapResult available to the mode.
     * @param inputFilenames Input files for alignments associated with this slice.
     */
    public void parseIncludeReferenceArgument(final JSAPResult jsapResult, String[] inputFilenames) {

        if (!genomicRange.setTargetIds(inputFilenames)) {
           // System.err.println("Unable to obtain target reference identifiers from alignment files.");
            throw new IllegalArgumentException("Unable to obtain target reference identifiers from alignment files.");
        }
        startOffsetArgument = jsapResult.getString("start-position");
        endOffsetArgument = jsapResult.getString("end-position");
        if (startOffsetArgument != null && endOffsetArgument == null ||
                endOffsetArgument != null && startOffsetArgument == null) {
            System.err.println("Start (-s) and end offset (-e) arguments must be specified together or not at all.");
            System.exit(1);
        }
        if (StringUtils.isEmpty(startOffsetArgument) && StringUtils.isEmpty(endOffsetArgument)) {
            // OK: no slice
            startOffset = 0;
            endOffset = Long.MAX_VALUE;
            genomicRange=null;
        } else {
            try {
                startOffset = Long.parseLong(startOffsetArgument);
            } catch (NumberFormatException e) {
                // Unable to parse as long offset
                if (!hasDelimiter(startOffsetArgument)) {

                    System.err.println("start offset must contain a coma or semi-colon delimiter.");
                    System.exit(1);
                }
                genomicRange.startChromosome = getRefId(startOffsetArgument);
                genomicRange.startReferenceIndex = getRefIndex(startOffsetArgument, genomicRange.getTargetIds());
                genomicRange.startPosition = getPosition(startOffsetArgument);
            }
            try {
                endOffset = Long.parseLong(endOffsetArgument);

            } catch (NumberFormatException e) {
                // Unable to parse as long offset
                if (!hasDelimiter(endOffsetArgument)) {

                    System.err.println("end offset must contain a coma or semi-colon delimiter.");
                    System.exit(1);
                }
                genomicRange.endChromosome = getRefId(endOffsetArgument);
                genomicRange.endReferenceIndex = getRefIndex(endOffsetArgument, genomicRange.getTargetIds());
                genomicRange.endPosition = getPosition(endOffsetArgument);
            }


        }
    }

    private boolean hasDelimiter(String startOffsetArgument) {
        if (startOffsetArgument.contains(",")) return true;
        if (startOffsetArgument.contains(":")) return true;
        return false;
    }

    private String getRefId(String offsetArgument) {
        final String[] tokens = offsetArgument.split("[,]");

        return tokens[0];
    }

    private int getRefIndex(String offsetArgument, DoubleIndexedIdentifier referenceIds) {

        final String[] tokens = offsetArgument.split("[,:]");

        int referenceIndex = referenceIds.getIndex(tokens[0]);

        if (referenceIndex == -1) {
            String message = String.format("One of the reference identifier specified for start and end limits does not exist %s. ",
                    tokens[0]);
            LOG.error(message);
            throw new IllegalArgumentException(message);
        }
        return referenceIndex;
    }


    private int getPosition(String offsetArgument) {

        final String[] tokens = offsetArgument.split("[,:]");

        int position = Integer.parseInt(tokens[1]);

        return position;
    }

    public GenomicRange getGenomicRange() {
        return genomicRange;
    }
}
