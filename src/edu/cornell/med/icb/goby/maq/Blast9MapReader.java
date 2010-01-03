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

import edu.cornell.med.icb.io.TsvToFromMap;
import it.unimi.dsi.lang.MutableString;

import java.io.DataInput;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Map;
import java.util.Properties;

/**
 * Read an Eland MAP file with MAQ data structures.
 * TODO: NOT DONE!!
 * @author Kevin Dorff
 */
public class Blast9MapReader extends MaqMapReader {

    public static final double DEFAULT_MIN_PERCENT_IDENTITY = 95.0;

    public static final int DEFAULT_MIN_ALIGNMENT_LENGTH = 150;

    private double minPercentIdentity;

    private double minAlignmentLength;

    /**
     * Helper object to read eland into a map of data.
     */
    private TsvToFromMap blast9ToMap = new TsvToFromMap();

    /**
     * Constructor.
     * @param filenameVal the filename to read the eland map data from
     * @param maxReadLenVal the max read length value (64 or 128, generally 128).
     * @param header the previous MaqMapHeader to append to or null to start with a new header
     * @throws java.io.IOException error opening file
     */
    public Blast9MapReader(final String filenameVal, final int maxReadLenVal,
                          final MaqMapHeader header) throws IOException {
        this(new DataInputStream(new FileInputStream(filenameVal)), maxReadLenVal, header);
        minPercentIdentity = DEFAULT_MIN_PERCENT_IDENTITY;
        minAlignmentLength = DEFAULT_MIN_ALIGNMENT_LENGTH;
    }

    /**
     * Configure.
     * @param props properties to configure with
     */
    @Override
    public void configure(final Properties props) {
        if (props == null) {
            return;
        }
        if (props.contains("min-percent-identity")) {
            minPercentIdentity = (Double) props.get("min-percent-identity");
        }
        if (props.contains("min-alignment-length")) {
            minAlignmentLength = (Integer) props.get("min-alignment-length");
        }

    }

    /**
     * Constructor.
     * @param dataInput the data input for reading the eland map data. This should be
     * a RandomAccessFile as eland needs to be able to seek() to reset the file pointer
     * during the call to readHeader().
     * @param maxReadLenVal the max read length value (64 or 128, generally 128).
     * @param header the previous MaqMapHeader to append to or null to start with a new header
     * @throws java.io.IOException error reading
     */
    public Blast9MapReader(final DataInput dataInput, final int maxReadLenVal,
                          final MaqMapHeader header) throws IOException {
        super(maxReadLenVal, dataInput, header);
        blast9ToMap.setLenientColumnCount(true);
        blast9ToMap.addColumn("query-id");  //keep
        blast9ToMap.addColumn("subject-id"); //keep
        blast9ToMap.addColumn("percent-identity"); // filter if >= 95.0
        blast9ToMap.addColumn("alignment-length"); // filter if >= 150
        blast9ToMap.addColumn("mismatches"); // keep as # of mismatches
        blast9ToMap.addColumn("gap-openings");
        blast9ToMap.addColumn("q-start");
        blast9ToMap.addColumn("q-end");
        blast9ToMap.addColumn("s-start"); // keep start position
        blast9ToMap.addColumn("s-end");
        blast9ToMap.addColumn("e-value");
        blast9ToMap.addColumn("bit-score");
    }

    /**
     * Read the header information. For an Eland file, there is nothing to do.
     * @return the header
     * @throws java.io.IOException error reading
     */
    @Override
    protected MaqMapHeader readHeader() throws IOException {
        return new MaqMapHeader();
    }

    /**
     * Read the next mapped Eland entry (marked with Ux). This code is based on
     * MAQs eland2maq.cc from MAQ 0.71, the method eland2maq_core(...). The MAQ
     * code splits on whitepsace, this splits on TAB (which is a better method).
     * @return the next Eland entry
     * @throws java.io.IOException error reading
     */
    @Override
    protected MaqMapEntry readEntry() throws IOException {
        // Read the next Ux entry, ignoring all other entries
        Map<String, String> blast9DataMap;
        while (true) {
            final String line = is.readLine();
            if (line == null) {
                throw new EOFException();
            }
            if (line.startsWith("#")) {
                continue;
            }
            blast9DataMap = blast9ToMap.readDataToMap(line);
            final double percentIdentity = textToDouble(blast9DataMap.get("percent-identity"));
            if (percentIdentity < minPercentIdentity) {
                continue;
            }
            final double alignmentLength = textToInt(blast9DataMap.get("alignment-length"));
            if (alignmentLength < minAlignmentLength) {
                continue;
            }
            break;
        }
        final MaqMapEntry entry = new MaqMapEntry(maxReadLen);
        entry.setReadName(new MutableString(blast9DataMap.get("query-id")));
        final int nMis = textToInt(blast9DataMap.get("mismatches"));
        final int refNameIndex = getMaqMapHeader().addRefName(blast9DataMap.get("subject-id"));
        entry.setSeqId(refNameIndex);

        // We don't have a sequence, so we just leave it as 0's.
        // and set size to be the mininum of actual sequence size
        // and maxReadLen. Also, C and quality values are just 0.
        // info1 and info2 encode number of misses just as is done
        // with eland.
        final short[] seqArray = new short[maxReadLen];
        entry.setSeq(seqArray);

        final long seqStart = textToInt(blast9DataMap.get("s-start"));
        final long seqEnd = textToInt(blast9DataMap.get("s-end"));
        final boolean forward;
        final short seqSize;
        // We cannot store more than maxReadLen (128) but we can store
        // that we have a seq size of up to 255. No biggie since with Blast9
        // we aren't storing the sequence in the MAQ Map file, anyway
        if (seqStart > seqEnd) {
            forward = false;
            seqSize = (short) Math.min(seqStart - seqEnd, 255);
        } else {
            forward = true;
            seqSize = (short) Math.min(seqEnd - seqStart, 255);
        }
        entry.setSize(seqSize);
        final short[] cArray = new short[2];
        entry.setC(cArray);
        entry.setFlag((short) 0);
        entry.setDist(0);
        entry.setPos((seqStart - 1) << 1 | (forward ? 0 : 1));

        entry.setInfo1((short) (nMis << 4 | nMis));
        entry.setInfo2((short) (nMis * MaqConstants.DEFAULT_QUAL));
        entry.setMapQual((short) 0);
        entry.setAltQual((short) 0);
        return entry;
    }

    private static double textToDouble(final String text) {
        return Double.parseDouble(text);
    }

    private static int textToInt(final String text) {
        return Integer.parseInt(text);
    }

}
