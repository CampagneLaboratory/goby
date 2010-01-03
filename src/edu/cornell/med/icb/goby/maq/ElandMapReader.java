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

import java.io.DataInput;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Map;

/**
 * Read an Eland MAP file with MAQ data structures.
 * TODO: NOT DONE!!
 * @author Kevin Dorff
 */
public class ElandMapReader extends MaqMapReader {

    /**
     * Temporary sequence array for when "r-f" column is "R".
     */
    private final short[] tempSeqArray;

    /**
     * Helper object to read eland into a map of data.
     */
    private TsvToFromMap elandToMap = new TsvToFromMap();

    /**
     * Constructor.
     * @param filenameVal the filename to read the eland map data from
     * @param maxReadLenVal the max read length value (64 or 128, generally 128).
     * @param header the previous MaqMapHeader to append to or null to start with a new header
     * @throws IOException error opening file
     */
    public ElandMapReader(final String filenameVal, final int maxReadLenVal,
                          final MaqMapHeader header) throws IOException {
        this(new DataInputStream(new FileInputStream(filenameVal)), maxReadLenVal, header);
    }

    /**
     * Constructor.
     * @param dataInput the data input for reading the eland map data. This should be
     * a RandomAccessFile as eland needs to be able to seek() to reset the file pointer
     * during the call to readHeader().
     * @param maxReadLenVal the max read length value (64 or 128, generally 128).
     * @param header the previous MaqMapHeader to append to or null to start with a new header
     * @throws IOException error reading
     */
    public ElandMapReader(final DataInput dataInput, final int maxReadLenVal,
                          final MaqMapHeader header) throws IOException {
        super(maxReadLenVal, dataInput, header);
        elandToMap.setLenientColumnCount(true);
        elandToMap.addColumn("name");
        elandToMap.addColumn("sequence");
        elandToMap.addColumn("match-info");
        elandToMap.addColumn("c0"); // Number of exact matches
        elandToMap.addColumn("c1"); // Number of 1-error matches
        elandToMap.addColumn("c2"); // Number of 2-error matches
        elandToMap.addColumn("ref-name");
        elandToMap.addColumn("pos");
        elandToMap.addColumn("r-f"); // Forward / Reverse
        tempSeqArray = new short[maxReadLen];
    }

    /**
     * Read the header information. For an Eland file, there is nothing to do.
     * @return the header
     * @throws IOException error reading
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
     * @throws IOException error reading
     */
    @Override
    protected MaqMapEntry readEntry() throws IOException {
        // Read the next Ux entry, ignoring all other entries
        Map<String, String> elandDataMap;
        while (true) {
            final String line = is.readLine();
            if (line == null) {
                throw new EOFException();
            }
            elandDataMap = elandToMap.readDataToMap(line);
            if (elandDataMap.get("match-info").charAt(0) == 'U') {
                break;
            }
        }
        final MaqMapEntry entry = new MaqMapEntry(maxReadLen);
        entry.setReadName(elandDataMap.get("name"));
        final short[] seqArray = new short[maxReadLen];
        entry.setSeq(seqArray);
        final int nMis = elandDataMap.get("match-info").charAt(1) - '0';
        final int refNameIndex = getMaqMapHeader().addRefName(
                fixElandName(elandDataMap.get("ref-name")));
        entry.setSeqId(refNameIndex);
        final String sequence = elandDataMap.get("sequence");
        final short seqSize = (short) sequence.length();
        entry.setSize(seqSize);
        for (int i = 0; i < seqSize; i++) {
            final short seqNum = MaqConstants.NST_NT4_TABLE[sequence.charAt(i)];
            if (seqNum > 3) {
                seqArray[i] = 0;
            } else {
                // Does this really make sense?
                seqArray[i] = (short) (seqNum << 6 | MaqConstants.DEFAULT_QUAL);
            }
        }
        if (elandDataMap.get("r-f").equals("R")) {
            for (int i = seqSize - 1; i >= 0; i--) {
                final short tempValue;
                if (seqArray[i] == 0) {
                    tempValue = 0;
                } else {
                    tempValue = (short) ((0xc0 - (seqArray[i] & 0xc0)) | (seqArray[i] & 0x3f));
                }
                tempSeqArray[seqSize - i - 1] = tempValue;
            }
            System.arraycopy(tempSeqArray, 0, seqArray, 0, seqSize);
        }
        final short[] cArray = new short[2];
        entry.setC(cArray);
        cArray[0] = (short) textToMaxInt(elandDataMap.get("c0"), 255);
        cArray[1] = (short) textToMaxInt(elandDataMap.get("c1"), 255);
        entry.setFlag((short) 0);
        entry.setDist(0);
        final long rawPos = Long.parseLong(elandDataMap.get("pos"));
        entry.setPos((rawPos - 1) << 1 | (elandDataMap.get("r-f").equals("F") ? 0 : 1));

        entry.setInfo1((short) (nMis << 4 | nMis));
        entry.setInfo2((short) (nMis * MaqConstants.DEFAULT_QUAL));
        final short mapQual = calculateMapQuality(
                cArray[0], cArray[1], (short) textToMaxInt(elandDataMap.get("c2"), 255));
        seqArray[maxReadLen - 1] = mapQual;
        entry.setMapQual(mapQual);
        entry.setAltQual(mapQual);
        return entry;
    }

    private static int textToMaxInt(final String text, final int maxInt) {
        final int value = Integer.parseInt(text);
        if (value > maxInt) {
            return maxInt;
        } else {
            return value;
        }
    }

    /**
     * Calculate map quality.
     * @param c0 c0 value from eland file
     * @param c1 c1 value from eland file
     * @param c2 c2 value from eland file
     * @return the map quality value
     */
    private static short calculateMapQuality(final short c0, final short c1, final short c2) {
            if (c0 == 1) {
                    if (c1 == 0 && c2 == 0) {
                        return (short) (3 * MaqConstants.DEFAULT_QUAL);
                    }
                    if (c1 == 0) {
                        return (short) (2 * MaqConstants.DEFAULT_QUAL - MaqConstants.LOG_N[c2]);
                    }
                    return (short) (MaqConstants.DEFAULT_QUAL - MaqConstants.LOG_N[c1]);
            }
            if (c1 == 1) {
                    if (c2 == 0) {
                        return (short) (2 * MaqConstants.DEFAULT_QUAL);
                    }
                    return (short) (MaqConstants.DEFAULT_QUAL - MaqConstants.LOG_N[c2]);
            }
            if (c2 == 1) {
                return (short) (MaqConstants.DEFAULT_QUAL - 3);
            }
            return MaqConstants.DEFAULT_QUAL;
    }

    /**
     * If the eland name is SOMESTUFF/ENSmorestuff just return the
     * ENSmorestuff part.
     * @param brokenRefName the identifier
     * @return the fixed identifer
     */
    private String fixElandName(final String brokenRefName) {
        if (brokenRefName == null) {
            return brokenRefName;
        }
        final int slashIndex = brokenRefName.indexOf("/ENS");
        if (slashIndex == -1) {
            return brokenRefName;
        } else {
            return brokenRefName.substring(slashIndex + 1);
        }
    }

}
