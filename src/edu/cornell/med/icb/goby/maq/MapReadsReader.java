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
import edu.cornell.med.icb.goby.fastxReaders.FastXEntry;
import edu.cornell.med.icb.goby.fastxReaders.FastXReader;
import it.unimi.dsi.fastutil.ints.Int2ObjectArrayMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import org.apache.commons.lang.StringUtils;

import java.io.EOFException;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.zip.GZIPInputStream;

/**
 * Read an Eland MAP file with MAQ data structures.
 * TODO: NOT DONE!!
 *
 * @author Kevin Dorff
 */
public class MapReadsReader extends MaqMapReader {

    private final Int2ObjectMap<String> indexToChromosomeMap;

    private FastXReader csFastaReader;

    /**
     * Constructor.
     *
     * @param csFastaFile the .csfasta mapping output from mapreads / corona-lite.
     * @param cmapFile used when mapping to create the .csfasta file
     * @param maxReadLenVal the max read length value (64 or 128, generally 128).
     * @param header the previous MaqMapHeader to append to or null to start with a new header
     * @throws java.io.IOException error opening file
     */
    public MapReadsReader(
            final String csFastaFile, final String cmapFile, final int maxReadLenVal,
            final MaqMapHeader header) throws IOException {
        this((csFastaFile.endsWith(".gz") ? new GZIPInputStream(new FileInputStream(csFastaFile)) :
                new FileInputStream(csFastaFile)),
                new FileInputStream(cmapFile), maxReadLenVal, header);
    }

    /**
     * Constructor.
     *
     * @param csFastaInputStream the .csfasta mapping output from mapreads / corona-lite.
     * @param maxReadLenVal the max read length value (64 or 128, generally 128).
     * @param header the previous MaqMapHeader to append to or null to start with a new header
     * @throws java.io.IOException error reading
     */
    public MapReadsReader(final InputStream csFastaInputStream, final InputStream cmapFileInputStream,
                          final int maxReadLenVal, final MaqMapHeader header) throws IOException {
        super(maxReadLenVal, null, header);
        csFastaReader = new FastXReader(csFastaInputStream);
        indexToChromosomeMap = new Int2ObjectArrayMap<String>();
        for (final String line : new TextFileLineIterator(cmapFileInputStream)) {
            final String[] parts = StringUtils.split(line, "\t");
            indexToChromosomeMap.put(Integer.parseInt(parts[0]), parts[1]);
        }
        cmapFileInputStream.close();
    }

    /**
     * Read the header information. For an Eland file, there is nothing to do.
     *
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
     *
     * @return the next Eland entry
     * @throws java.io.IOException error reading or EOFException at the end of file
     */
    @Override
    protected MaqMapEntry readEntry() throws IOException {
        // Read the next csfasta entry
        String[] fields = null;
        FastXEntry csFastaEntry;

        final String matchName;
        final String refName;
        final String sequence;
        String matchType;
        long bestMatchPosition;
        String sign;
        short[] counts;

        while (true) {
            if (!csFastaReader.hasNext()) {
                throw new EOFException();
            }
            csFastaEntry = csFastaReader.next();
            fields = StringUtils.split(csFastaEntry.getEntryHeader().toString(), ",");
            if (fields.length == 1) {
                continue;
            }

            final String[] keys = new String[fields.length - 1];
            final int[] values = new int[fields.length - 1];
            for (int i = 1; i < fields.length; i++) {
                final String[] hitParts = StringUtils.split(fields[i], ".");
                keys[i - 1] = hitParts[0];
                values[i - 1] = Integer.parseInt(hitParts[1]);
            }

            String bestmatch = "_";
            matchType = "NM";
            sign = "F";
            int ref = 0;

            counts = collectCounts(values, 3);

            if (counts[0] == 1) {
                matchType = "U0";
                bestmatch = keys[indexForValue(values, 0)];
            } else if (counts[0] > 1) {
                matchType = "R0";
            } else if (counts[1] == 1) {
                matchType = "U1";
                bestmatch = keys[indexForValue(values, 1)];
            } else if (counts[1] > 1) {
                matchType = "R1";
            } else if (counts[2] == 1) {
                matchType = "U2";
                bestmatch = keys[indexForValue(values, 2)];
            } else if (counts[2] > 1) {
                matchType = "R2";
            }
            if (matchType.charAt(0) != 'U') {
                continue;
            }

            if (bestmatch.indexOf("_") != -1) {
                final String[] bestMatchParts = StringUtils.split(bestmatch, "_");
                ref = Integer.parseInt(bestMatchParts[0]);
                bestmatch = bestMatchParts[1];
            }
            bestMatchPosition = Long.parseLong(bestmatch);
            if (bestMatchPosition < 0) {
                sign = "R";
                bestMatchPosition = Math.abs(bestMatchPosition);
            }

            // Effectively convert the data to Eland and use the same process to convert to MAQ MAP
            // as we used with ElandMapReader
            matchName = fields[0].substring(1);
            sequence = csFastaEntry.getSequence().toString();
            refName = indexToChromosomeMap.get(ref);
            break;
        }

        final MaqMapEntry entry;
        if (mutableReads) {
            entry = mutableEntry;
        } else {
            entry = new MaqMapEntry(maxReadLen, readNameIdentifiers);
        }

        entry.setReadName(matchName);
        final short[] seqArray = new short[maxReadLen];
        entry.setSeq(seqArray);
        final int nMis = matchType.charAt(1) - '0';
        final int refNameIndex = getMaqMapHeader().addRefName(refName);
        entry.setSeqId(refNameIndex);

        // Remove one character from the sequence because it is color space
        final short seqSize = (short) Math.min(sequence.length() - 1, 255);
        entry.setSize(seqSize);

        final short[] cArray = new short[2];
        entry.setC(cArray);
        cArray[0] = counts[0];
        cArray[1] = counts[1];
        entry.setFlag((short) 0);
        entry.setDist(0);
        entry.setPos((bestMatchPosition - 1) << 1 | (sign.equals("F") ? 0 : 1));

        entry.setInfo1((short) (nMis << 4 | nMis));
        entry.setInfo2((short) (nMis * MaqConstants.DEFAULT_QUAL));
        final short mapQual = calculateMapQuality(counts[0], counts[1], counts[2]);
        seqArray[maxReadLen - 1] = mapQual;
        entry.setMapQual(mapQual);
        entry.setAltQual(mapQual);
        return entry;
    }

    /**
     * Search the counts array. Counts the number of '0' values and retuns this
     * count in return[0]... This is done for 0..numToCollect-1.
     * Will not return a value in return[x] > 255.
     *
     * @param counts the input array of counts
     * @return the number of times counts of 0..numToCollect-1 are found in counts
     */
    private short[] collectCounts(final int[] counts, final int numToCollect) {
        final short[] results = new short[numToCollect];
        for (final int count : counts) {
            if (count < numToCollect) {
                if (results[count] < 255) {
                    results[count] = (short) (results[count] + 1);
                }
            }
        }
        return results;
    }

    private int indexForValue(final int[] values, final int valueToFind) {
        for (int index = 0; index < values.length; index++) {
            if (values[index] == valueToFind) {
                return index;
            }
        }
        return -1;
    }

    @Override
    public void close() throws IOException {
        csFastaReader.close();
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
     *
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
}
