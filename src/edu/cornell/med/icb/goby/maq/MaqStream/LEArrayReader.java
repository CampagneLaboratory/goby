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

package edu.cornell.med.icb.goby.maq.MaqStream;

import it.unimi.dsi.lang.MutableString;

/**
 * A class to read little endian (LE) data from a byte array. This is useful when you have
 * a record of data where every record is the same size and can be read in one go into
 * a byte[]. This class will then read the LE data from that byte[].  This is designed so
 * the byte array can be filled multiple times and read over and over (just call reset()
 * after you have read a new record / re-filled the array).
 *
 * @author Kevin Dorff
 */
public class LEArrayReader {
    /**
     * The position we are currently reading from.
     */
    public int position;

    /**
     * The data we are reading.
     */
    public final byte[] data;

    /**
     * The data we are reading.
     */
    public char[] readCharBuffer;

    /**
     * Constructor.
     *
     * @param data the byte[] that will be repeatedly used to read LE data
     */
    public LEArrayReader(final byte[] data) {
        this.data = data;
        this.position = 0;
    }

    /**
     * Reset the read position back to the start of the array.
     */
    public void reset() {
        position = 0;
    }

    /**
     * Skip a specific number of bytes forward.
     * @param numBytes number of bytes to skip
     */
    public void skipBytes(final int numBytes) {
        position += numBytes;
    }

    /**
     * Read a little endian unsigned 8 bit number. As the number is unsigned 8 bits it
     * cannot be returned in a byte, but is instead returned in a short (16 bits, signed).
     *
     * @return the 8 bit unsigned number
     */
    public short readUInt8() {
        return data[position++];
    }


    /**
     * Read a little endian unsigned 32 bit number. As the number is unsigned 32 bits it
     * cannot be returned in a int, but is instead returned in a long (64 bits, signed).
     *
     * @return the 32 bit unsigned number
     */
    public long readUInt32() {
        return ((data[position++] & 0xff)
                | (data[position++] & 0xff) << 8
                | (data[position++] & 0xff) << 16
                | (data[position++] & 0xff) << 24);
    }

    /**
     * Read an int, 32-bits. Like DataInputStream.readInt except little endian.
     *
     * @return little-endian binary int from the datastream
     */
    public int readInt() {
        return ((data[position++] & 0xff)
                | ((data[position++] & 0xff) << 8)
                | ((data[position++] & 0xff) << 16)
                | ((data[position++]) << 24));
    }

    /**
     * Read a varaible length string. This is a single (binary, little endian) integer
     * followed by that many bytes of ASCII characters. This will return a java String.
     * Any trailing NULL characters ('\0') will be stripped.
     *
     * @return the String
     */
    public MutableString readVariableLengthString() {
        return readFixedLengthString(readInt());
    }

    /**
     * Read a fixed length length string. This will read the specified number of ASCII
     * characters and return as a MutableString. Any trailing NULL characters ('\0')
     * will be stripped.
     *
     * @param length the number of 8 bit ASCII chars to read from the data stream.
     * @return the String
     */
    public MutableString readFixedLengthString(final int length) {
        final MutableString temp = decodeToMutableStringAscii(data, position, length);
        int omit = 0;
        for (int i = length - 1; i >= 0; i--) {
            // Removed padding '\0's at the end of the string
            if (temp.charAt(i) == 0) {
                omit++;
            } else {
                break;
            }
        }
        temp.length(length - omit);
        position += length;
        return temp;
    }

    /**
     * Decode a string from offset off with length len.
     *
     * This version ONLY works if the string being read is normal ascii but in our
     * case this should work fine.
     *
     * @param ba the byte array (source)
     * @param off the offset where the string starts
     * @param len the length of the data to be read
     * @return the string value
     */
    private synchronized MutableString decodeToMutableStringAscii(
            final byte[] ba, final int off, final int len) {
        final MutableString dest;
        if (len == 0) {
            return new MutableString(0);
        } else {
            dest = new MutableString(len);
        }
        final int end = off + len;
        for (int pos = off; pos < end; pos++) {
            dest.append((char)(ba[pos] & 0xFF));
        }
        return dest;
    }
}
