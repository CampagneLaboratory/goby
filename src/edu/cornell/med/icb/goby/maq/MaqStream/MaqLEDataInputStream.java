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

import com.mindprod.ledatastream.LEDataInputStream;

import java.io.IOException;
import java.io.InputStream;

/**
 * Slight augmentation / extention of LEDataInputStream. Aliases
 * for reading and writing UInt8 and UInt32 plus reading and writing
 * variable and fixed length strings.
 * @author Kevin Dorff
 */
public class MaqLEDataInputStream extends LEDataInputStream {

    /**
     * constructor.
     *
     * @param in binary inputstream of little-endian data.
     */
    public MaqLEDataInputStream(final InputStream in) {
        super(in);
    }

    /**
     * Read an unsigned integer, 32 bits. Like DataInputStream.readUnsignedShort except 32bits and
     * little endian. Note, returns long even though it reads a int.
     *
     * @return little-endian int from the stream.
     * @throws IOException if read fails.
     */
    public long readUnsignedInteger() throws IOException {
        readFully(work, 0, 4);
        return ((work[3] & 0xff) << 24
                | (work[2] & 0xff) << 16
                | (work[1] & 0xff) << 8
                | (work[0] & 0xff));
    }

    /**
     * Read a little endian unsigned 8 bit number. As the number is unsigned 8 bits it
     * cannot be returned in a byte, but is instead returned in a short (16 bits, signed).
     * @return the 8 bit unsigned number
     * @throws java.io.IOException error reading
     */
    public short readUInt8() throws IOException {
        return (short) readUnsignedByte();
    }

    /**
     * Read a little endian unsigned 32 bit number. As the number is unsigned 32 bits it
     * cannot be returned in a int, but is instead returned in a long (64 bits, signed).
     * @return the 32 bit unsigned number
     * @throws IOException error reading
     */
    public long readUInt32() throws IOException {
        return readUnsignedInteger();
    }

    /**
     * Read a varaible length string. This is a single (binary, little endian) integer
     * followed by that many bytes of ASCII characters. This will return a java String.
     * Any trailing NULL characters ('\0') will be stripped.
     * @return the String
     * @throws IOException error reading
     */
    public String readVariableLengthString() throws IOException {
        final int size = readInt();
        return readFixedLengthString(size);
    }

    /**
     * Read a fixed length length string. This will read the specified number of ASCII
     * characters and return as a Java String. Any trailing NULL characters ('\0')
     * will be stripped.
     * @param length the number of 8 bit ASCII chars to read from the data stream.
     * @return the String
     * @throws IOException error reading
     */
    public String readFixedLengthString(final int length) throws IOException {
        final byte[] data = new byte[length];
        readFully(data);
        final String temp = new String(data);
        int omit = 0;
        for (int i = data.length - 1; i >= 0; i--) {
            // Removed padding '\0's at the end of the string
            if (data[i] == 0) {
                omit++;
            } else {
                break;
            }
        }
        return temp.substring(0, temp.length() - omit);
    }

}
