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

import com.mindprod.ledatastream.LEDataOutputStream;

import java.io.IOException;
import java.io.OutputStream;

/**
 * Describe class here.
 *
 * @author Kevin Dorff
 */
public class MaqLEDataOutputStream extends LEDataOutputStream {
    /**
     * constructor.
     *
     * @param out the outputstream we write little endian binary data onto.
     */
    public MaqLEDataOutputStream(final OutputStream out) {
        super(out);
    }

    /**
     * Write an unsigned 32-bit integer.
     *
     * @param v the long to write
     * @throws IOException if write fails.
     */
    public void writeUnsignedInteger(final long v) throws IOException {
        work[0] = (byte) v;
        work[1] = (byte) (v >> 8);
        work[2] = (byte) (v >> 16);
        work[3] = (byte) (v >> 24);
        write(work, 0, 4);
    }

    /**
     * Write a little endian unsigned 8 bit number. As the number is unsigned 8 bits it
     * cannot be returned in a byte, but is instead returned in a short (16 bits, signed).
     * @param value the 8 bit unsigned number
     * @throws java.io.IOException error writing
     */
    public void writeUInt8(final int value) throws IOException {
        write(value);
    }

    /**
     * Write a little endian unsigned 32 bit number. As the number is unsigned 32 bits it
     * cannot be returned in a int, but is instead returned in a long (64 bits, signed).
     * @param value the 32 bit unsigned number
     * @throws java.io.IOException error writing
     */
    public void writeUInt32(final long value) throws IOException {
        writeUnsignedInteger(value);
    }

    /**
     * Write a varaible length string. This will write a (binary, little endian, signed) integer
     * which is the size of the string followed by that many bytes of ASCII characters. This WILL
     * append a '\0' to the end of the data that is written and the length that is written will be
     * increased by one to account for the '\0' that is written.
     * followed by that many bytes of ASCII characters.
     * @param value the string to write
     * @throws java.io.IOException error writing
     */
    public void writeVariableLengthString(final String value) throws IOException {
        final int length = value.length() + 1;
        writeInt(length);
        writeFixedLengthString(value, length);
    }

    /**
     * Write a fixed length string. Will be tail padded with '\0'. If valueStr cannot
     * fit within length and still include at least one trailing '\0' the string will
     * be shorted.
     * @param valueStr the number of 8 bit ASCII chars to write from the data stream.
     * @param length the number of 8 bit ASCII chars to write from the data stream.
     * @throws java.io.IOException error writing
     */
    public void writeFixedLengthString(final String valueStr, final int length) throws IOException {
        final byte[] out = new byte[length];
        final byte[] outBytes;
        if (valueStr.length() > (length - 1)) {
            // Make sure the string fits within length and can still have a trailing '\0'
            outBytes = valueStr.substring(0, length - 1).getBytes();
        } else {
            // It fits
            outBytes = valueStr.getBytes();
        }
        System.arraycopy(outBytes, 0, out, 0, outBytes.length);
        write(out);
    }
}
