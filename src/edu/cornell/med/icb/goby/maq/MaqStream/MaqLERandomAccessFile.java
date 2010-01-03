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

import com.mindprod.ledatastream.LERandomAccessFile;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

/**
 * Describe class here.
 *
 * @author Kevin Dorff
 */
public class MaqLERandomAccessFile extends LERandomAccessFile {

    /**
     * constructor.
     *
     * @param file file to read/write.
     * @param rw   like {@link java.io.RandomAccessFile} where "r" for read "rw" for read and
     * write, "rws" for read-write sync, and "rwd" for read-write dsync. Sync ensures the
     * physical i/o has completed before the method returns.
     * @throws java.io.FileNotFoundException if open fails.
     */
    public MaqLERandomAccessFile(final File file, final String rw) throws FileNotFoundException {
        super(file, rw);
    }

    /**
     * constructors.
     *
     * @param file name of file.
     * @param rw   string "r" or "rw" depending on read or read/write.
     * @throws java.io.FileNotFoundException if open fails.
     * @noinspection SameParameterValue
     */
    public MaqLERandomAccessFile(final String file, final String rw) throws FileNotFoundException {
        super(file, rw);
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
