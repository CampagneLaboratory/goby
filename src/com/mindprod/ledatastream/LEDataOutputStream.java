package com.mindprod.ledatastream;

import java.io.DataOutput;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.Closeable;

/**
 * This class is the same as the one distributed by Roedy Green, Canadian Mind* Products
 * EXCEPT I have removed the "final" declaration on the class and methods to make it
 * extendable and I have reformatted the code. I also added Closable. Kevin Dorff, ICB-WMC.
 *
 * <pre>
 *       LEDataOutputStream.java
 * <p/>
 *        copyright (c) 1998-2009 Roedy Green,
 *        Canadian Mind Products
 *       #101 - 2536 Wark Street
 *       Victoria,  BC Canada V8T 4G8
 *       hel: (250) 361-9093
 *       roedy g at mindprod dotcom
 *       http://mindprod.com
 * <p/>
 *        Version 1.0 1998-01-6
 * <p/>
 *       1.1 1998-01-07-officially implements DataInput
 * <p/>
 *       1.2 1998-01-09- add LERandomAccessFile
 * <p/>
 *       1.3 1998-08-28 1.4 1998-11-10 - add new address and phone.
 * <p/>
 *       1.5 1999-10-08- use com.mindprod.ledatastream
 *       package name. Very similar to DataOutputStream except it writes
 * little-endian
 *       instead of big-endian binary data. We can't extend DataOutputStream
 * directly
 *       since it has only final methods. This forces us implement
 * LEDataOutputStream
 *       with a DataOutputStream object, and use wrapper methods.
 * </pre>
 *
 * @noinspection WeakerAccess
 */
public class LEDataOutputStream implements DataOutput, Closeable {
    // ------------------------------ FIELDS ------------------------------

    /**
     * undisplayed copyright notice.
     *
     * @noinspection UnusedDeclaration
     */
    private static final String EMBEDDED_COPYRIGHT =
            "copyright (c) 1999-2009 Roedy Green, Canadian Mind Products, http://mindprod.com";

    /**
     * to get at big-Endian write methods of DataOutPutStream.
     *
     * @noinspection WeakerAccess
     */
    protected final DataOutputStream dos;

    /**
     * work array for composing output.
     *
     * @noinspection WeakerAccess
     */
    protected final byte[] work;

    // -------------------------- PUBLIC INSTANCE  METHODS --------------------------
    /**
     * constructor.
     *
     * @param out the outputstream we write little endian binary data onto.
     */
    public LEDataOutputStream(final OutputStream out) {
        this.dos = new DataOutputStream(out);
        work = new byte[8]; // work array for composing output
    }

    /**
     * Close stream.
     *
     * @throws IOException if close fails.
     */
    public void close() throws IOException {
        dos.close();
    }

    /**
     * Flush stream without closing.
     *
     * @throws IOException if flush fails.
     */
    public void flush() throws IOException {
        dos.flush();
    }

    /**
     * Get size of stream.
     *
     * @return bytes written so far in the stream. Note this is a int, not a long as you
     * would exect. This because the underlying DataInputStream has a design flaw.
     */
    public int size() {
        return dos.size();
    }

    /**
     * This method writes only one byte, even though it says int (non-Javadoc).
     *
     * @param ib the byte to write.
     * @throws IOException if write fails.
     * @see java.io.DataOutput#write(int)
     */
    public synchronized void write(final int ib) throws IOException {
        dos.write(ib);
    }

    /**
     * Write out an array of bytes.
     *
     * @param ba the data
     * @throws IOException if write fails.
     * @see java.io.DataOutput#write(byte[])
     */
    public void write(final byte[] ba) throws IOException {
        dos.write(ba, 0, ba.length);
    }

    /**
     * Writes out part of an array of bytes.
     *
     * @param ba the data.
     * @param off the start offset in the data.
     * @param len the number of bytes to write.
     * @throws IOException if write fails.
     * @see java.io.DataOutput#write(byte[],int,int)
     */
    public synchronized void write(
            final byte[] ba, final int off, final int len) throws IOException {
        dos.write(ba, off, len);
    }

    /**
     * Write a booleans as one byte. Only writes one byte.
     *
     * @param v boolean to write.
     * @throws IOException if write fails.
     * @see java.io.DataOutput#writeBoolean(boolean)
     */
    public void writeBoolean(final boolean v) throws IOException {
        dos.writeBoolean(v);
    }

    /**
     * write a byte.
     *
     * @param v the byte to write.
     * @throws IOException if write fails.
     * @see java.io.DataOutput#writeByte(int)
     */
    public void writeByte(final int v) throws IOException {
        dos.writeByte(v);
    }

    /**
     * Write a string.
     *
     * @param s the string to write.
     * @throws IOException if write fails.
     * @see java.io.DataOutput#writeBytes(java.lang.String)
     */
    public void writeBytes(final String s) throws IOException {
        dos.writeBytes(s);
    }

    /**
     * Write a char. Like DataOutputStream.writeChar. Note the parm is an int even though
     * this as a writeChar
     *
     * @param v the char to write
     * @throws IOException if write fails.
     */
    public void writeChar(final int v) throws IOException {
        // same code as writeShort
        work[0] = (byte) v;
        work[1] = (byte) (v >> 8);
        dos.write(work, 0, 2);
    }

    /**
     * Write a string, not a char[]. Like DataOutputStream.writeChars, flip endianness of each char.
     *
     * @param s the String containing the characters to write
     * @throws IOException if write fails.
     */
    public void writeChars(final String s) throws IOException {
        int len = s.length();
        for (int i = 0; i < len; i++) {
            writeChar(s.charAt(i));
        }
    } // end writeChars

    /**
     * Write a double.
     *
     * @param v the double to write. Like DataOutputStream.writeDouble.
     * @throws IOException if write fails.
     */
    public void writeDouble(final double v) throws IOException {
        writeLong(Double.doubleToLongBits(v));
    }

    /**
     * Write a float. Like DataOutputStream.writeFloat.
     *
     * @param v the float to write.
     * @throws IOException if write fails.
     */
    public void writeFloat(final float v) throws IOException {
        writeInt(Float.floatToIntBits(v));
    }

    /**
     * Write an int, 32-bits.  Like DataOutputStream.writeInt.
     *
     * @param v the int to write
     * @throws IOException if write fails.
     */
    public void writeInt(final int v) throws IOException {
        work[0] = (byte) v;
        work[1] = (byte) (v >> 8);
        work[2] = (byte) (v >> 16);
        work[3] = (byte) (v >> 24);
        dos.write(work, 0, 4);
    }

    /**
     * Write a long, 64-bits. like DataOutputStream.writeLong.
     *
     * @param v the long to write
     * @throws IOException if write fails.
     */
    public void writeLong(final long v) throws IOException {
        work[0] = (byte) v;
        work[1] = (byte) (v >> 8);
        work[2] = (byte) (v >> 16);
        work[3] = (byte) (v >> 24);
        work[4] = (byte) (v >> 32);
        work[5] = (byte) (v >> 40);
        work[6] = (byte) (v >> 48);
        work[7] = (byte) (v >> 56);
        dos.write(work, 0, 8);
    }

    /**
     * Write short, 16-bits. Like DataOutputStream.writeShort. also acts as a writeUnsignedShort
     *
     * @param v the short you want written in little endian binary format
     * @throws IOException if write fails.
     */
    public void writeShort(final int v) throws IOException {
        work[0] = (byte) v;
        work[1] = (byte) (v >> 8);
        dos.write(work, 0, 2);
    }

    /**
     * Write a string as a UTF counted string.
     *
     * @param s the string to write.
     * @throws IOException if write fails.
     * @see java.io.DataOutput#writeUTF(java.lang.String)
     */
    public void writeUTF(final String s) throws IOException {
        dos.writeUTF(s);
    }

} // end LEDataOutputStream
