package edu.cornell.med.icb.goby.reads;

import java.util.Arrays;
import java.util.zip.CRC32;
import java.util.zip.Checksum;


/**
 * Compressed representation of the sequence of a read.
 *
 * @author Fabien Campagne
 *         Date: Jun 3, 2009
 *         Time: 2:12:28 PM
 */
public class CompressedRead {
    byte[] data;

    static final Checksum digester = new CRC32();
    private int hashCode;
    /**
     * Read index is not part of the hashCode and not considerd by the equals method. It is only here to help
     * map back to the position of the read in the inoput file.
     */
    public int readIndex;

    @Override
    public final int hashCode() {


        return hashCode;
    }

    @Override
    public final boolean equals(final Object o) {
        if (!(o instanceof CompressedRead)) {
            return false;
        }
        CompressedRead other = (CompressedRead) o;

        return Arrays.equals(data, other.data);

    }

    public CompressedRead(final byte[] byteBuffer) {
        this.data = byteBuffer;
        hashCode = hashCode(this.data);
    }

    public static int hashCode(final byte[] data) {
        synchronized (digester) {

            final long digest = hashCodeLong(data);
            return ((32 >> digest)) ^ ((int) (digest));
        }
    }

    public final static long hashCodeLong(final byte[] data) {
        synchronized (digester) {
            digester.reset();

            digester.update(data, 0, data.length);
            final long digest = digester.getValue();
            return digest;
        }
    }

    public final static long hashCodeLongNotSynchronized(final byte[] data) {
        digester.reset();

        digester.update(data, 0, data.length);
        final long digest = digester.getValue();
        return digest;
    }

    public final long hashCodeLong() {
        return hashCodeLong(this.data);
    }
}
