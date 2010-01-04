package edu.cornell.med.icb.goby.reads;

/**
 * @author Fabien Campagne
 *         Date: Jun 20, 2009
 *         Time: 10:28:37 AM
 */
public final class SequenceEncoder {
   public void convertToByteBuffer(final byte[] sequence, final int referencePosition, final int readLength, final byte[] byteBuffer) {
        for (int i = 0; i < readLength; ++i) {
            byteBuffer[i] = (byte) sequence[i + referencePosition];
        }
    }
}
