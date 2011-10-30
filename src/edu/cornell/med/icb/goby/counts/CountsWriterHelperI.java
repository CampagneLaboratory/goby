package edu.cornell.med.icb.goby.counts;

import java.io.Closeable;
import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: 10/30/11
 *         Time: 12:20 PM
 */
public interface CountsWriterHelperI extends Closeable {
    void appendCountAtPosition(int count, int position) throws IOException;

    void close() throws IOException;
}
