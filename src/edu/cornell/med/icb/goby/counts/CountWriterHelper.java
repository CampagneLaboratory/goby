package edu.cornell.med.icb.goby.counts;

import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: Jun 12, 2009
 *         Time: 4:44:06 PM
 */
public class CountWriterHelper {
    CountsWriter delegate;

    public CountWriterHelper(CountsWriter delegate) {
        this.delegate = delegate;
    }

    int previousPosition = -1;

    int previousCount = 0;
    int lengthConstant = 1;
    int previousPositionNotWritten;

    public void appendCountAtPosition(final int count, final int position) throws IOException {

        System.out.printf("// count=%d position=%d previousCount=%d %n", count, position, previousCount);
        lengthConstant++;
        if (count == previousCount) {
            lengthConstant+=position-previousPositionNotWritten;
            previousPositionNotWritten=position;
        } else {

            delegate.appendCount(previousCount, lengthConstant);
            previousCount = count;
            lengthConstant =0 ;
            previousPosition = position;
            previousPositionNotWritten=position;
        }

    }

    public void close() throws IOException {
        //  if (previousCount != 0) delegate.appendCount(0, 1);
        delegate.close();
    }
}
