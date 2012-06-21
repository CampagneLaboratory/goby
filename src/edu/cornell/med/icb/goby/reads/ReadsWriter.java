/*
 * Copyright (C) 2009-2012 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.cornell.med.icb.goby.reads;

import java.io.Closeable;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Properties;

/**
 * @author Fabien Campagne
 *         Date: 6/21/12
 *         Time: 1:50 PM
 */
public interface ReadsWriter extends Closeable {
    void setQualityScores(byte[] qualityScores);

    void setDescription(CharSequence description);

    void setSequence(CharSequence sequence);

    void setPairSequence(CharSequence sequence);

    void appendEntry(CharSequence description,
                     CharSequence sequence,
                     byte[] qualityScores) throws IOException;

    void appendEntry(CharSequence description,
                     CharSequence sequence) throws IOException;

    void appendEntry(CharSequence sequence) throws IOException;

    void appendEntry(Reads.ReadEntry.Builder entryBuilder) throws IOException;

    void close() throws IOException;

    void appendEntry() throws IOException;

    void appendEntry(int readIndex) throws IOException;

    void setNumEntriesPerChunk(int numEntriesPerChunk);

    void setIdentifier(CharSequence identifier);

    long getSequenceBasesWritten();

    void printStats(PrintStream out);

    void setBarcodeIndex(int barcodeIndex);

    void setQualityScoresPair(byte[] qualityScores);

    void appendMetaData(String key, String value);

    void setMetaData(Properties keyValuePairs);

    void setCodec(ReadCodec codec);
}
