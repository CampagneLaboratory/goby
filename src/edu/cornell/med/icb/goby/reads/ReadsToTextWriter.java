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

import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import it.unimi.dsi.lang.MutableString;

import java.io.IOException;
import java.io.PrintStream;
import java.util.Properties;

/**
 * Writes reads to text. Useful for debugging and for writing more intelligible Junit tests.
 *
 * @author Fabien Campagne
 *         Date: 6/21/12
 *         Time: 1:51 PM
 */
public class ReadsToTextWriter implements ReadsWriter {
    private MutableString textOutput = new MutableString();

    public String getTextOutput() {
        return textOutput.toString();
    }

    @Override
    public void setQualityScores(final byte[] qualityScores) {
        final String fieldName = "setQualityScores";
        final String value = ByteArrayList.wrap(qualityScores).toString();
        printSetField(fieldName, value);
    }

    private void printSetField(final String fieldName, final CharSequence value) {
        textOutput.append(fieldName);
        textOutput.append(": ");
        textOutput.append(value.toString());
        textOutput.append("\n");
    }

    @Override

    public void setDescription(final CharSequence description) {
        printSetField("description", description);
    }

    @Override
    public void setSequence(final CharSequence sequence) {
        printSetField("sequence", sequence);
    }

    @Override
    public void setPairSequence(final CharSequence sequence) {
        printSetField("sequence", sequence);
    }

    @Override
    public void appendEntry(final CharSequence description, final CharSequence sequence, final byte[] qualityScores) throws IOException {
        printAppend(description, sequence, qualityScores);
    }

    private void printAppend(Object... o) {
        textOutput.append("append: ");
        for (int i = 0; i < o.length; i += 2) {
            String fieldName = (String) o[i];
            Object argument = o[i + 1];
            printSetField(fieldName, argument.toString());
        }
        textOutput.append("\n");
    }

    @Override
    public void appendEntry(final CharSequence description, final CharSequence sequence) throws IOException {
        printAppend("description", description, "sequence", sequence);
    }

    @Override
    public void appendEntry(final CharSequence sequence) throws IOException {
        printAppend("sequence", sequence);
    }

    @Override
    public void appendEntry(final Reads.ReadEntry.Builder entryBuilder) throws IOException {
        printAppend("entryBuilder", entryBuilder);
    }

    @Override
    public void close() throws IOException {
       textOutput.append("close()\n");
    }

    @Override
    public void appendEntry() throws IOException {
        printAppend();
    }

    @Override
    public void appendEntry(final int readIndex) throws IOException {
        printAppend("readIndex", readIndex);
    }

    @Override
    public void setNumEntriesPerChunk(final int numEntriesPerChunk) {
        printSetField("numEntriesPerChunk", new Integer(numEntriesPerChunk).toString());
    }

    @Override
    public void setIdentifier(final CharSequence identifier) {
        printSetField("identifier", identifier);
    }

    @Override
    public long getSequenceBasesWritten() {
        printGetField("SequenceBasesWritten");
        return -1;
    }

    private void printGetField(final String getterName) {
        this.textOutput.append("get");
        this.textOutput.append(getterName);
        this.textOutput.append("\n");
    }

    @Override
    public void printStats(final PrintStream out) {
        textOutput.append("printStats\n");
    }

    @Override
    public void setBarcodeIndex(final int barcodeIndex) {
        printSetField("barcodeIndex", new Integer(barcodeIndex).toString());
    }

    @Override
    public void setQualityScoresPair(final byte[] qualityScores) {
        printSetField("barcodeIndex", ByteArrayList.wrap(qualityScores).toString());
    }

    @Override
    public void appendMetaData(final String key, final String value) {
        textOutput.append("appendMetaData: ");
        textOutput.append(key);
        textOutput.append("=");
        textOutput.append(value);
        textOutput.append("\n");
    }

    @Override
    public void setMetaData(final Properties keyValuePairs) {
        printSetField("MetaData", keyValuePairs.toString());
    }

    @Override
    public void setCodec(final ReadCodec codec) {
        printSetField("codec", codec.name());

    }
}
