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

package edu.cornell.med.icb.goby.util;

import net.sf.samtools.util.BlockCompressedInputStream;

import java.io.*;
import java.nio.CharBuffer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

/**
 * This class takes a reader and a pattern and removes lines that don't match the pattern.
 * Line terminators are converted to a \n.
 * Date: 1/2/12
 * Time: 11:26 AM
 */

public class GrepReader extends FilterReader {
    // This variable holds the current line.
    // If null and emitNewline is false, a newline must be fetched.
    String curLine;

    // This is the index of the first unread character in curLine.
    // If at any time curLineIx == curLine.length, curLine is set to null.
    int curLineIx;

    // If true, the newline at the end of curLine has not been returned.
    // It would have been more convenient to append the newline
    // onto freshly fetched lines. However, that would incur another
    // allocation and copy.
    boolean emitNewline;

    // Matcher used to test every line
    Matcher matcher;

    public GrepReader(String filename, String patternStr) throws IOException {
        super(filename.endsWith(".gz") ?
                new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream(filename)))) :
               new BufferedReader( new FileReader(filename)));


        Pattern pattern = Pattern.compile(patternStr);
        matcher = pattern.matcher("");
    }

    // This overridden method fills cbuf with characters read from in.
    public int read(char cbuf[], int off, int len) throws IOException {
        // Fetch new line if necessary
        if (curLine == null && !emitNewline) {
            getNextLine();
        }

        // Return characters from current line
        if (curLine != null) {
            int num = Math.min(len, Math.min(cbuf.length - off,
                    curLine.length() - curLineIx));
            // Copy characters from curLine to cbuf
            for (int i = 0; i < num; i++) {
                cbuf[off++] = curLine.charAt(curLineIx++);
            }

            // No more characters in curLine
            if (curLineIx == curLine.length()) {
                curLine = null;

                // Is there room for the newline?
                if (num < len && off < cbuf.length) {
                    cbuf[off++] = '\n';
                    emitNewline = false;
                    num++;
                }
            }

            // Return number of character read
            return num;
        } else if (emitNewline && len > 0) {
            // Emit just the newline
            cbuf[off] = '\n';
            emitNewline = false;
            return 1;
        } else if (len > 0) {
            // No more characters left in input reader
            return -1;
        } else {
            // Client did not ask for any characters
            return 0;
        }
    }

    // Get next matching line
    private void getNextLine() throws IOException {
        curLine = ((BufferedReader) in).readLine();
        while (curLine != null) {
            matcher.reset(curLine);
            if (!matcher.find()) {
                emitNewline = true;
                curLineIx = 0;
                return;
            }
            curLine = ((BufferedReader) in).readLine();
        }
        return;
    }

    public boolean ready() throws IOException {
        return curLine != null || emitNewline || in.ready();
    }

    public boolean markSupported() {
        return false;

    }
}
