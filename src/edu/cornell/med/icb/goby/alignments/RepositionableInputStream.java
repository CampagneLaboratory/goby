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

package edu.cornell.med.icb.goby.alignments;

import it.unimi.dsi.fastutil.io.RepositionableStream;
import org.apache.log4j.Logger;

import java.io.*;
import java.net.*;
import java.nio.channels.FileChannel;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Fabien Campagne
 *         Date: 2/22/12
 *         Time: 11:48 AM
 */
public class RepositionableInputStream extends InputStream implements RepositionableStream {

    private InputStream delegate;
    private String resource;
    private long currentPosition = 0;
    private long markedPosition;
    private boolean isLocalFile;
    // Only populated if the resource is a local file.
    private FileChannel channel;

    public RepositionableInputStream(final String resource) throws IOException {
        this.resource = resource;

        delegate = getStream(resource);
        if (RepositionableInputStream.isLocalFile(resource)) {
            isLocalFile = true;
            channel = ((FileInputStream) delegate).getChannel();
        }
    }

    @Override
    public int read() throws IOException {
        final int read = delegate.read();
        currentPosition += 1;
        return read;
    }

    @Override
    public int read(byte[] bytes) throws IOException {
        int read = delegate.read(bytes);
        currentPosition += read;
        return read;
    }

    @Override
    public int read(byte[] bytes, int i, int i1) throws IOException {

        int read = delegate.read(bytes, i, i1);
        currentPosition += read;
        return read;
    }

    @Override
    public long skip(long l) throws IOException {
        long skip = delegate.skip(l);
        currentPosition += l;

        return skip;
    }

    @Override
    public int available() throws IOException {
        return delegate.available();
    }

    @Override
    public void close() throws IOException {
        if (delegate != null) {
            delegate.close();
            delegate = null;
        }
        if (isLocalFile && channel != null) {
            channel.close();
            channel = null;

        }

    }

    @Override
    public synchronized void mark(int i) {

        markedPosition = currentPosition;
    }

    @Override
    public synchronized void reset() throws IOException {
        delegate.close();
        delegate = null;
        delegate = getStream(resource, markedPosition);
        currentPosition = markedPosition;
    }

    @Override
    public boolean markSupported() {
        return true;
    }

    @Override
    public void position(long p) throws IOException {
        if (isLocalFile) {
            channel.position(p);
        } else {
            if (p != currentPosition) {

                if (LOG.isTraceEnabled()) {

                    LOG.trace(String.format("repositioning URL stream to byte position= %d%n", p));
                }
                currentPosition = p;
                delegate.close();
                delegate = null;
                delegate = getStream(resource, p);
                markedPosition = p;

            }
        }
    }

    @Override
    public long position() throws IOException {
        return currentPosition;

    }

    public static InputStream getStream(String resource) throws IOException {
        return getStream(resource, 0);
    }

    private static Logger LOG = Logger.getLogger(RepositionableInputStream.class);

    public static InputStream getStream(String resource, long startOffset) throws IOException {

        try {
            final URL url = new URL(resource);
            if ("file".equals(url.getProtocol())) {
                resource = resource.replaceFirst("file://", "");
            } else {

                HttpURLConnection urlConnection = (HttpURLConnection) url.openConnection();
                urlConnection.setRequestMethod("GET");
                urlConnection.setRequestProperty("Range", "bytes="+Long.toString(startOffset)+"-");
                LOG.debug(String.format("Opening URL=%s at startOffset=%d", resource, startOffset));
                return new BufferedInputStream(urlConnection.getInputStream());
            }

        } catch (MalformedURLException e) {

        }
        // not a URL, try opening the file directly:
        return new FileInputStream(resource);

    }

    public static boolean isLocalFile(final String resource) {
        try {
            final URL url = new URL(resource);
            return "file".equals(url.getProtocol());
        } catch (MalformedURLException e) {
            return true;
        }

    }

    /**
     * Determine if a File.URL exists.
     *
     * @param resource
     */
    public static boolean resourceExist(String resource) {
        try {
            URL url = new URL(resource);
            if ("file".equals(url.getProtocol())) {
                resource = resource.replaceFirst("file://", "");
            } else {
                HttpURLConnection huc = (HttpURLConnection) url.openConnection();
                huc.setRequestMethod("GET");  //OR  huc.setRequestMethod ("HEAD");
                //HttpURLConnection.setFollowRedirects(false);    // may crash applets, not sure about JNLP
                huc.connect();
                final int code = huc.getResponseCode();

                return code == 200;
            }
        } catch (MalformedURLException e) {
            // not a URL
        } catch (ProtocolException e) {
            throw new RuntimeException("Unable to determe if URL exists, likely a non HTTP or file url: " + resource);
        } catch (IOException e) {

        }
        return new File(resource).exists();
    }

    public static RepositionableStream get(String resource) throws IOException {
        return new RepositionableInputStream(resource);
    }
}