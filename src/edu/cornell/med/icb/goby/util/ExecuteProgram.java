/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                         Weill Medical College of Cornell University
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

import edu.cornell.med.icb.goby.exception.GobyRuntimeException;
import edu.cornell.med.icb.util.ICBStringUtils;
import org.apache.commons.lang.RandomStringUtils;
import org.apache.commons.lang.SystemUtils;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

/**
 * @author Fabien Campagne
 *         Date: Jul 9, 2009
 *         Time: 1:49:35 PM
 */
public class ExecuteProgram {
    /**
     * Used to log debug and informational messages.
     */
    private static final Logger log = Logger.getLogger(ExecuteProgram.class);

    public void execute(final String command) throws InterruptedException, IOException {
        execute(command, System.out);
    }

    public void executeToLog(
            final String command, final Class logAsClass, final Level logLevel, final String prefix)
            throws InterruptedException, IOException {
        final LoggingOutputStream loggingOutputStream =
                new LoggingOutputStream(logAsClass, logLevel, prefix);
        execute(command, loggingOutputStream, loggingOutputStream, SystemUtils.IS_OS_UNIX);
        loggingOutputStream.close();
    }

    public void executeStderrToLog(
            final String command, final OutputStream stdOutOutput, final Class logAsClass,
            final Level logLevel, final String prefix) throws InterruptedException, IOException {
        final LoggingOutputStream loggingOutputStream =
                new LoggingOutputStream(logAsClass, logLevel, prefix);
        execute(command, stdOutOutput, loggingOutputStream, SystemUtils.IS_OS_UNIX);
        loggingOutputStream.close();
    }

    boolean ignoreFailure;

    public void setIgnoreFailure(final boolean ignoreFailure) {
        this.ignoreFailure = ignoreFailure;
    }

    boolean silent;

    public void setSilent(final boolean silent) {
        this.silent = silent;
    }

    public void execute(final String command, final OutputStream stdOutOutput)
            throws InterruptedException, IOException {
        execute(command, stdOutOutput, System.err, SystemUtils.IS_OS_UNIX);
    }

    public void execute(
            String command, final OutputStream stdOutOutput, final OutputStream stdErrOutput,
            final boolean nice) throws InterruptedException, IOException {
        Process proc = null;
        try {
            // nice each job to preserve system responsiveness
            if (nice) {
                command = "nice " + command;
            }

            final String execTag = RandomStringUtils.randomAlphabetic(7);
            if (stdOutOutput instanceof LoggingOutputStream) {
                ((LoggingOutputStream) stdOutOutput).setOutputTag(execTag);
            }
            if (stdErrOutput instanceof LoggingOutputStream) {
                ((LoggingOutputStream) stdErrOutput).setOutputTag(execTag);
            }

            if (!silent) {
                log.info(String.format("[%s] About to execute command %s", execTag, command));
            }
            final long start = System.currentTimeMillis();
            proc = Runtime.getRuntime().exec(command);

            final Thread[] threads;
            threads = consumeProcessOutput(proc, stdOutOutput, stdErrOutput);
            proc.waitFor();
            // we wait for the output and error thread dumpers to finish:
            threads[1].join();
            threads[0].join();

            final String duration = ICBStringUtils.millis2hms(System.currentTimeMillis() - start);
            if (!silent) {
                log.info(String.format("[%s] Finshed executing command %s in %s",
                        execTag, command, duration));
            }

            final int error = proc.exitValue();
            if (error != 0) {
                if (ignoreFailure) {
                    log.warn("Received status $error executing command '" + command + "', ignoring");
                } else {
                    final String message = String.format("Error %d executing command '%s'", error, command);
                    log.error(message);
                    throw new GobyRuntimeException(message);
                }
            }
        } catch (InterruptedException e) {
            // We just got signaled this JVM is being killed. Kill whichever process we were executing..
            if (proc != null) {
                proc.destroy();
            }
            throw e;
        } catch (IOException e) {
            log.error(e);
            throw e;
        }
    }


    /**
     * Gets the output and error streams from a process and reads them
     * to keep the process from blocking due to a full output buffer.
     * The processed stream data is appended to the supplied OutputStream.
     * For this, two Threads are started, so this method will return immediately.
     *
     * @param self   a Process
     * @param output an OutputStream to capture the process stdout
     * @param error  an OutputStream to capture the process stderr
     * @since 1.5.2
     */
    public static Thread[] consumeProcessOutput(final Process self, final OutputStream output, final OutputStream error) {
        final Thread[] threads = new Thread[2];
        threads[0] = consumeProcessOutputStream(self, output);
        threads[1] = consumeProcessErrorStream(self, error);
        return threads;
    }

    /**
     * Gets the error stream from a process and reads it
     * to keep the process from blocking due to a full buffer.
     * The processed stream data is appended to the supplied OutputStream.
     * A new Thread is started, so this method will return immediately.
     *
     * @param self a Process
     * @param err  an OutputStream to capture the process stderr
     * @since 1.5.2
     */
    public static Thread consumeProcessErrorStream(final Process self, final OutputStream err) {
        final Thread thread = new Thread(new ByteDumper(self.getErrorStream(), err));
        thread.start();
        return thread;
    }

    /**
     * Gets the output stream from a process and reads it
     * to keep the process from blocking due to a full output buffer.
     * The processed stream data is appended to the supplied StringBuffer.
     * A new Thread is started, so this method will return immediately.
     *
     * @param self   a Process
     * @param output a StringBuffer to capture the process stdout
     * @since 1.5.2
     */
    public static Thread consumeProcessOutputStream(final Process self, final OutputStream output) {
        final Thread thread = new Thread(new ByteDumper(self.getInputStream(), output));
        thread.start();
        return thread;
    }

    /**
     * From DetaultGroovyMethods.
     */
    private static class ByteDumper implements Runnable {
        InputStream in;
        OutputStream out;

        public ByteDumper(final InputStream in, final OutputStream out) {
            this.in = new BufferedInputStream(in);
            this.out = out;
        }

        public ByteDumper(final InputStream in) {
            this.in = new BufferedInputStream(in);
        }

        public void run() {
            final byte[] buf = new byte[1024 * 1024];
            int next;
            try {
                while ((next = in.read(buf)) != -1) {
                    if (out != null) {
                        out.write(buf, 0, next);
                        out.flush();
                    }
                }
                if (out != null) {
                    out.flush();
                }
            } catch (IOException e) {
                throw new GobyRuntimeException("exception while dumping process stream", e);
            }
        }
    }
}
