package edu.cornell.med.icb.goby.reads;

import it.unimi.dsi.fastutil.booleans.BooleanArrayList;
import it.unimi.dsi.fastutil.booleans.BooleanList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.logging.ProgressLogger;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: Jun 20, 2009
 *         Time: 10:26:58 AM
 */
public class ReadsLoader {
    private ReadSet readIndexFilter;
    private MutableString sequence;

    SequenceDigests[] digests = null;
    private ObjectArrayList<byte[]> compressedReads;
    private int numReads;
    private File readsFile;
    private static final int MAX_PROCESS_READS = Integer.MAX_VALUE;

    int readLength;
    private byte[] byteBuffer = null;
    private int numberOfMismaches = 0;


    public int getReadLength() {
        return readLength;
    }

    public ReadsLoader(ReadSet readIndexFilter, File readsFile) {
        this.readIndexFilter = readIndexFilter;

        compressedReads = new ObjectArrayList<byte[]>(1000000);
        this.readsFile = readsFile;
        this.sequence = new MutableString();
        this.encoder = new SequenceEncoder();

    }

    boolean colorSpace;
    SequenceEncoder encoder;


    public void setColorSpace(boolean colorSpace) {
        this.colorSpace = colorSpace;
    }

    public int read() throws IOException {

        final ProgressLogger progress = new ProgressLogger();
        progress.displayFreeMemory = true;
        progress.start("parsing reads");
        final ReadsReader readsReader = new ReadsReader(new FileInputStream(readsFile));

        int maxReadLength = 0;
        int minReadLength = Integer.MAX_VALUE;
        for (final Reads.ReadEntry readEntry : readsReader) {

            if (readIndexFilter == null || readIndexFilter.contains(readEntry.getReadIndex())) {

                maxReadLength = Math.max(maxReadLength, readEntry.getReadLength());
                minReadLength = Math.min(minReadLength, readEntry.getReadLength());
                if (minReadLength != maxReadLength) {
                    System.err.println("The read length must be fixed.");
                    System.exit(1);
                }
                if (digests == null) {
                    if (numberOfMismaches == 0) {
                        digests = new SequenceDigests[2];
                        digests[0] = new SequenceDigests(minReadLength, true);     // forward strand
                        digests[1] = new SequenceDigests(minReadLength, false);    // reverse strand

                    } else if (numberOfMismaches == 1) {
                        digests = new SequenceDigests[4];
                        BooleanList mask1 = new BooleanArrayList();
                        BooleanList mask2 = new BooleanArrayList();
                        for (int i = 0; i < readLength; i++) {

                            mask1.set(i, i % 2 == 1);
                            mask2.set(i, !mask1.get(i));
                        }
                        digests[0] = new SequenceDigests(mask1, minReadLength, true);
                        digests[1] = new SequenceDigests(mask1, minReadLength, false);
                        digests[2] = new SequenceDigests(mask2, minReadLength, true);
                        digests[3] = new SequenceDigests(mask2, minReadLength, false);
                    }

                }


                byteBuffer = new byte[maxReadLength]; //2 bits per base require four times less space.

                final int readIndex = readEntry.getReadIndex();
                final byte[] bytes = readEntry.getSequence().toByteArray();
                for (SequenceDigests digest : digests) {

                    digest.digestAndStore(bytes, 0, readIndex);
                }
                System.arraycopy(bytes, 0, byteBuffer, 0, maxReadLength);

                if (readIndex > compressedReads.size()-1) {

                    compressedReads.size((readIndex + 1)*3/2);
                }
                compressedReads.set(readIndex, byteBuffer);

                progress.lightUpdate();

                if (numReads > MAX_PROCESS_READS) break;
            }
            numReads++;
        }

        readLength = minReadLength;
        compressedReads.size(numReads+1);
        readsReader.close();
        progress.stop("Finished parsing reads.");
        return numReads;
    }

    public ObjectList<byte[]> getCompressedReads() {
        return compressedReads;
    }


    public byte[] getByteBuffer() {
        return byteBuffer;
    }

    public SequenceDigests[] getDigests() {
        return digests;
    }
}
