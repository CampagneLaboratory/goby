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

package edu.cornell.med.icb.goby.reads;

import com.google.protobuf.ByteString;
import it.unimi.dsi.fastutil.booleans.BooleanList;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;

/**
 * Calculates and keeps track of sequence digests.
 *
 * @author Fabien Campagne
 *         Date: Jun 20, 2009
 *         Time: 2:29:43 PM
 */
public class SequenceDigests {
    private final Int2ObjectMap<int[]> readDigestsToSequenceIndex;
    private final boolean matchesPositiveStrand;
    private final boolean[] mask;
    private final int readLength;
    private final int step = 0;
    private static final int[] NIL = new int[0];


    /**
     * Construct digests over the bases active in mask.
     * Mask is a list of booleans, where each boolean indicates if a base of the sequence (at the corresponding position of the mask)
     * should be included in the digest.
     *
     * @param mask
     * @param readLength
     */
    public SequenceDigests(final BooleanList mask, final int readLength, final boolean matchesPositiveStrand) {
        this(mask.toBooleanArray(), readLength, matchesPositiveStrand, new Int2ObjectOpenHashMap<int[]>());


    }

    private SequenceDigests(final boolean[] mask, final int readLength, final boolean matchesPositiveStrand, final Int2ObjectMap<int[]> readDigests) {
        readDigestsToSequenceIndex = readDigests;
        readDigestsToSequenceIndex.defaultReturnValue(NIL);
        this.mask = mask;
        this.readLength = readLength;
        this.matchesPositiveStrand = matchesPositiveStrand;
    }


    /**
     * Construct digests over the entire length of the sequence.
     *
     * @param readLength
     */
    public SequenceDigests(final int readLength, final boolean matchesPositiveStrand) {
        this(null, readLength, matchesPositiveStrand, new Int2ObjectOpenHashMap<int[]>());

    }


    public final int digest(final byte[] sequence, final int offset, final int length) {

        return mask == null ?
                (matchesPositiveStrand ?
                        crc32(sequence, offset, length) :
                        crc32ReverseComplement(sequence, offset, length)) :

                crc32WithMask(sequence, offset, length);


    }

    public final int digestDirectStrandOnly(final byte[] sequence, final int offset, final int length) {

        return mask == null ?
                crc32(sequence, offset, length) :
                crc32WithMask(sequence, offset, length);


    }

    public final int digest(final byte[] sequence, final int offset) {

        return mask == null ?
                (matchesPositiveStrand ?
                        crc32(sequence, offset, readLength) :
                        crc32ReverseComplement(sequence, offset, readLength)) :

                crc32WithMask(sequence, offset, readLength);

    }

    public final int digest(final ByteString sequence, final int offset, final int length) {
        return mask == null ?
                (matchesPositiveStrand ?
                        crc32(sequence, offset, length) :
                        crc32ReverseComplement(sequence, offset, length)) :

                crc32WithMask(sequence, offset, length);
    }


    public final void digestAndStore(final byte[] sequence, final int offset, final int sequenceIndex) {
        final int digest = digest(sequence, offset);
        final int[] previous = readDigestsToSequenceIndex.get(digest);

        final int[] newArray = new int[previous.length + 1];
        int i = 0;
        for (; i < previous.length; ++i) {
            if (newArray[i] == sequenceIndex) {
                // sequence index is already associated with this digest. We are done.
                return;
            }
            newArray[i] = previous[i];

        }
        newArray[i] = sequenceIndex;
        readDigestsToSequenceIndex.put(digest, newArray);

    }

    /**
     * ***********************************************************************
     * Using table lookup
     * Reference: http://snippets.dzone.com/tag/crc32
     * ************************************************************************
     */

    private static final int[] table = {
            0x00000000, 0x77073096, 0xee0e612c, 0x990951ba, 0x076dc419, 0x706af48f, 0xe963a535, 0x9e6495a3,
            0x0edb8832, 0x79dcb8a4, 0xe0d5e91e, 0x97d2d988, 0x09b64c2b, 0x7eb17cbd, 0xe7b82d07, 0x90bf1d91,
            0x1db71064, 0x6ab020f2, 0xf3b97148, 0x84be41de, 0x1adad47d, 0x6ddde4eb, 0xf4d4b551, 0x83d385c7,
            0x136c9856, 0x646ba8c0, 0xfd62f97a, 0x8a65c9ec, 0x14015c4f, 0x63066cd9, 0xfa0f3d63, 0x8d080df5,
            0x3b6e20c8, 0x4c69105e, 0xd56041e4, 0xa2677172, 0x3c03e4d1, 0x4b04d447, 0xd20d85fd, 0xa50ab56b,
            0x35b5a8fa, 0x42b2986c, 0xdbbbc9d6, 0xacbcf940, 0x32d86ce3, 0x45df5c75, 0xdcd60dcf, 0xabd13d59,
            0x26d930ac, 0x51de003a, 0xc8d75180, 0xbfd06116, 0x21b4f4b5, 0x56b3c423, 0xcfba9599, 0xb8bda50f,
            0x2802b89e, 0x5f058808, 0xc60cd9b2, 0xb10be924, 0x2f6f7c87, 0x58684c11, 0xc1611dab, 0xb6662d3d,
            0x76dc4190, 0x01db7106, 0x98d220bc, 0xefd5102a, 0x71b18589, 0x06b6b51f, 0x9fbfe4a5, 0xe8b8d433,
            0x7807c9a2, 0x0f00f934, 0x9609a88e, 0xe10e9818, 0x7f6a0dbb, 0x086d3d2d, 0x91646c97, 0xe6635c01,
            0x6b6b51f4, 0x1c6c6162, 0x856530d8, 0xf262004e, 0x6c0695ed, 0x1b01a57b, 0x8208f4c1, 0xf50fc457,
            0x65b0d9c6, 0x12b7e950, 0x8bbeb8ea, 0xfcb9887c, 0x62dd1ddf, 0x15da2d49, 0x8cd37cf3, 0xfbd44c65,
            0x4db26158, 0x3ab551ce, 0xa3bc0074, 0xd4bb30e2, 0x4adfa541, 0x3dd895d7, 0xa4d1c46d, 0xd3d6f4fb,
            0x4369e96a, 0x346ed9fc, 0xad678846, 0xda60b8d0, 0x44042d73, 0x33031de5, 0xaa0a4c5f, 0xdd0d7cc9,
            0x5005713c, 0x270241aa, 0xbe0b1010, 0xc90c2086, 0x5768b525, 0x206f85b3, 0xb966d409, 0xce61e49f,
            0x5edef90e, 0x29d9c998, 0xb0d09822, 0xc7d7a8b4, 0x59b33d17, 0x2eb40d81, 0xb7bd5c3b, 0xc0ba6cad,
            0xedb88320, 0x9abfb3b6, 0x03b6e20c, 0x74b1d29a, 0xead54739, 0x9dd277af, 0x04db2615, 0x73dc1683,
            0xe3630b12, 0x94643b84, 0x0d6d6a3e, 0x7a6a5aa8, 0xe40ecf0b, 0x9309ff9d, 0x0a00ae27, 0x7d079eb1,
            0xf00f9344, 0x8708a3d2, 0x1e01f268, 0x6906c2fe, 0xf762575d, 0x806567cb, 0x196c3671, 0x6e6b06e7,
            0xfed41b76, 0x89d32be0, 0x10da7a5a, 0x67dd4acc, 0xf9b9df6f, 0x8ebeeff9, 0x17b7be43, 0x60b08ed5,
            0xd6d6a3e8, 0xa1d1937e, 0x38d8c2c4, 0x4fdff252, 0xd1bb67f1, 0xa6bc5767, 0x3fb506dd, 0x48b2364b,
            0xd80d2bda, 0xaf0a1b4c, 0x36034af6, 0x41047a60, 0xdf60efc3, 0xa867df55, 0x316e8eef, 0x4669be79,
            0xcb61b38c, 0xbc66831a, 0x256fd2a0, 0x5268e236, 0xcc0c7795, 0xbb0b4703, 0x220216b9, 0x5505262f,
            0xc5ba3bbe, 0xb2bd0b28, 0x2bb45a92, 0x5cb36a04, 0xc2d7ffa7, 0xb5d0cf31, 0x2cd99e8b, 0x5bdeae1d,
            0x9b64c2b0, 0xec63f226, 0x756aa39c, 0x026d930a, 0x9c0906a9, 0xeb0e363f, 0x72076785, 0x05005713,
            0x95bf4a82, 0xe2b87a14, 0x7bb12bae, 0x0cb61b38, 0x92d28e9b, 0xe5d5be0d, 0x7cdcefb7, 0x0bdbdf21,
            0x86d3d2d4, 0xf1d4e242, 0x68ddb3f8, 0x1fda836e, 0x81be16cd, 0xf6b9265b, 0x6fb077e1, 0x18b74777,
            0x88085ae6, 0xff0f6a70, 0x66063bca, 0x11010b5c, 0x8f659eff, 0xf862ae69, 0x616bffd3, 0x166ccf45,
            0xa00ae278, 0xd70dd2ee, 0x4e048354, 0x3903b3c2, 0xa7672661, 0xd06016f7, 0x4969474d, 0x3e6e77db,
            0xaed16a4a, 0xd9d65adc, 0x40df0b66, 0x37d83bf0, 0xa9bcae53, 0xdebb9ec5, 0x47b2cf7f, 0x30b5ffe9,
            0xbdbdf21c, 0xcabac28a, 0x53b39330, 0x24b4a3a6, 0xbad03605, 0xcdd70693, 0x54de5729, 0x23d967bf,
            0xb3667a2e, 0xc4614ab8, 0x5d681b02, 0x2a6f2b94, 0xb40bbe37, 0xc30c8ea1, 0x5a05df1b, 0x2d02ef8d,
    };

    public final int crc32(final byte[] bytes, final int offset, final int length) {

        int crc = 0xffffffff;
        final int end = offset + length;
        for (int i = offset; i < end; ++i) {
            final byte b = bytes[i];
            crc = (crc >>> 8) ^ table[(crc ^ b) & 0xff];

        }

        // flip bits
        crc = ~crc;
        return crc;
    }

    public final int crc32(final ByteString sequence, final int offset, final int length) {
        int crc = 0xffffffff;
        final int end = offset + length;
        for (int i = offset; i < end; ++i) {
            final byte b = sequence.byteAt(i);
            crc = (crc >>> 8) ^ table[(crc ^ b) & 0xff];

        }

        // flip bits
        crc = ~crc;
        return crc;
    }

    private int crc32ReverseComplement(final byte[] sequence, final int offset, final int length) {

        int crc = 0xffffffff;
        final int end = offset + length - 1;
        for (int i = end; i >= offset; --i) {
            final byte b = complement(sequence[i]);
            crc = (crc >>> 8) ^ table[(crc ^ b) & 0xff];

        }

        // flip bits
        crc = ~crc;
        return crc;
    }

    private int crc32ReverseComplement(final ByteString sequence, final int offset, final int length) {

        int crc = 0xffffffff;
        final int end = offset + length - 1;
        for (int i = end; i >= offset; --i) {
            final byte b = complement(sequence.byteAt(i));
            crc = (crc >>> 8) ^ table[(crc ^ b) & 0xff];

        }

        // flip bits
        crc = ~crc;
        return crc;
    }

    public final int crc32WithMask(final byte[] bytes, final int offset, final int length) {

        int crc = 0xffffffff;
        final int end = offset + length;
        for (int i = offset; i < end; ++i) {
            final byte b = bytes[i];
            if (mask[i]) {
                crc = (crc >>> 8) ^ table[(crc ^ b) & 0xff];

            }
        }

        // flip bits
        crc = ~crc;
        return crc;
    }
    public final int crc32WithMask(final ByteString sequence, final int offset, final int length) {

        int crc = 0xffffffff;
        final int end = offset + length;
        for (int i = offset; i < end; ++i) {
            final byte b =sequence.byteAt(i);
            if (mask[i]) {
                crc = (crc >>> 8) ^ table[(crc ^ b) & 0xff];

            }
        }

        // flip bits
        crc = ~crc;
        return crc;
    }


    /**
     * Return the possible sequence indices for a segment of sequence.  A digest is calculated from the bytes
     * of the sequence between in the range [offset, offset+length[.
     *
     * @param sequence Complete sequence
     * @param offset   Offset from start of complete sequence where the segment starts.
     * @return
     */
    public final int[] lookup(final byte[] sequence, final int offset, final int length) {
        final int digest = digest(sequence, offset, length);
        return readDigestsToSequenceIndex.get(digest);
    }

    /**
     * Get the sequence indices associated with a digest.
     *
     * @param digest
     * @return
     */
    public final int[] getReadIndices(final int digest) {
        return readDigestsToSequenceIndex.get(digest);
    }

    /**
     * Remove a digest from the digest to read index map. getReadIndex and lookup will no longer be able to retrive the
     * sequence index for this digest.
     *
     * @param digest
     */
    public final void remove(final int digest) {
        readDigestsToSequenceIndex.remove(digest);
    }

    /**
     * Remove the association between a digest and a read index. If the digest is associated with no other read index,
     * removes the digest entirely.
     *
     * @param digest
     */
    public void remove(final int digest, final int readIndex) {
        final int[] readIndices = readDigestsToSequenceIndex.get(digest);
        int removed = 0;
        for (int i = 0; i < readIndices.length; ++i) {
            if (readIndices[i] == readIndex) {
                readIndices[i] = -1;
                removed++;
            }
        }
        if (removed == readIndices.length) {
            // this digest no longer references a sequence.
            remove(digest);
        }
    }

    public final boolean confirmMatch(final byte[] readCompressedSeq,
                                      final byte[] byteBuffer,
                                      final int referencePosition,
                                      final int readLength) {
        if (matchesPositiveStrand) {

            final int end = referencePosition + readLength;
            int j = 0;
            for (int i = referencePosition; i < end; ++i) {
                if (readCompressedSeq[j] != byteBuffer[i]) {
                    return false;
                }
                ++j;
            }
            return true;
        } else {
            int j = referencePosition;

            for (int i = this.readLength - 1; i >= 0; --i) {
                if (complement(readCompressedSeq[i]) != byteBuffer[j++]) {
                    return false;
                }
            }
            return true;
        }
    }

    /**
     * Return the complement of base b.
     *
     * @param b
     * @return
     */
    private byte complement(final byte b) {
        switch (b) {
            case 'A':
                return 'T';
            case 'C':
                return 'G';
            case 'T':
                return 'A';
            case 'G':
                return 'C';
            default:
                return b;
        }
    }

    public boolean isMatchingPositiveStrand() {
        return matchesPositiveStrand;
    }


}
