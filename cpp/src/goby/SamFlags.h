/*
 * Definition of SAM flags, used by C_Alignments when converting from SAM data.
 */

#ifndef SAM_FLAGS_H_
#define SAM_FLAGS_H_

#define SAM_FLAGS_PAIRED_READ        0x0001
#define SAM_FLAGS_PAIRED_MAPPING     0x0002
#define SAM_FLAGS_QUERY_UNMAPPED     0x0004
#define SAM_FLAGS_MATE_UNMAPPED      0x0008
#define SAM_FLAGS_QUERY_MINUSP       0x0010
#define SAM_FLAGS_MATE_MINUSP        0x0020
#define SAM_FLAGS_FIRST_READ_P       0x0040
#define SAM_FLAGS_SECOND_READ_P      0x0080
#define SAM_FLAGS_NOT_PRIMARY        0x0100
#define SAM_FLAGS_BAD_READ_QUALITY   0x0200
#define SAM_FLAGS_DUPLICATE_READ     0x0400

#endif /* SAM_FLAGS_H_ */
