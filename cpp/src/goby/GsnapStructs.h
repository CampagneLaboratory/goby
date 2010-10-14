/*
 * GsnapSequence.h
 *
 *  Created on: Sep 13, 2010
 *      Author: kdorff
 */

#ifndef GSNAPSTRUCTS_H_
#define GSNAPSTRUCTS_H_

#define GSNAP
#undef PMAP

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/**
 * The following definition comes from GMap/GSnap Sequence.c. They are ONLY
 * used by the C++ code as the C code will use the standard definitions in GMap/GSnap.
 */

struct Sequence_T {
    char *acc; /* Accession */
    char *restofheader; /* Rest of header */
    char *contents; /* Original sequence, ends with '\0' */
    char *contents_alloc; /* Allocation */
    int fulllength; /* Full length (not including chopped sequence) */

#ifdef GSNAP
    char *contents_uc; /* Original sequence, ends with '\0' */
    char *contents_uc_alloc; /* Allocation */
    char *chop;
    int choplength;
    char *quality; /* For Illumina short reads read via FASTQ */
    char *quality_alloc; /* Allocation */
#endif

    int trimstart; /* Start of trim */
    int trimend; /* End of trim */
#ifdef PMAP
    int fulllength_given; /* Full length minus implicit stop codon at end */
#endif
    int subseq_offset; /* Used only for subsequences */
    int skiplength; /* Used only for sequences longer than MAXSEQLEN */
};

/* The following types (UINT2, UINT4, UINT8) come from GMap/GSnap types.h. */

/* A 2-byte word */
typedef unsigned short UINT2;

/* A 4-byte word */
typedef unsigned int UINT4;
typedef int INT4;

/* An 8-byte word */
#if (SIZEOF_UNSIGNED_LONG == 8)
#define HAVE_64_BIT
typedef unsigned long UINT8;
#elif (SIZEOF_UNSIGNED_LONG_LONG == 8)
#define HAVE_64_BIT
typedef unsigned long long UINT8;
#endif

#endif /* GSNAPSTRUCTS_H_ */
