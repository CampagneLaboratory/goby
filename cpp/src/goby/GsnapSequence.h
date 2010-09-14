/*
 * GsnapSequence.h
 *
 *  Created on: Sep 13, 2010
 *      Author: kdorff
 */

#ifndef GSNAPSEQUENCE_H_
#define GSNAPSEQUENCE_H_

#define GSNAP
#define PMAP

/**
 * The following definition comes from GMap/GSnap Sequence.c.
 */
#define T Sequence_T
struct T {
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

#endif /* GSNAPSEQUENCE_H_ */
