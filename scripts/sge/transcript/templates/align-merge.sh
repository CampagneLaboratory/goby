#!/bin/sh

#
# Script to submit a SGE job that will run an alignment using goby
#

ALIGN_JOBS=

for TRANSCRIPT_FILE_NAME in `cat transcript-list.txt`; do
    TRANSCRIPT_NAME=${TRANSCRIPT_FILE_NAME%.*}
    # Run the aligner and pass in the name of the transcript file to run with
    ALIGN=`qsub -terse -v TRANSCRIPT_FILE_NAME=${TRANSCRIPT_FILE_NAME} \
        -N %SGE_JOB_NAME%-${TRANSCRIPT_NAME}-align goby-align.qsub`
    echo $ALIGN
    if [ -z ${ALIGN_JOBS} ]; then
        ALIGN_JOBS=${ALIGN}
    else
        ALIGN_JOBS=${ALIGN_JOBS},${ALIGN}
    fi
done

# Merge results from the alignment
MERGE=`qsub -hold_jid ${ALIGN_JOBS} -terse goby-merge.qsub`
echo $MERGE
