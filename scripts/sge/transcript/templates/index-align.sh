#!/bin/sh

#
# Script to submit a SGE job that will build a reference index and
# run an alignment using goby

# Build the reference index
INDEX=`qsub -terse goby-index.qsub`
echo $INDEX

# Run the aligner if the index was built successfully
ALIGN=`qsub -terse -hold_jid $INDEX -v TRANSCRIPT_FILE_NAME=${TRANSCRIPT_FILE_NAME} \
    -N %SGE_JOB_NAME%-${TRANSCRIPT_NAME}-align goby-align.qsub`
echo $ALIGN
