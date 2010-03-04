#!/bin/sh

#
# Script to submit a SGE job that will run an alignment using goby
#

for TRANSCRIPT_FILE_NAME in `cat transcript-list.txt`; do
    TRANSCRIPT_NAME=${TRANSCRIPT_FILE_NAME%.*}
    # Run the aligner and pass in the name of the transcript file to run with
    ALIGN=`qsub -terse -v TRANSCRIPT_FILE_NAME=${TRANSCRIPT_FILE_NAME} \
        -N %SGE_JOB_NAME%-${TRANSCRIPT_FILE}-align goby-align.qsub`
    echo $ALIGN
done
