#!/bin/sh

#
# Script to submit a SGE job that will build a reference index,
# run an alignment and merge the results using goby#

# Build the reference index
INDEX=`qsub -terse goby-index.qsub`
echo $INDEX

# if the index job is an array job then the id returned
# is of the form N.I-J:K where N is the job id and the
# I,J,K are the array indices and step value
# we only care about the job id (N)
INDEX=${INDEX%%.*}
echo $INDEX

# Run the aligner if the index was built successfully
ALIGN=`qsub -terse -hold_jid $INDEX -N %SGE_JOB_NAME%-align goby-align.qsub`
echo $ALIGN

# if the align job is an array job then the id returned
# is of the form N.I-J:K where N is the job id and the
# I,J,K are the array indices and step value
# we only care about the job id (N)
ALIGN=${ALIGN%%.*}
echo $ALIGN

# Merge results from the alignment
MERGE=`qsub -hold_jid $ALIGN -terse goby-merge.qsub`
echo $MERGE
