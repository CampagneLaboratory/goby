#!/bin/sh

#
# Script to submit a SGE job that will run an alignment using goby
#

ALIGN=`qsub -terse -N %SGE_JOB_NAME%-align goby-align.qsub`
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
