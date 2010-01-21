#!/bin/sh

#
# Script to submit a SGE job that will run an alignment and
# concatenate the results using goby
#

# Run the aligner
ALIGN=`qsub -terse goby-align.qsub`
echo $ALIGN

# if the align job is an array job then the id returned
# is of the form N.I-J:K where N is the job id and the
# I,J,K are the array indices and step value
# we only care about the job id (N)
ALIGN=${ALIGN%%.*}
echo $ALIGN

# Concatenate results from the alignment results
CONCAT=`qsub -terse -hold_jid $ALIGN goby-concat.qsub`
echo $CONCAT
