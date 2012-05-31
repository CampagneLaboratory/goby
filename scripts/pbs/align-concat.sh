#!/bin/sh

#
# Script to submit a PBS job that will run an alignment and
# concatenate the results using goby
#

# Run the aligner
ALIGN=`qsub goby-align.qsub`
echo $ALIGN

# Concatenate results from the alignment results
CONCAT=`qsub -W depend=afterany:$ALIGN goby-concat.qsub`
echo $CONCAT
