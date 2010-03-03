#!/bin/sh

#
# Script to submit a SGE job that will concatenate the results
# of multiple alignments using goby
#

# Concatenate results from the alignment results
CONCAT=`qsub -terse goby-concat.qsub`
echo $CONCAT
