#!/bin/sh

#
# Script to submit a SGE job that will merge the results
# of multiple alignments using goby
#

# Merge results from the alignment
MERGE=`qsub -terse goby-merge.qsub`
echo $MERGE
