#!/bin/sh

#
# Script to submit a SGE job that will run an alignment and
# concatenate the results using goby
#

# Run the aligner
ALIGN=`qsub -terse goby-align.qsub`
echo $ALIGN

# Concatenate results from the alignment results
CONCAT=`qsub -terse -hold_jid $ALIGN goby-concat.qsub`
echo $CONCAT
