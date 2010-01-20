#!/bin/sh

#
# Script to submit a SGE job that will build a reference index,
# run an alignment and concatenate the results using goby
#

# Build the reference index
INDEX=`qsub -terse goby-index.qsub`
echo $INDEX

# Run the aligner if the index was built successfully
ALIGN=`qsub -terse -hold_jid $INDEX goby-align.qsub`
echo $ALIGN

# Concatenate results from the alignment results
CONCAT=`qsub -terse -hold_jid $ALIGN goby-concat.qsub`
echo $CONCAT
