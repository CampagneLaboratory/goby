#!/bin/sh

#
# Script to submit a SGE job that will build a reference index and
# run an alignment using goby

# Build the reference index
INDEX=`qsub -terse goby-index.qsub`
echo $INDEX

# Run the aligner if the index was built successfully
ALIGN=`qsub -terse -hold_jid $INDEX goby-align.qsub`
echo $ALIGN
