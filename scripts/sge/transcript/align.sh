#!/bin/sh

#
# Script to submit a SGE job that will run an alignment using goby
#

# Run the aligner
ALIGN=`qsub -terse goby-align.qsub`
echo $ALIGN
