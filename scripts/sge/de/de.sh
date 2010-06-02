#!/bin/sh

#
# Script to submit a SGE job that will run a differential expression analysis using Goby
#

# Run the aligner
DE=`qsub -terse goby-de.qsub`
echo $DE
