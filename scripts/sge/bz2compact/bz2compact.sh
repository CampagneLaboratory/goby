#!/bin/sh

#
# Script to submit a SGE job that will run a differential expression analysis using Goby
#

# Run the aligner
BZ2COMPACT=`qsub -terse bz2compact.qsub`
echo $BZ2COMPACT
