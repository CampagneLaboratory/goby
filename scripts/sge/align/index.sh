#!/bin/sh

#
# Script to submit a SGE job that will build a reference index using goby

# Build the reference index
INDEX=`qsub -terse goby-index.qsub`
echo $INDEX
