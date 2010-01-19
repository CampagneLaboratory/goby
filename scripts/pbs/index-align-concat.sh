#!/bin/sh

#
# Script to submit a PBS job that will build a reference index,
# run an alignment and concatenate the results using goby
#

# Build the reference index
INDEX=`qsub goby-index.qsub`
echo $INDEX

# Run the aligner if the index was built successfully
ALIGN=`qsub -W depend=afterok:$INDEX goby-align.qsub`
echo $ALIGN

# Concatenate results from the alignment results
CONCAT=`qsub -W depend=afterany:$ALIGN goby-concat.qsub`
echo $CONCAT
