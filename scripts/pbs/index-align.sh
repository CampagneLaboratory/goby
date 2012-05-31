#!/bin/sh

#
# Script to submit a PBS job that will build a reference index and
# run an alignment using goby

# Build the reference index
INDEX=`qsub goby-index.qsub`
echo $INDEX

# Run the aligner if the index was built successfully
ALIGN=`qsub -W depend=afterok:$INDEX goby-align.qsub`
echo $ALIGN
