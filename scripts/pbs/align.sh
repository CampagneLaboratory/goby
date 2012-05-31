#!/bin/sh

#
# Script to submit a PBS job that will run an alignment using goby
#

# Run the aligner
ALIGN=`qsub goby-align.qsub`
echo $ALIGN
