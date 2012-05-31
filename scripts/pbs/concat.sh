#!/bin/sh

#
# Script to submit a PBS job that will concatenate the results
# of multiple alignments using goby
#

# Concatenate results from the alignment results
CONCAT=`qsub goby-concat.qsub`
echo $CONCAT
