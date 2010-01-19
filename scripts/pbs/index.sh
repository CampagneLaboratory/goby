#!/bin/sh

#
# Script to submit a PBS job that will build a reference index using goby

# Build the reference index
INDEX=`qsub goby-index.qsub`
echo $INDEX
