#!/bin/sh

#
# Script to submit a SGE job that will run a differential expression analysis using Goby
#

# Run the aligner

export WORKING_DIR=`dirname "$0"`
export GOBY=`which goby`
export BASENAME=`basename $1 .fastq.bz2`
export DATA_DIR=`dirname $1`

CMD=`qsub -terse -v GOBY=$GOBY -v BASENAME=$BASENAME -v DATA_DIR=$DATA_DIR ${WORKING_DIR}/bz2compact.qsub `
echo $CMD


