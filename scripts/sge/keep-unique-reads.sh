#!/bin/sh

#
# Script to create a uniqueset file for a compact-reads file.
# Example of use:
# foreach file (full-path/*.compact-reads)
#  bz2compact.sh $file
# end
# The previous lines will analyze each file with extension .compact-reads to create a uniqset (indices of reads that
# are unique and should be kept when removing duplicate reads), running in parallel on an SGE grid.
# The file keep-unique-reads.sh must be on the path.


export WORKING_DIR=`dirname "$0"`
export GOBY=`which goby`
export BASENAME=`basename $1 .compact-reads`
export DATA_DIR=`dirname $1`

CMD=`qsub -terse -v GOBY=$GOBY -v BASENAME=$BASENAME -v DATA_DIR=$DATA_DIR ${WORKING_DIR}/keep-unique-reads.qsub `
echo $CMD


