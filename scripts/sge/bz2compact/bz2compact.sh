#!/bin/sh

#
# Script to convert a fastq.bz2 file to compact format.
# Example of use:
# foreach file (*.fastq.bz2)
#  bz2compact.sh `pwd`/$file
# end
# The previous lines will convert each file with extension .fastq.bz2 to the compact-reads format, running
# in parallel on an SGE grid. The file bz2compact.sh must be on the path.

#

# Run the aligner

export WORKING_DIR=`dirname "$0"`
export GOBY=`which goby`
export BASENAME=`basename $1 .fastq.bz2`
export DATA_DIR=`dirname $1`

CMD=`qsub -terse -v GOBY=$GOBY -v BASENAME=$BASENAME -v DATA_DIR=$DATA_DIR ${WORKING_DIR}/bz2compact.qsub `
echo $CMD


