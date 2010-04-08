#!/bin/sh

#
# Script to submit a SGE job that will run an alignment using goby
#

ALIGN=`qsub -terse -N %SGE_JOB_NAME%-align goby-align.qsub`
echo $ALIGN
