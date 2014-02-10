#!/bin/bash
revision=$1
destination_directory=$2

# Export the content of the current Git repository at ${revision}  to the  ${destination_directory}
git checkout  ${revision}

git archive ${revision} | tar -x -C ${destination_directory}
