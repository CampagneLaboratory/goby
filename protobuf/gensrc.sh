#!/bin/bash

#
# Generate Goby source code from protocol buffer definition files.
# Assumes that the ".proto" files reside alongside the script.
#

# Absolute path to this script.
SCRIPT="$(cd "${0%/*}" 2>/dev/null; echo "$PWD"/"${0##*/}")"

# Absolute path this script is in.
SCRIPT_DIR=`dirname $SCRIPT` 

# Destination directory for C++ code
CPP_DEST_DIR=../cpp/src

# Destination directory for Java code
JAVA_DEST_DIR=../src

# Destination directory for Python code
PYTHON_DEST_DIR=../python/goby

pushd ${SCRIPT_DIR} > /dev/null
protoc --cpp_out="${CPP_DEST_DIR}" --java_out="${JAVA_DEST_DIR}" --python_out="${PYTHON_DEST_DIR}" *.proto
popd > /dev/null