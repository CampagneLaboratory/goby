@echo off

rem Generate Goby source code from protocol buffer definition files.
rem Assumes that the ".proto" files reside alongside the script.
rem

rem Destination directory for C++ code
set CPP_DEST_DIR=..\cpp

rem Destination directory for Java code
set JAVA_DEST_DIR=..\src

rem Destination directory for Python code
set PYTHON_DEST_DIR=..\python\goby

protoc --cpp_out="%CPP_DEST_DIR%" --java_out="%JAVA_DEST_DIR%" --python_out="%PYTHON_DEST_DIR%" *.proto
