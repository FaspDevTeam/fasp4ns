#!/bin/sh
# 
# usage is fhead.sh ${CMAKE_CURRENT_SOURCE_DIR} where we assume that
# ${CMAKE_CURRENT_SOURCE_DIR} has subdirectories src and include and
# the function names from src/*.c are put into include/fasp4ns_functs.h
set +x
/bin/cat $1/src/*.c \
        | awk -v name="fasp4ns_functs.h" -f mkheaders.awk \
        > $1/include/fasp4ns_functs.h
set -x
