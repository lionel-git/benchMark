#!/bin/bash
TARGET=_defines
g++ -dM -E - < /dev/null               > ${TARGET}.txt
g++               -Q --help=target     > ${TARGET}.default.txt
g++ -march=native -Q --help=target     > ${TARGET}.native.txt

