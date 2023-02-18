#!/bin/bash
TARGET=_defines
g++ -dM -E - < /dev/null           > ${TARGET}.txt         2>&1
g++               -Q --help=target > ${TARGET}.default.txt 2>&1
g++ -march=native -Q --help=target > ${TARGET}.native.txt  2>&1

PROC=`uname -p`
if [ $PROC == "riscv64" ]; then
    g++ -march=rv64imafdcsuh -Q --help=target > ${TARGET}.riscv.txt 2>&1
fi
