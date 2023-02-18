#!/bin/bash
TARGET=_defines
g++ -dM -E - < /dev/null               > ${TARGET}.txt
g++               -Q --help=target     > ${TARGET}.default.txt
g++ -march=native -Q --help=target     > ${TARGET}.native.txt

PROC=`uname -p`
if [ $PROC == "riscv64" ]; then
	    g++ -march=rv64imafdcsuh -Q --help=target     > ${TARGET}.riscv.txt
fi


#g++ -march=rv64imafdcsuh -Q --help=target     > ${TARGET}.riscv.txt

