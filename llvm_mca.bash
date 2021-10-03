#!/bin/bash
clang++ -O3 benchMark.cpp  -S -o - | llvm-mca -mcpu=btver2 > analysis.txt 2>&1
