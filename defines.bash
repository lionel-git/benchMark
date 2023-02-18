#!/bin/bash
TARGET=defines_g++.txt
g++ -dM -E - < /dev/null > ${TARGET}
echo "=========="       >> ${TARGET}
g++ -Q --help=target    >> ${TARGET}
