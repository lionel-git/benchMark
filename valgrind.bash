#!/bin/bash
valgrind --tool=callgrind --dump-instr=yes  --collect-jumps=yes ./benchMark Pi
