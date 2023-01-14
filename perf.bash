#!/bin/bash
EVT=instructions
for ev in `cat perf.evt`
do
    EVT=${EVT},${ev}
done

perf stat -e $EVT ./benchMark Pi
