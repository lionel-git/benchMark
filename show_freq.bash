#!/bin/bash
CMD="cat /sys/devices/system/cpu/cpu*/cpufreq/cpuinfo_cur_freq"
VCGENCMD=/usr/bin/vcgencmd
if [ -f ${VCGENCMD} ]
then
   CMD="${CMD}; ${VCGENCMD} measure_temp"
fi
#echo $CMD
sudo watch -n 1 ${CMD}

