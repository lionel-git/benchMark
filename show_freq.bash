#!/bin/bash
CPUFREQFILE=scaling_cur_freq
CMD="cat /sys/devices/system/cpu/cpu*/cpufreq/${CPUFREQFILE}"
VCGENCMD=/usr/bin/vcgencmd
if [ -f ${VCGENCMD} ]
then
   CMD="${CMD}; ${VCGENCMD} measure_temp"
fi
#echo $CMD
sudo watch -n 1 ${CMD}

