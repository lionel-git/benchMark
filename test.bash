#!/bin/bash

CMD=/opt/vc/bin/vcgencmd

#$CMD commands

i="0"

while [ $i -lt 2000 ]
do
    $CMD measure_temp
    $CMD measure_clock arm
    $CMD get_throttled
    sleep 1
    i=$[$i+1]
done
