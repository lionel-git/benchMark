#!/bin/bash
VCGENCMD="sudo vcgencmd"
COMMANDS=`${VCGENCMD} commands`
for command in ${COMMANDS}
do
    command1=${command::-1}
    echo ========== ${command1} ==============
    ${VCGENCMD} ${command1}
done
