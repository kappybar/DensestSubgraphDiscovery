#!/bin/sh


for data in `ls ./data`
do
  if [ !`echo "${data}" | grep weighted` ]; then
    ./src/main "./data/${data}" "dsp-greedy"
    ./src/main "./data/${data}" "dsp-greedy++"
    echo "\n"
  fi
done




echo DONE