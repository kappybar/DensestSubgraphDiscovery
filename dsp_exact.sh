#!/bin/sh


for data in `ls ./data`
do
  if [ !`echo "${data}" | grep weighted` ]; then
    ./src/main "./data/${data}" "dsp-exact"
    echo "\n"
  fi
done




echo DONE