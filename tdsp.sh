#!/bin/sh

for data in `ls ./data`
do
  if [ !`echo "${data}" | grep weighted` ]; then
    ./src/main "./data/${data}" "tdsp-greedy"
    ./src/main "./data/${data}" "tdsp-greedy++"
    echo "\n"
  fi
done

echo DONE