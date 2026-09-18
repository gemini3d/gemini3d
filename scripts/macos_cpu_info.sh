#!/bin/bash
# for macOS only

Nperf=$(sysctl -n hw.nperflevels)
echo "Number of performance levels: $Nperf"

for i in $(seq 0 $((Nperf-1))); do
  N=$(sysctl -n hw.perflevel${i}.physicalcpu)
  Nname=$(sysctl -n hw.perflevel${i}.name)
  echo "$Nname: $N"
done
