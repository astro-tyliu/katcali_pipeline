#!/bin/bash

input_file5=$1
fname=$2
pix_deg=$3  # 0.3

if [ $# -lt 3 ]; then
    echo "Error: 3 arguments are required!"
    echo "Usage: ./script.sh <input_file5> <fname> <pix_deg>"
    exit 1
fi

dir="./results/py_results/${input_file5}"

mkdir -p "$dir/map"

for i in $(seq -w 000 063); do

    ant="m$i"
    cp "$dir/${fname}_$ant/F_${fname}_${ant}_pix${pix_deg}d_map.png" "$dir/map/F_${fname}_${ant}_pix${pix_deg}d_map.png"

done
