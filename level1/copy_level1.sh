#!/bin/bash

fname=$1

dir="./results/py_results/20251015_090000"

mkdir -p "$dir/map/$fname"
mkdir -p "$dir/clean_vis/$fname"

for i in $(seq -w 000 063); do
    ant="m$i"
    mv "$dir/${fname}_${ant}_clean_vis.png" "$dir/clean_vis/$fname/${ant}_clean_vis.png"
    mv "$dir/${fname}_${ant}h_map.png" "$dir/map/$fname/${ant}h_map.png"
    mv "$dir/${fname}_${ant}v_map.png" "$dir/map/$fname/${ant}v_map.png"
done