#!/bin/bash

# 1678122565 level4_1678122565_20250609_150000
# 1676657789 level4_1676657789_20250609_150000
# 1675816512 level4_1675816512_20250609_150000
# 1679247986 level4_1679247986_20250714_050000
# 1675632179 level4_1675632179_20250714_050000
# 1677011008 level4_1677011008_20250714_050000

fname=$1
input_file4=$2

dir="./results/py_results/${input_file4}"

mkdir -p "$dir/filter_final"
# mkdir -p "$dir/Tres"

for i in $(seq -w 000 063); do

    ant="m$i"
    cp "$dir/${fname}_$ant/F_${fname}_${ant}_ch_filter_final.png" "$dir/filter_final/F_${fname}_${ant}_ch_filter_final.png"
    # cp "$dir/${fname}_$ant/F_${fname}_${ant}_Tres.png" "$dir/Tres/F_${fname}_${ant}_Tres.png"

done