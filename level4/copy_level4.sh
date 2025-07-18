#!/bin/bash

# level4_1675021905_20250609_150000 1675021905
# level4_1675210948_20250609_150000 1675210948
# level4_1675623808_20250609_150000 1675623808
# level4_1675632179_20250714_050000 1675632179

input_file4=$1
fname=$2

dir="./results/py_results/${input_file4}"

mkdir -p "$dir/filter_final"
# mkdir -p "$dir/Tres"

for i in $(seq -w 000 063); do

    ant="m$i"
    cp "$dir/${fname}_$ant/F_${fname}_${ant}_ch_filter_final.png" "$dir/filter_final/F_${fname}_${ant}_ch_filter_final.png"
    # cp "$dir/${fname}_$ant/F_${fname}_${ant}_Tres.png" "$dir/Tres/F_${fname}_${ant}_Tres.png"

done