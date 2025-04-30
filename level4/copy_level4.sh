#!/bin/bash

# level4_20250317_050515 1675632179
# level4_20250317_050616 1679247986
# level4_20250317_050646 1679615321
# level4_20250317_050717 1680644082

input_file4=$1
fname=$2

dir="./results/py_results/${input_file4}"

# mkdir -p "$dir/filter_final"
mkdir -p "$dir/Tres"

for i in $(seq -w 000 063); do

    ant="m$i"
    # cp "$dir/${fname}_$ant/F_${fname}_${ant}_ch_filter_final.png" "$dir/filter_final/F_${fname}_${ant}_ch_filter_final.png"
    cp "$dir/${fname}_$ant/F_${fname}_${ant}_Tres.png" "$dir/Tres/F_${fname}_${ant}_Tres.png"

done