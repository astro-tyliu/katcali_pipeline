#!/bin/bash

fname=1679592842

dir="./results/py_results/clean_20250311_144739"

#mkdir -p ./results/py_results/scheme1/1684087370/ratio_imgs/  # 创建目标目录（如果不存在）

mkdir -p "$dir/ratio_imgs"
# mkdir -p "$dir/ratio_imgs_x"
mkdir -p "$dir/clean_vis"
mkdir -p "$dir/maps"

for i in $(seq -w 000 063); do
    ant="m$i"
    cp "$dir/${fname}_$ant/masked_ratio.png" "$dir/ratio_imgs/$ant.png"
    # cp "$dir/${fname}_$ant/h_masked_ratio.png" "$dir/ratio_imgs_x/${ant}_h.png"
    # cp "$dir/${fname}_$ant/v_masked_ratio.png" "$dir/ratio_imgs_x/${ant}_v.png"
    cp "$dir/${fname}_$ant/clean_vis.png" "$dir/clean_vis/$ant.png"
    cp "$dir/${fname}_$ant/h_map.png" "$dir/maps/${ant}_h.png"
    cp "$dir/${fname}_$ant/v_map.png" "$dir/maps/${ant}_v.png"

    # cp "$dir/${fname}_$ant/h_map0.png" "$dir/maps0/${ant}_h.png"
    # cp "$dir/${fname}_$ant/v_map0.png" "$dir/maps0/${ant}_v.png"
    
#    fi
done