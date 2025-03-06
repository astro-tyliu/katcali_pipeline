#!/bin/bash

#fname=1684087370
fname=1675643846

if [ "$para1" -eq 1684087370 ]; then
    bad_antennas=("m053" "m055")
elif [ "$para1" -eq 1675643846 ]; then
    bad_antennas=("m017" "m020" "m036")
else
    echo "Invalid value for a"
fi

dir="/mnt/results/katcali_results/MeerKLASS_UHF/level1/py_results/$fname"

#mkdir -p ./results/py_results/scheme1/1684087370/ratio_imgs/  # 创建目标目录（如果不存在）

for i in $(seq -w 000 063); do
    ant="m$i"
    cp "$dir/$ant/h500300105100_v500300105100_100100_masked_ratio.png" "$dir/ratio_imgs/100100_$ant.png"
    cp "$dir/$ant/h500300105100_100100/masked_ratio.png" "$dir/ratio_imgs_x/100100_${ant}_h.png"
    cp "$dir/$ant/v500300105100_100100/masked_ratio.png" "$dir/ratio_imgs_x/100100_${ant}_v.png"
    cp "$dir/$ant/h500300105100_v500300105100_100100_clean_vis.png" "$dir/clean_vis/100100_$ant.png"
    cp "$dir/$ant/h500300105100_100100/map.png" "$dir/maps/100100_${ant}_h.png"
    cp "$dir/$ant/v500300105100_100100/map.png" "$dir/maps/100100_${ant}_v.png"
#    fi
done