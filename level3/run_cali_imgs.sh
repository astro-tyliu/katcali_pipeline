#!/bin/bash

# 1675021905 level3_1675021905_20250609_150000
# 1675210948 level3_1675210948_20250609_150000
# 1675623808 level3_1675623808_20250609_150000
# 1675632179 level3_1675632179_20250609_150000

logs_dir="/scratch3/users/liutianyang/katcali_pipeline/level3/logs"

#SBATCH --job-name=cali_imgs
#SBATCH --cpus-per-task=1
#SBATCH --mem=8GB
#SBATCH --time=2:00:00
#SBATCH --output=${logs_dir}/job_%j.out
#SBATCH --error=${logs_dir}/job_%j.err

if [ "$#" -ne 2 ]; then
    echo "Usage: sbatch $0 <fname> <cali_output>"
    exit 1
fi

fname=$1
cali_output=$2

num_ants=64
ants=()
for i in $(seq 0 $((num_ants-1))); do
    ants+=("m$(printf '%03d' $i)")
done
echo "Generated ants: ${ants[@]}"

pols=("h" "v")

# 循环处理每一个ant和pol组合
for ant in "${ants[@]}"; do
    for pol in "${pols[@]}"; do
        echo "Submitting job for cali_output: $cali_output, fname: $fname, ant: $ant, pol: $pol"
        # sbatch --export=ant=$ant,pol=$pol $0 $fname $cali_output  # 使用 --export 传递参数
        singularity exec /data/exp_soft/containers/katcal.sif python3 ./cali_imgs.py "$cali_output" "$fname" "$ant" "$pol"
    done
done

# # 运行 Python 脚本并传入参数
# singularity exec /data/exp_soft/containers/katcal.sif python3 ./cali_imgs.py "$cali_output" "$fname"
