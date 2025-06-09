#!/bin/bash

# 1675623808 level3_1675623808_20250605_184237

logs_dir="/scratch3/users/liutianyang/katcali_pipeline/level3/logs"

#SBATCH --job-name=cali_imgs
#SBATCH --cpus-per-task=1
#SBATCH --mem=128GB
#SBATCH --time=2:00:00
#SBATCH --output=${logs_dir}/job_%j.out
#SBATCH --error=${logs_dir}/job_%j.err

# 检查是否提供了两个参数
if [ "$#" -ne 2 ]; then
    echo "Usage: sbatch $0 <fname> <cali_output>"
    exit 1
fi

fname=$1
cali_output=$2

# 运行 Python 脚本并传入参数
singularity exec /data/exp_soft/containers/katcal.sif python3 ./cali_imgs.py "$cali_output" "$fname"
