#!/bin/bash

# clean_20250306_092316 cali_20250306_163332 cali_20250308_081931 1675632179
# clean_20250311_143125 cali_20250311_171432 cali_20250313_140745 1679247986
# clean_20250311_162405 cali_20250312_091713 cali_20250314_164547 1679615321
# clean_20250311_144810 cali_20250312_091714 cali_20250314_164718 1680644082

logs_dir="/scratch3/users/liutianyang/katcali_pipeline/level3/logs"

#SBATCH --job-name=cali_imgs
#SBATCH --cpus-per-task=1
#SBATCH --mem=128GB
#SBATCH --time=2:00:00
#SBATCH --output=${logs_dir}/job_%j.out
#SBATCH --error=${logs_dir}/job_%j.err

# 检查是否提供了两个参数
if [ "$#" -ne 2 ]; then
    echo "Usage: sbatch $0 <cali_output> <fname>"
    exit 1
fi

cali_output=$1
fname=$2

# 运行 Python 脚本并传入参数
singularity exec /data/exp_soft/containers/katcal.sif python3 ./cali_imgs.py "$cali_output" "$fname"
