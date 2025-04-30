#!/bin/bash

logs_dir="/scratch3/users/liutianyang/katcali_pipeline/level2/logs"

#SBATCH --job-name=cali_imgs
#SBATCH --cpus-per-task=1
#SBATCH --mem=64GB
#SBATCH --time=1:00:00
#SBATCH --output=${logs_dir}/job_%j.out
#SBATCH --error=${logs_dir}/job_%j.err

# 运行 Python 脚本
singularity exec /data/exp_soft/containers/katcal.sif python3 ./cali_imgs.py
