#!/bin/bash

# 生成时间戳
file_timestamp=$(date +"%Y%m%d_%H%M%S")

# 创建目录（存放本次运行的数据）
OUTPUT_DIR="/scratch3/users/liutianyang/katcali_pipeline/level1/py_results/clean_${file_timestamp}"
logs_dir="/scratch3/users/liutianyang/katcali_pipeline/level1/logs"
mkdir -p ${OUTPUT_DIR}

fname=$1  # 1709832691 1710869782 1712685146 1715012489 
if [ $# -lt 1 ]; then
    echo "Error: 1 argument is required!"
    echo "Usage: ./sub_level4.sh <fname>"
    exit 1
fi

Threshold_factor1="8."
Threshold_factor2="4."
Threshold_factor11="16."
Threshold_factor22="8."

# 遍历参数二（m000 到 m063）
for i in {000..063}; do
    ant="m${i}"

    # 生成 Slurm 任务脚本
    JOB_SCRIPT="KATcali_${ant}.sh"

    echo "#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH --time=2:00:00
#SBATCH --error=${logs_dir}/job_${fname}_${ant}.%J.err
#SBATCH --output=${logs_dir}/job_${fname}_${ant}.%J.out

export SINGULARITY_SHELL=/bin/bash

singularity exec /data/exp_soft/containers/katcal.sif python3 KATcali_level1_py3.py ${fname} ${ant} ${Threshold_factor1} ${Threshold_factor2} ${Threshold_factor11} ${Threshold_factor22} ${file_timestamp}" > ${JOB_SCRIPT}

    # 提交任务
    sbatch ${JOB_SCRIPT}

    # 删除临时任务脚本
    rm -f ${JOB_SCRIPT}
done