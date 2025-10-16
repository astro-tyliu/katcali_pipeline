#!/bin/bash

# BOX 13: 1750887085 1753129121 1753559424 1753730543 1754679689 1755196048 1755369199 1756059088

# 生成时间戳
# file_timestamp=$(date +"%Y%m%d_%H%M%S")
file_timestamp="20251015_090000"

fname=$1

# 创建目录（存放本次运行的数据）

OUTPUT_DIR="/scratch3/users/liutianyang/katcali_pipeline/level1/py_results/${file_timestamp}"
logs_dir="/scratch3/users/liutianyang/katcali_pipeline/level1/logs/${file_timestamp}"
mkdir -p ${OUTPUT_DIR}
mkdir -p ${logs_dir}

# '1678743988', '1682448988', '1678295187', '1676313206', '1677020482', '1679592842', '1675643846', '1678734987', '1678122565', '1689003684', '1689176790', '1676657789', '1677777992', '1688399183', '1677183387', '1679247986', '1677011008', '1675623808', '1680626188', '1677002481', '1677195529', '1681143685', '1680798562', '1684087370', '1683492604', '1681920680', '1685641589', '1679419886', '1689090392', '1675210948', '1678726283', '1684781618', '1678467685', '1679605292', '1679615321', '1679333668', '1680644082', '1677174749', '1678381591', '1675021905', '1675632179', '1681229848', '1677795989', '1678899080', '1675816512', '1675106912'

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
    JOB_SCRIPT="KATcali_${fname}_${ant}.sh"

    echo "#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH --time=2:00:00
#SBATCH --error=${logs_dir}/job_${fname}_${ant}.%J.err
#SBATCH --output=${logs_dir}/job_${fname}_${ant}.%J.out
#SBATCH --exclude=compute-103

export SINGULARITY_SHELL=/bin/bash

singularity exec /data/exp_soft/containers/katcal.sif python3 KATcali_level1_py3.py ${fname} ${ant} ${Threshold_factor1} ${Threshold_factor2} ${Threshold_factor11} ${Threshold_factor22} ${OUTPUT_DIR}" > ${JOB_SCRIPT}

    # 提交任务
    sbatch ${JOB_SCRIPT}

    # 删除临时任务脚本
    rm -f ${JOB_SCRIPT}
done