#! /bin/bash

# BOX 13: 1750887085 1753129121 1753559424 1753730543 1754679689 1755196048 1755369199 1756059088

# 定义 fname 列表
fname=$1

input_file3="BOX13_20251015_090000"
file_timestamp="BOX13_20251015_090000_"
output_dir="/scratch3/users/liutianyang/katcali_pipeline/level4/py_results/${file_timestamp}"
logs_dir="/scratch3/users/liutianyang/katcali_pipeline/level4/logs/${file_timestamp}"
mkdir -p ${output_dir}
mkdir -p ${logs_dir}

for i in {000..063}; do
    ant="m${i}"

    echo "${fname} ${ant}"

    script_name="level4_${fname}_${ant}"
    echo "#! /bin/bash
    
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16GB
#SBATCH --time=1:00:00
#SBATCH --error=${logs_dir}/job_${fname}_${ant}_%J.err
#SBATCH --output=${logs_dir}/job_${fname}_${ant}_%J.out
#SBATCH --exclude=compute-103

export SINGULARITY_SHELL=/bin/bash" > ${script_name}

    echo "singularity exec /data/exp_soft/containers/katcal.sif ./KATcali_UHF_level4.py ${fname} ${ant} ${input_file3} ${file_timestamp}" >> ${script_name}

    sbatch ${script_name}

    rm -f ${script_name}
done
