#! /bin/bash

# BOX 13: 1750887085 1753129121 1753559424 1753730543 1754679689 1755196048 1755369199 1756059088

# 定义 fname 列表
fname=$1
check_finish=$2  # Only "True" or "False"

input_file="desi1_20251109_160000"
file_timestamp="desi1_20251109_160000"
# file_timestamp=$(date +"%Y%m%d_%H%M%S")
output_dir="/scratch3/users/liutianyang/katcali_pipeline/level2/py_results/${file_timestamp}"
logs_dir="/scratch3/users/liutianyang/katcali_pipeline/level2/logs/${file_timestamp}"
mkdir -p ${output_dir}
mkdir -p ${logs_dir}

# MaxMem: 
# fname 循环
# for fname in "${fnames[@]}"; do

# ant 从 m000 到 m063 循环
for i in {000..063}; do
    ant="m${i}"

    # pol 循环
    for pol in v h; do 

        echo "${fname} ${ant} ${pol}"

        script_name="level2_${fname}_${ant}_${pol}"
        echo "#! /bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=14GB
#SBATCH --time=72:00:00
#SBATCH --error=${logs_dir}/job_${fname}_${ant}${pol}_%J.err
#SBATCH --output=${logs_dir}/job_${fname}_${ant}${pol}_%J.out
#SBATCH --exclude=compute-103

export SINGULARITY_SHELL=/bin/bash" > ${script_name}

        echo "singularity exec /data/exp_soft/containers/katcal.sif python3 ./KATcali_UHF_level2.py ${fname} ${ant} ${pol} ${input_file} ${output_dir} ${check_finish}" >> ${script_name}
        # echo "singularity exec /data/exp_soft/containers/katcal.sif python3 ./tmp/KATcali_UHF_level2.py ${fname} ${ant} ${pol} ${input_file} ${output_dir}" >> ${script_name}

        sbatch ${script_name}

        rm -f ${script_name}

    done
done
# done
