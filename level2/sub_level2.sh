#! /bin/bash

# 1676313206 1676313206_20250520_161045
# 1678295187 1678295187_20250520_161045
# 1678743988 1678743988_20250520_161045
# 1682448988 1682448988_20250520_161045

# 定义 fname 列表
# 1679592842 is a bad block
fname="1682448988"  # 1675632179 1679247986 1679592842 1679615321 1680644082
input_file="1682448988_20250520_161045"

file_timestamp=$(date +"%Y%m%d_%H%M%S")
output_dir="/scratch3/users/liutianyang/katcali_pipeline/level2/py_results/level2_${fname}_${file_timestamp}"
logs_dir="/scratch3/users/liutianyang/katcali_pipeline/level2/logs/job_${file_timestamp}"
mkdir -p ${output_dir}
mkdir -p ${logs_dir}

# MaxMem: 10GB
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
#SBATCH --time=48:00:00
#SBATCH --error=${logs_dir}/job_${fname}_${ant}_%J.err
#SBATCH --output=${logs_dir}/job_${fname}_${ant}_%J.out

export SINGULARITY_SHELL=/bin/bash" > ${script_name}

        echo "singularity exec /data/exp_soft/containers/katcal.sif python3 ./KATcali_UHF_level2.py ${fname} ${ant} ${pol} ${input_file} ${output_dir}" >> ${script_name}

        sbatch ${script_name}

        rm -f ${script_name}

    done
done
# done
