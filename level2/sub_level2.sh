#! /bin/bash

file_timestamp=$(date +"%Y%m%d_%H%M%S")
output_dir="/scratch3/users/liutianyang/katcali_pipeline/level2/py_results/cali_${file_timestamp}"
logs_dir="/scratch3/users/liutianyang/katcali_pipeline/level2/logs/job_${file_timestamp}"
mkdir -p ${output_dir}
mkdir -p ${logs_dir}

# clean_20250306_092316 1675632179 cali_20250306_163332  done
# clean_20250311_143125 1679247986 cali_20250311_171432  done
# clean_20250311_162405 1679615321 cali_20250312_091713  done
# clean_20250311_144810 1680644082 cali_20250312_091714  done
# 1709832691 clean_20250428_092837 cali_20250429_073534
# 1710869782 clean_20250428_101452 cali_20250429_073640
# 1712685146 clean_20250428_101554 cali_20250429_073714
# 1715012489 clean_20250428_101657 cali_20250429_073818

# 定义 fname 列表
# 1679592842 is a bad block
fnames=("1715012489")  # 1675632179 1679247986 1679592842 1679615321 1680644082
input_file="clean_20250428_101657"

# MaxMem: 10GB
# fname 循环
for fname in "${fnames[@]}"; do

    # ant 从 m000 到 m063 循环
    for i in {000..063}; do
        ant="m${i}"

        # pol 循环
        for pol in v h; do 

            echo "${fname} ${ant} ${pol}"

            script_name="level2_${fname}_${ant}_${pol}"
            echo "#! /bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12GB
#SBATCH --time=40:00:00
#SBATCH --error=${logs_dir}/job_${fname}_${ant}_%J.err
#SBATCH --output=${logs_dir}/job_${fname}_${ant}_%J.out

export SINGULARITY_SHELL=/bin/bash" > ${script_name}

            echo "singularity exec /data/exp_soft/containers/katcal.sif python3 ./KATcali_UHF_level2.py ${fname} ${ant} ${pol} ${input_file} ${file_timestamp}" >> ${script_name}

            sbatch ${script_name}

            rm -f ${script_name}

        done
    done
done
