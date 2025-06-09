#! /bin/bash

file_timestamp=$(date +"%Y%m%d_%H%M%S")

# 1675623808 level3_1675623808_20250605_184237
# 1675643846 level3_1675643846_20250605_184237
# 1676313206 level3_1676313206_20250605_184237
# 1678295187 level3_1678295187_20250605_184309
# 1678743988 level3_1678743988_20250605_184309
# 1675210948 level3_1675210948_20250605_184341

if [ $# -lt 2 ]; then
    echo "Error: 2 arguments are required!"
    echo "Usage: ./script.sh <fname> <input_file3>"
    exit 1
fi

# Assign input parameters to variables
fname=$1  # 1675632179 1679247986 1679592842 1679615321 1680644082
input_file3=$2

echo "output directory and block name: $fname level4_${fname}_${file_timestamp}"
echo "Input file level3: $input_file3"


output_dir="/scratch3/users/liutianyang/katcali_pipeline/level4/py_results/level4_${fname}_${file_timestamp}"
logs_dir="/scratch3/users/liutianyang/katcali_pipeline/level4/logs/job_${fname}_${file_timestamp}"
mkdir -p ${output_dir}
mkdir -p ${logs_dir}

echo ${fname} ${ant}

for i in {000..063}; do
    ant="m${i}"

    echo "${fname} ${ant}"

    script_name="level4_${fname}_${ant}"
    echo "#! /bin/bash
    
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16GB
#SBATCH --time=1:00:00
#SBATCH --error=${logs_dir}/job_${fname}_${ant}${pol}_%J.err
#SBATCH --output=${logs_dir}/job_${fname}_${ant}${pol}_%J.out
#SBATCH --exclude=compute-103

export SINGULARITY_SHELL=/bin/bash" > ${script_name}

    echo "singularity exec /data/exp_soft/containers/katcal.sif ./KATcali_UHF_level4.py ${fname} ${ant} ${input_file3} ${file_timestamp}" >> ${script_name}
    sbatch ${script_name}

    rm -f ${script_name}
        
done
