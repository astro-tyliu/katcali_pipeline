#! /bin/bash

file_timestamp=$(date +"%Y%m%d_%H%M%S")

# desi 1: x_cen x_half y_cen y_half; 148 20 -3 11
# 1675623808 level4_1675623808_20250607_120732 148 20 -3 11
# 1675643846 level4_1675643846_20250607_124604 148 20 -3 11
# 1676313206 level4_1676313206_20250607_124632 148 20 -3 11
# 1678295187 level4_1678295187_20250607_124632 148 20 -3 11
# 1678743988 level4_1678743988_20250607_124704 148 20 -3 11
# 1675210948 level4_1675210948_20250607_124705 148 20 -3 11

if [ $# -lt 6 ]; then
    echo "Error: 6 arguments are required!"
    echo "Usage: ./script.sh <fname> <input_file4> <x_cen> <x_half> <y_cen> <y_half>"
    exit 1
fi

# Assign input parameters to variables
fname=$1
input_file4=$2
x_cen=$3
x_half=$4
y_cen=$5
y_half=$6

echo "block name and output directory: $fname level5_${fname}_${file_timestamp}"
echo "Input file level4: $input_file4"

output_dir="/scratch3/users/liutianyang/katcali_pipeline/level5/py_results/level5_${fname}_${file_timestamp}"
logs_dir="/scratch3/users/liutianyang/katcali_pipeline/level5/logs/job_${fname}_${file_timestamp}"
mkdir -p ${output_dir}
mkdir -p ${logs_dir}

echo ${fname} ${ant}

for i in {000..063}; do
    ant="m${i}"

    echo "${fname} ${ant}"

    script_name="level5_${fname}_${ant}"
    echo "#! /bin/bash
    
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16GB
#SBATCH --time=1:00:00
#SBATCH --error=${logs_dir}/job_${fname}_${ant}${pol}_%J.err
#SBATCH --output=${logs_dir}/job_${fname}_${ant}${pol}_%J.out
#SBATCH --exclude=compute-103

export SINGULARITY_SHELL=/bin/bash" > ${script_name}

    echo "singularity exec /data/exp_soft/containers/katcal.sif ./KATcali_UHF_level5.py ${fname} ${ant} ${input_file4} ${file_timestamp} ${x_cen} ${x_half} ${y_cen} ${y_half}" >> ${script_name}
    sbatch ${script_name}

    rm -f ${script_name}
        
done
