#! /bin/bash

# file_timestamp=$(date +"%Y%m%d_%H%M%S")

# desi 1: x_cen x_half y_cen y_half; 149 26 -4 12
# 1675623808 level4_1675623808_20250607_120732 149 26 -4 12
# 1675643846 level4_1675643846_20250607_124604 149 26 -4 12
# 1676313206 level4_1676313206_20250607_124632 149 26 -4 12
# 1678295187 level4_1678295187_20250607_124632 149 26 -4 12
# 1678743988 level4_1678743988_20250607_124704 149 26 -4 12
# 1675210948 level4_1675210948_20250607_124705 149 26 -4 12

file_timestamp="20250609_150000"

if [ $# -lt 4 ]; then
    echo "Error: 4 arguments are required!"
    echo "Usage: ./script.sh <x_cen> <x_half> <y_cen> <y_half>"
    exit 1
fi

sed -n '501,1248p' /scratch3/users/liutianyang/katcali_pipeline/level2/py_results/others/level2_desi1_result_sort.txt | while read line

# Assign input parameters to variables
do
    fname=`echo $line | awk '{print $1}'`
    ant=`echo $line | awk '{print $2}'`
    input_file4="level4_${fname}_20250609_150000"

    # fname=$1
    # input_file4=$2
    x_cen=$1
    x_half=$2
    y_cen=$3
    y_half=$4
    
    echo "block name and output directory: $fname level5_${fname}_${file_timestamp}"
    echo "Input file level4: $input_file4"
    
    output_dir="/scratch3/users/liutianyang/katcali_pipeline/level5/py_results/sigma40/level5_${fname}_${file_timestamp}"
    logs_dir="/scratch3/users/liutianyang/katcali_pipeline/level5/logs/sigma40/job_${fname}_${file_timestamp}"
    mkdir -p ${output_dir}
    mkdir -p ${logs_dir}
    
    echo ${fname} ${ant}

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
