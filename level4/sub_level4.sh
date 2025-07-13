#! /bin/bash

# 1675623808 level3_1675623808_20250605_184237
# 1675643846 level3_1675643846_20250605_184237
# 1676313206 level3_1676313206_20250605_184237
# 1678295187 level3_1678295187_20250605_184309
# 1678743988 level3_1678743988_20250605_184309
# 1675210948 level3_1675210948_20250605_184341

file_timestamp="20250609_150000"

# Assign input parameters to variables
sed -n '501,1248p' /scratch3/users/liutianyang/katcali_pipeline/level2/py_results/others/level2_desi1_result_sort.txt | while read line  # counting from 1 instead of 0.

do
    fname=`echo $line | awk '{print $1}'`
    ant=`echo $line | awk '{print $2}'`
    input_file3="level3_${fname}_20250609_150000"
    
    echo "output directory and block name: $fname level4_${fname}_${file_timestamp}"
    echo "Input file level3: $input_file3"
    
    output_dir="/scratch3/users/liutianyang/katcali_pipeline/level4/py_results/sigma40/level4_${fname}_${file_timestamp}"
    logs_dir="/scratch3/users/liutianyang/katcali_pipeline/level4/logs/sigma40/job_${fname}_${file_timestamp}"
    mkdir -p ${output_dir}
    mkdir -p ${logs_dir}
    
    echo ${fname} ${ant}

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
