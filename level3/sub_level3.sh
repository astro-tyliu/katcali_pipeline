#! /bin/bash

# file_timestamp=$(date +"%Y%m%d_%H%M%S")
file_timestamp="20250609_150000"

# Assign input parameters to variables
# sed -n '1301,1600p' /scratch3/users/liutianyang/katcali_pipeline/level2/py_results/others/level2_desi1_result_sort.txt | while read line  # 1248
sed -n '1001,1115p' /scratch3/users/liutianyang/katcali_pipeline/level2/py_results/others/level2_desi2_result_sort.txt | while read line

do
    fname=`echo $line | awk '{print $1}'`
    ant=`echo $line | awk '{print $2}'`
    # input_line=$(grep "^${fname} " /scratch3/users/liutianyang/katcali_pipeline/level2/py_results/others/desi1_fname_mapping.txt)
    input_line=$(grep "^${fname} " /scratch3/users/liutianyang/katcali_pipeline/level2/py_results/others/desi2_fname_mapping.txt)
    input_file1=$(echo $input_line | awk '{print $2}')
    input_file2=$(echo $input_line | awk '{print $3}')
    echo "block name: $fname"
    echo "Input file 1: $input_file1"
    echo "Input file 2: $input_file2"
    
    for pol in h v
    do   
        echo ${fname} ${ant} ${pol}
        
        output_dir="/scratch3/users/liutianyang/katcali_pipeline/level3/py_results/level3_${fname}_${file_timestamp}"
        logs_dir="/scratch3/users/liutianyang/katcali_pipeline/level3/logs/job_${fname}_${file_timestamp}"
        mkdir -p ${output_dir}
        mkdir -p ${logs_dir}
        
        echo "${fname} ${ant} ${pol}"
        
        script_name="level3_${fname}_${ant}_${pol}"
        echo "#! /bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12GB
#SBATCH --time=36:00:00
#SBATCH --error=${logs_dir}/job_${fname}_${ant}${pol}_%J.err
#SBATCH --output=${logs_dir}/job_${fname}_${ant}${pol}_%J.out
#SBATCH --exclude=compute-103
    
        export SINGULARITY_SHELL=/bin/bash" > ${script_name}
        
        # export SINGULARITY_SHELL=/bin/bash" > sub_${fname}_${ant}_${pol}
        
        echo "singularity exec /data/exp_soft/containers/katcal.sif ./KATcali_UHF_level3.py ${fname} ${ant} ${pol} ${input_file1} ${input_file2} ${file_timestamp}" >> ${script_name}
        sbatch ${script_name}
        
        # sbatch sub_${fname}_${ant}_${pol}
        
        rm -f ${script_name}
        # rm -f sub_${fname}_${ant}_${pol}
    done
done
