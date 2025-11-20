#! /bin/bash

# desi 1: x_cen x_half y_cen y_half; 149 26 -4 12 pix_deg
# desi 2: x_cen x_half y_cen y_half; 168 26 -4 12 pix_deg
# BOX 13: x_cen x_half y_cen y_half; -18 15 -31 15 pix_deg

fname=$1
x_cen=$2
x_half=$3
y_cen=$4
y_half=$5
pix_deg=$6

# manually modifying
file_timestamp_input="BOX13_20251015_090000"
file_timestamp="BOX13_20251015_090000"

output_dir="/scratch3/users/liutianyang/katcali_pipeline/level5/py_results/${file_timestamp_input}"
logs_dir="/scratch3/users/liutianyang/katcali_pipeline/level5/logs/${file_timestamp}"

mkdir -p ${output_dir}
mkdir -p ${logs_dir}

if [ $# -lt 6 ]; then
    echo "Error: 6 arguments are required!"
    echo "Usage: ./script.sh <fname> <x_cen> <x_half> <y_cen> <y_half> <pix_deg>"
    exit 1
fi

for i in {000..063}; do
    ant="m${i}"

    echo "block name and output directory: $fname ${file_timestamp}"
    echo "Input file level4: ${file_timestamp_input}"    
    echo ${fname} ${ant}

    script_name="level5_${fname}_${ant}"
    echo "#! /bin/bash
        
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16GB
#SBATCH --time=1:00:00
#SBATCH --error=${logs_dir}/job_${fname}_${ant}_%J.err
#SBATCH --output=${logs_dir}/job_${fname}_${ant}_%J.out
#SBATCH --exclude=compute-103

export SINGULARITY_SHELL=/bin/bash" > ${script_name}
    
    echo "singularity exec /data/exp_soft/containers/katcal.sif ./KATcali_UHF_level5.py ${fname} ${ant} ${file_timestamp_input} ${file_timestamp} ${x_cen} ${x_half} ${y_cen} ${y_half} ${pix_deg}" >> ${script_name}
    sbatch ${script_name}

    rm -f ${script_name}
            
done
