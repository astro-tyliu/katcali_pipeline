#! /bin/bash

file_timestamp=$(date +"%Y%m%d_%H%M%S")
output_dir="/scratch3/users/liutianyang/katcali_pipeline/level5/py_results/level5_${file_timestamp}"
logs_dir="/scratch3/users/liutianyang/katcali_pipeline/level5/logs/cali_${file_timestamp}"
mkdir -p ${output_dir}
mkdir -p ${logs_dir}

# clean_20250306_092316 cali_20250306_163332 cali_20250308_081931 level4_20250315_092306 1675632179  level4 done
# clean_20250311_143125 cali_20250311_171432 cali_20250313_140745 level4_20250315_080433 1679247986  level4 done
# clean_20250311_162405 cali_20250312_091713 cali_20250314_164547 level4_20250315_170600 1679615321  level4 done
# clean_20250311_144810 cali_20250312_091714 cali_20250314_164718 level4_20250315_170733 1680644082  level4 done

if [ $# -lt 2 ]; then
    echo "Error: 2 arguments are required!"
    echo "Usage: ./script.sh <input_file4> <fname>"
    exit 1
fi

# Assign input parameters to variables
input_file4=$1
fname=$2

echo "Input file level4: $input_file4"
echo "output directory and block name: level5_${file_timestamp} $fname"

echo ${fname} ${ant}

for i in {000..063}; do
    ant="m${i}"

    echo "${fname} ${ant}"

    script_name="level5_${fname}_${ant}"
    echo "#! /bin/bash
    
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16GB
#SBATCH --time=1:00:00
#SBATCH --error=${logs_dir}/job_%J_${ant}.err
#SBATCH --output=${logs_dir}/job_%J_${ant}.out

export SINGULARITY_SHELL=/bin/bash" > ${script_name}

    echo "singularity exec /data/exp_soft/containers/katcal.sif ./KATcali_UHF_level5.py ${fname} ${ant} ${input_file4} ${file_timestamp}" >> ${script_name}
    sbatch ${script_name}

    rm -f ${script_name}
        
done
