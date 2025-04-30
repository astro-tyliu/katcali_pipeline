#! /bin/bash

file_timestamp=$(date +"%Y%m%d_%H%M%S")
output_dir="/scratch3/users/liutianyang/katcali_pipeline/level4/py_results/level4_${file_timestamp}"
logs_dir="/scratch3/users/liutianyang/katcali_pipeline/level4/logs/cali_${file_timestamp}"
mkdir -p ${output_dir}
mkdir -p ${logs_dir}

# clean_20250306_092316 cali_20250306_163332 cali_20250308_081931 1675632179  level3 done
# clean_20250311_143125 cali_20250311_171432 cali_20250313_140745 1679247986  level3 done
# clean_20250311_162405 cali_20250312_091713 cali_20250314_164547 1679615321  level3 done
# clean_20250311_144810 cali_20250312_091714 cali_20250314_164718 1680644082  level3 done

if [ $# -lt 2 ]; then
    echo "Error: 2 arguments are required!"
    echo "Usage: ./script.sh <input_file3> <fname>"
    exit 1
fi

# Assign input parameters to variables
input_file3=$1
fname=$2  # 1675632179 1679247986 1679592842 1679615321 1680644082

echo "Input file level3: $input_file3"
echo "output directory and block name: level4_${file_timestamp} $fname"

# fname="1680644082" 

echo ${fname} ${ant}

for i in {000..063}; do
    ant="m${i}"

    echo "${fname} ${ant}"

    script_name="level4_${fname}_${ant}"
    echo "#! /bin/bash
    
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16GB
#SBATCH --time=1:00:00
#SBATCH --error=${logs_dir}/job_%J_${ant}.err
#SBATCH --output=${logs_dir}/job_%J_${ant}.out

export SINGULARITY_SHELL=/bin/bash" > ${script_name}

    echo "singularity exec /data/exp_soft/containers/katcal.sif ./KATcali_UHF_level4.py ${fname} ${ant} ${input_file3} ${file_timestamp}" >> ${script_name}
    sbatch ${script_name}

    rm -f ${script_name}
        
done
