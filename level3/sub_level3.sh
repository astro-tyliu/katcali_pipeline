#! /bin/bash

file_timestamp=$(date +"%Y%m%d_%H%M%S")
output_dir="/scratch3/users/liutianyang/katcali_pipeline/level3/py_results/cali_${file_timestamp}"
logs_dir="/scratch3/users/liutianyang/katcali_pipeline/level3/logs/cali_${file_timestamp}"
mkdir -p ${output_dir}
mkdir -p ${logs_dir}

# clean_20250306_092316 cali_20250306_163332 1675632179  level2 done
# clean_20250311_143125 cali_20250311_171432 1679247986  level2 done
# clean_20250311_162405 cali_20250312_091713 1679615321  level2 done
# clean_20250311_144810 cali_20250312_091714 1680644082  level2 done

if [ $# -lt 3 ]; then
    echo "Error: 3 arguments are required!"
    echo "Usage: ./script.sh <input_file1> <input_file2> <fname>"
    exit 1
fi

# Assign input parameters to variables
input_file1=$1
input_file2=$2
fname=$3  # 1675632179 1679247986 1679592842 1679615321 1680644082

echo "Input file 1: $input_file1"
echo "Input file 2: $input_file2"
echo "block name: $fname"

# fname="1680644082" 

echo ${fname} ${ant} ${pol}

for i in {000..063}; do
    ant="m${i}"
    for pol in v h; do 

        echo "${fname} ${ant} ${pol}"

        script_name="level3_${fname}_${ant}_${pol}"
        echo "#! /bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12GB
#SBATCH --time=36:00:00
#SBATCH --error=${logs_dir}/job_%J_${ant}.err
#SBATCH --output=${logs_dir}/job_%J_${ant}.out

export SINGULARITY_SHELL=/bin/bash" > ${script_name}

# export SINGULARITY_SHELL=/bin/bash" > sub_${fname}_${ant}_${pol}

        echo "singularity exec /data/exp_soft/containers/katcal.sif ./KATcali_UHF_level3.py ${fname} ${ant} ${pol} ${input_file1} ${input_file2} ${file_timestamp}" >> ${script_name}
        sbatch ${script_name}

        # sbatch sub_${fname}_${ant}_${pol}

        rm -f ${script_name}
        # rm -f sub_${fname}_${ant}_${pol}
        
   done
done
