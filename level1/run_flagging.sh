#!/bin/bash


function run_function1 {
#    para1=1684087370
    para1=1675643846
    para3=5
    para4=3
    para5=1.05
    para6=1
    para7=1.0
    para8=1.0

    if [ "$para1" -eq 1684087370 ]; then
        bad_antennas=("m053" "m055")
    elif [ "$para1" -eq 1675643846 ]; then
        bad_antennas=("m017" "m020" "m036")
    else
        echo "Invalid value for a"
    fi

    # Generate the range of para2 (from m000 to m063, excluding m055)
    for i in $(seq -w 000 063); do
        para2="m$i"
        if [[ ! " ${bad_antennas[@]} " =~ " $para2 " ]]; then
            # Combine all parameters into an array
            params=($para1 $para2 $para3 $para4 $para5 $para6 $para7 $para8)
            echo "nohup python3 KATcali_level1_py3.py ${params[@]} > /dev/null 2>&1 &"

            # Execute commands in parallel with xargs (up to 50 processes), run in the background without logs
            echo "nohup python3 KATcali_level1_py3.py ${params[@]} > /dev/null 2>&1 &" | xargs -I {} -P 50 bash -c '{}'
        fi
    done
}


# Function 2: Change para3 and para5, with other parameters fixed
function run_function2 {
    para1=1684087370
#    para1=1675643846
    para2=m000
    start_para3=2.5
    end_para3=7
    step_para3=0.5
#    para3=5
#    start_para4=1.8
#    end_para4=4.2
#    step_para4=0.6
    para4=3
    start_para5=0.55
    end_para5=1.45
    step_para5=0.1
#    para5=1.05
#    start_para6=0.6
#    end_para6=1.4
#    step_para6=0.2
    para6=1
    para7='scheme2'

    # Generate parameter ranges using seq and execute with xargs in parallel
    seq $start_para3 $step_para3 $end_para3 | \
        while read para3; do
            seq $start_para5 $step_para5 $end_para5 | \
            while read para5; do
                # Combine all parameters into an array
                params=($para1 $para2 $para3 $para4 $para5 $para6 $para7)

                # Execute commands in parallel with xargs (up to 50 processes), run in the background without logs
                echo "nohup python3 KATcali_level1_py3.py ${params[@]} > /dev/null 2>&1 &" | xargs -I {} -P 50 bash -c '{}'
            done
        done
}


# Function 3: Fix para3 and para5, change para2
function run_function3 {
#    para1 = 1684087370
    para1=1675643846
    para3=5
    para4=3
    para5=1.05
    para6=1
    para7='scheme3'

    if [ "$para1" -eq 1684087370 ]; then
        bad_antennas=("m055")
    elif [ "$para1" -eq 1675643846 ]; then
        bad_antennas=("m017" "m020" "m036")
    else
        echo "Invalid value for a"
    fi

    # Generate the range of para2 (from m000 to m063, excluding m055)
    for i in $(seq -w 000 063); do
        para2="m$i"
        if [[ ! " ${bad_antennas[@]} " =~ " $para2 " ]]; then
            # Combine all parameters into an array
            params=($para1 $para2 $para3 $para4 $para5 $para6 $para7)
            echo "nohup python3 KATcali_level1_py3.py ${params[@]} > /dev/null 2>&1 &"

            # Execute commands in parallel with xargs (up to 50 processes), run in the background without logs
            echo "nohup python3 KATcali_level1_py3.py ${params[@]} > /dev/null 2>&1 &" | xargs -I {} -P 50 bash -c '{}'
        fi
    done
}


# Function 4: Fix para3 and para5, change para2
function run_function4 {
#    para1 = 1684087370
    para1=1675643846
    para3=5
    para4=3
    para5=1.05
    para6=1
    para7='scheme4'
    para8=0.4
    para9=0.5

    if [ "$para1" -eq 1684087370 ]; then
        bad_antennas=("m055")
    elif [ "$para1" -eq 1675643846 ]; then
        bad_antennas=("m017" "m020" "m036")
    else
        echo "Invalid value for a"
    fi

    # Generate the range of para2 (from m000 to m063, excluding m055)
    for i in $(seq -w 000 063); do
        para2="m$i"
        if [[ ! " ${bad_antennas[@]} " =~ " $para2 " ]]; then
            # Combine all parameters into an array
            params=($para1 $para2 $para3 $para4 $para5 $para6 $para7 $para8 $para9)
            echo "nohup python3 KATcali_level1_py3.py ${params[@]} > /dev/null 2>&1 &"

            # Execute commands in parallel with xargs (up to 50 processes), run in the background without logs
            echo "nohup python3 KATcali_level1_py3.py ${params[@]} > /dev/null 2>&1 &" | xargs -I {} -P 50 bash -c '{}'
        fi
    done
}


# Display options and select functionality
echo "Select functionality:"
echo "1) Function 1 (save all results for all antennas)"
echo "2) Function 2 (change para3 and para5 and just save map2.png)"
echo "3) Function 3 (fix para3 and para5, change para2, and just save map2.png)"
read -p "Enter your choice (1, 2 or 3): " choice

# Run the corresponding function based on the user's choice
if [ "$choice" == "1" ]; then
    run_function1
elif [ "$choice" == "2" ]; then
    run_function2
elif [ "$choice" == "3" ]; then
    run_function3
else
    echo "Invalid choice. Exiting."
    exit 1
fi
