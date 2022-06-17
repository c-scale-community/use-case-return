#!/bin/bash

# Example of running python script in a batch mode

#SBATCH -p interactive                  # slurm partition
#SBATCH -c 1                            # one CPU core
#SBATCH --array=0-224                   # how many tasks in the array

# Get start date and host
date
hostname
echo ${USER}

# Load software and check libs
source /home/${USER}/.bashrc
#conda env list
#conda activate yeoda

ROW=$((SLURM_ARRAY_TASK_ID/15)) # A = [0-14]
COL=$((SLURM_ARRAY_TASK_ID%15)) # B = [0-14]
#echo $ROW $COL

srun python process_Equi7tile.py -t 'E078N066T3' -r $ROW -c $COL -s 1000 -o '/project/return/Share/mm/S1_SA_output/'

# Get end date
date
echo "done"
