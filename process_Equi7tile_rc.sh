#!/bin/bash

# Slum settings
#
#SBATCH -p normal                       # slurm partition
#SBATCH -c 3                            # one CPU core
#SBATCH --array=0-224                   # how many tasks in the array
#SBATCH --output=/home/return-mmilenkovic/mm_tools/use-case-return/slurm_out/%A_%a.out

# Start date, host, user
#
date
hostname
echo ${USER}

# Load software
#
source /home/${USER}/.bashrc
conda activate yeoda

# Set ROW and COL using the SLURM_ARRAY_TASK_ID
# - dividing 225 tasks into a matrix of 15 rows (ROW) x 15 columns (COL)
#
ROW=$((SLURM_ARRAY_TASK_ID/15+1)) # A = [1-15]
COL=$((SLURM_ARRAY_TASK_ID%15+1)) # B = [1-15]
#echo $ROW $COL

# Process one folder (E078N066T3) with chunksize (-s) 1000x1000
# - a file is 15000 x 15000 and hence an 15 x 15 array with chunksizes of 1000 is ok
#
echo "START PROCESSING process_Equi7tile.sh"
echo "srun python /home/return-mmilenkovic/mm_tools/use-case-return/process_Equi7tile.py -t 'E078N060T3' -r $ROW -c $COL -s 1000 -o '/project/return/Share/mm/S1_SA_output/E078N060T3/'"

srun python /home/return-mmilenkovic/mm_tools/use-case-return/process_Equi7tile.py -t 'E078N060T3' -r $ROW -c $COL -s 1000 -o '/project/return/Share/mm/S1_SA_output/E078N060T3/'

# End date
#
date
echo "done"
