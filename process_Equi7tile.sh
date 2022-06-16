#!/bin/bash

#SBATCH -p interactive -c 1

# Start date, host, user
date
hostname
echo ${USER}

# Load software
source /home/${USER}/.bashrc
conda activate yeoda
echo "HELLO"
srun python process_Equi7tile.py -t 'E078N066T3' -r 1 -c 1 -s 100 -o '/home/return-mmilenkovic/eq7/'

# End date
date
echo "done"
