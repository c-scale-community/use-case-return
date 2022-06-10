#!/bin/bash
#SBATCH -c 1
echo "start script"
hostname
date
conda init bash
source ~/.bashrc
conda activate yeoda
srun python /home/return-roonk/S1MM/test/test_ts_yeoda_processing.py
date
echo "end script"
