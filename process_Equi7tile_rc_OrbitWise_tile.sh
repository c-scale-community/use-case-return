#!/bin/bash

# ----------------------------------------------------------------------------------------------------
# Slum settings
#SBATCH -c 3                                 # one CPU core
#SBATCH --array=0-3599%200                   # how many tasks in the array
#SBATCH --output=/eodc/private/RETURN/upscale_congo/slurm_out/%A_%a.out 
#SBATCH --error=/eodc/private/RETURN/upscale_congo/slurm_out/%A_%a.err 
# ----------------------------------------------------------------------------------------------------
# Get tile from commandline
TILE=$1
echo $TILE
# Get tile from commandline
ODIR=$2
# Start date, host, user
# ----------------------------------------------------------------------------------------------------
#
date
hostname
echo ${USER}
# ----------------------------------------------------------------------------------------------------
# Load software
#
source /home/${USER}/.bashrc
conda activate yeoda4
# ----------------------------------------------------------------------------------------------------
# Set ROW and COL using the SLURM_ARRAY_TASK_ID
# - dividing 3600 tasks into a matrix of 60 rows (ROW) x 60 columns (COL)
#
ROW=$((SLURM_ARRAY_TASK_ID/60+1)) # A = [1-60]
COL=$((SLURM_ARRAY_TASK_ID%60+1)) # B = [1-60]
echo $ROW $COL
# ----------------------------------------------------------------------------------------------------
# Process one folder (specifed Equi7 tile) with chunksize (-s) 250x250
# - a file is 15000 x 15000 and hence an 60 x 60 array with chunksizes of 250 is ok
#

oDirChunk=${ODIR}/${ROW}_${COL}
echo $oDirChunk
mkdir $oDirChunk

echo "START PROCESSING process_Equi7tile_rc_OrbitWise_tile.sh"
echo "srun python /home/fs72029/bschumacher1/use-case-return/process_Equi7tile_OrbitWise.py -p 'DESCENDING' -t $TILE -r $ROW -c $COL -s 250 -o $oDirChunk"

srun python /home/fs72029/bschumacher1/use-case-return/process_Equi7tile_OrbitWise.py -p 'DESCENDING' -t $TILE -r $ROW -c $COL -s 250 -o ${oDirChunk}

sleep 1

# Do the rClone of the output (including the -I flag for overwriting and delete-enebled macaroon):
# rclone -I -v --config=/home/return-mmilenkovic/myMacaroons/maca_return_delete.conf copy $oDirChunk maca_return_delete:disk/S1_SA_UPSCALE/$TILE

sleep 1

#echo ""
#echo "   output copied to dCash storage"
#echo ""
#echo ""

# remove the output directory after the copy is finished:
# rm -rf $oDirChunk


# End date
#
date
echo "done"
