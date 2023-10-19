#!/bin/bash
# ------------------------------------------
#
# AUTHOR JBR OONK (SURF)
#
#  * 02/02/2023	:: created version 0.1
#  * 20/02/2023	:: created version 0.2 by Milutin
# ------------------------------------------
#
echo ""
echo "running submit_sbatch_tileWise.sh: " 
date
echo ""

# declaring array list and index iterator
declare -a array=()
i=0

# reading file in row mode, insert each line into array
while IFS= read -r line; do
    array[i]=$line
    let "i++"
    # reading from file path
done < "/home/return-mmilenkovic/mm_tools/use-case-return/equi7tile_list.txt"
#done < "/home/return-mmilenkovic/mm_tools/use-case-return/equi7tile_list_AOI_MatoGrosso.txt"
# done < "/home/return-mmilenkovic/mm_tools/use-case-return/equi7tile_list_AOI_Tapajos.txt" 


# iterate and print the output from array 
for line in "${array[@]}"
  do
    echo "$line"
  done

# create one white space
echo ""

# echo a specific index from array (0 to n_elements-1)
#echo "${array[105]}"


# set start job and number of jobs to submit
starttile=60
finaltile=69


# loop for sbatch submission (runs from 1 to and including ${njobs}
for i in `seq ${starttile} ${finaltile}`;
do

 tile=${array[i]}
 echo $tile
 
 # check if the out directory for the processing tile does not exsisct and create one:
 ODIR=/project/return/Share/mm/S1_SA_UPSCALE/${tile}
 if [ ! -d "$ODIR" ];
 then
	mkdir $ODIR
    mkdir /project/return/Share/mm/S1_SA_UPSCALE/slurm_out
    echo "$ODIR directory has been created."
 fi

 echo "submitting to squeue, job" $i

 # sbatch command
 sbatch /home/return-mmilenkovic/mm_tools/use-case-return/process_Equi7tile_rc_OrbitWise_tile.sh $tile $ODIR 

 echo ""
 echo "   queue check: squeue -u $USER "
 echo "   job   check: scontrol show job <ID>"
 echo ""
 echo ""

 sleep 1

done

