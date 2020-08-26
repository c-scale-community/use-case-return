# activate the 'Sentinel1_py37' python enviroment
conda activate Sentinel1_py37

# set the folder with python scripts:
pyModuleFolder=/home/milutin/PycharmProjects/tile128/

export tileID=$1

time(
echo "####################################################"
echo Start processing tile $tileID
echo "####################################################"
python ${pyModuleFolder}generalized_TS_categories_analysis.py -i $tileID -p VH
echo "--------------------------------------"
echo Start processing tile $tileID VV Band
echo "--------------------------------------"
python ${pyModuleFolder}generalized_TS_categories_analysis.py -i $tileID -p VV
echo "####################################################"
echo End processing tile $tileID
echo "####################################################"
) 2>&1 | tee /mnt/raid/milutin/upScaling/Rondonia/log_files/log_of_tileID_${tileID}.txt
