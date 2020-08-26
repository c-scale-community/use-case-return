# activate the 'Sentinel1_py37' python enviroment
conda activate Sentinel1_py37

# set the list with S1 tileId-s:
TileIDsList=/mnt/raid/milutin/upScaling/Rondonia/List_of_S1_TileIDs_v2.txt

# set the bash file to be run in paralel:
bashFile=/mnt/raid/milutin/upScaling/Rondonia/preprocess_tiles_part2_baseParallel.sh

# run the above bash script in paralel with 4 processors:
cat $TileIDsList | xargs -n 1 -P 1 bash $bashFile



