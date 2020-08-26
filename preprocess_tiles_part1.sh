# activate the 'Sentinel1_py37' python enviroment
# conda activate Sentinel1_py37

# set the list with S1 tileId-s:
TileIDsList=/mnt/raid/milutin/upScaling/Rondonia/List_of_S1_TileIDs.txt

# set the folder with python scripts:
pyModuleFolder=/home/milutin/PycharmProjects/tile128/

# loop trough the tiles:
for tileID in $(cat $TileIDsList)
do
    # clock it:
    tStart=$(date +%s)
    # convert the VV, VH and Hansen to xarray DataSet:
    python ${pyModuleFolder}collectionMultiband2xarray_ds.py -i $tileID
    # clock it:
    tEnd=$(date +%s)
    # report
    echo "####################################################"
    echo Tile $tileID converted to xarray DataSet! 
    echo "Processing time [min:sec]:"
    echo $((tEnd - tStart)) | awk '{print int($1/60)":"int($1%60)}'
    echo "####################################################"
done
