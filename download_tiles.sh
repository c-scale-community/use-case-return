# set the download folder
downloadFolder=/mnt/raid/milutin/upScaling/Rondonia/S1_all_tiles/GEE_preprocessed/
# set the folder with bash scripts:
bashFolder=/home/milutin/PycharmProjects/tile128/
# set the list with S1 tileId-s:
TileIDsList=/mnt/raid/milutin/upScaling/Rondonia/List_of_S1_TileIDs.txt

for tileID in $(cat $TileIDsList)
do
    # clock it:
    tStart=$(date +%s)
    # download the data from the drive:
    rclone copy remote:GEE_Rondonia_VV/Rondonia_S1_VV_TileID_$tileID.csv $downloadFolder
    rclone copy remote:GEE_Rondonia_VV/Rondonia_S1_VV_TileID_$tileID.tif $downloadFolder
    #
    rclone copy remote:GEE_Rondonia_VH/Rondonia_S1_VH_TileID_$tileID.csv $downloadFolder
    rclone copy remote:GEE_Rondonia_VH/Rondonia_S1_VH_TileID_$tileID.tif $downloadFolder
    #
    # prepare the Hansen data for the tile (convert to the UTM projection and cut to bbox of the S1 tile)
    sh ${bashFolder}prepare_Disturbance_Maps_Rondonia.sh $tileID
    # clock it:
    tEnd=$(date +%s)
    # report
    echo "####################################################"
    echo Tile $tileID downloaded and Hansen data preprocessed! 
    echo "Processing time [min:sec]:"
    echo $((tEnd - tStart)) | awk '{print int($1/60)":"int($1%60)}'
    echo "####################################################"
done
