#!/bin/bash
# ----------------------------------------------------------------------------------------------
# re-project, crop, and resample the Hansen lossyear data to match the S1-10m tile (WGS84 --> UTM - 21S)
# ----------------------------------------------------------------------------------------------
# set the current tile id
# i.e. get the first argument from the shell, e.g. sh prepare_Disturbance_Maps_Rondonia.sh 204 --> get 204
tileID=$1
# set the relevant folders
# S1_DIR="/mnt/nobackup/milutin/Para_upScaling/Rondonia/test_10_tiles/"
S1_DIR="/mnt/raid/milutin/upScaling/Rondonia/S1_all_tiles/GEE_preprocessed/"
# S1_DIR="/mnt/nobackup/milutin/Para_upScaling/ForestChangeMaps/Rondonia/MapBiomas/"
#
DistMap_DIR="/mnt/raid/milutin/upScaling/Rondonia/ForestChangeMaps/Hansen_GFC/"
# set the raster names:
S1_tifName="Rondonia_S1_VH_TileID_"$tileID".tif"
# S1_tifName="mapbiomas-brazil-collection-41-00000000000000000000-2017_2018_utm20S.tif"
#
DistMap_Name="Rondonia_Hansen_GFC-2018-v1.6_lossyear_mosaic"
# change the working dir:
cd $S1_DIR
# get the bbox of the S1 data:
x_LL=$(gdalinfo $S1_tifName | grep "Lower Left" | awk '{print substr($4,1, length($4)-1)}')
y_LL=$(gdalinfo $S1_tifName | grep "Lower Left" | awk '{print substr($5,1, length($4))}')
x_UR=$(gdalinfo $S1_tifName | grep "Upper Right" | awk '{print substr($4,1, length($4)-1)}')
y_UR=$(gdalinfo $S1_tifName | grep "Upper Right" | awk '{print substr($5,1, length($4))}')
# get the epsg code
epsg_code=$(gdalinfo $S1_tifName | grep "EPSG" | awk 'END{print}' | awk -F "," '{print substr($2,1, length($2)-2)}')

# warp the disturbance map to the S1:
cd $DistMap_DIR
gdalwarp -ot Byte -tr 30 30 -r near -co TILED=YES -co COMPRESS=LZW -t_srs EPSG:$epsg_code -te $x_LL $y_LL $x_UR $y_UR -te_srs EPSG:$epsg_code $DistMap_Name".tif" "Rondonia_Hansen_lossyear_TileID_"$tileID".tif"
# move file to the S1 data directory:
mv "Rondonia_Hansen_lossyear_TileID_"$tileID".tif" $S1_DIR





