cd /mnt/raid/milutin/upScaling/Rondonia/

# shp file with tiles:
mySHPDir=/mnt/ssd/milutin/Para_upScaling/ForestChangeMaps/Rondonia/

# get a list of all tile IDs:
ogrinfo -so -ro -al ${mySHPDir}Rondonia_Grid_36km.shp
#
ogrinfo -so -ro -al ${mySHPDir}Rondonia_Grid_36km.shp -dialect sqlite -sql "SELECT id FROM Rondonia_Grid_36km" 
# get the tileID list and save as txt file
ogrinfo -ro -al ${mySHPDir}Rondonia_Grid_36km.shp -dialect sqlite -sql "SELECT id FROM Rondonia_Grid_36km" | grep " id " | awk '{print $4}' > List_of_S1_TileIDs.txt




