#!/bin/bash
# -------------------------------------
# fetch the data from gDrive
# -------------------------------------
cd /project/return/Share/mm/forAgeLCLU_data
#
rclone copy remote:forAgeLCLU/ . -P

# -------------------------------------
# mosaic rasters
# -------------------------------------
conda activate yeoda
# make a virtual mosaic (supports the usage of the wildcard)
gdalbuildvrt mosaic.vrt *.tif
# check the metadata
gdalinfo mosaic.vrt
# create a geo tiff raster 
gdal_translate -of GTiff -co TILED=YES -co COMPRESS=LZW mosaic.vrt mosaic.tif

