# Author:       Milutin Milenkovic
# Created:      22/06/2020
# Copyright:    ???
# Licence:      ???

"""
This python function convert multiband S1 rasters (VV and VH) into xarray Dataset and stores it as the netCDF file
The input is the id number of the tile, whereas the other parameters are hardcoded.
"""

import xarray as xr
import os
import pandas as pd
import numpy as np
import argparse


def collectionMultiband2xarray_ds(tile_id):
    """
    This python function convert multiband S1 rasters (VV and VH) into xarray Dataset and stores it as the netCDF file
    The input is the id number of the tile, whereas the other parameters are hardcoded.
    """
    # -----------------------------------------------------------------------
    # set variables:
    # -----------------------------------------------------------------------
    # working directory with the input data:
    # data_dir=r'c:\Workspace\Data\GEE\test2_WUR_GEE_multiband'
    # data_dir = r'/mnt/nobackup/milutin/Para_upScaling/Rondonia/test_10_tiles/'
    data_dir = r'/mnt/raid/milutin/upScaling/Rondonia/S1_all_tiles/GEE_preprocessed/'
    # ---------
    # general processing flags:
    KEEP_ONLY_S1A = False
    LOSSYEAR_MAP = True
    MULTIPLE_DATES = True
    CALC_CR_RVI = False
    # ---------
    # specify the file names of the input rasters (Disturtbance Map, and S1 multiband tiff)
    # gfc_filename = r'fcl_2018_Hensen_roi_tapajos_utm21S.tif'
    gfc_filename = "Rondonia_Hansen_lossyear_TileID_" + str(tile_id) + ".tif"
    filename_vh = os.path.join(data_dir, "Rondonia_S1_VH_TileID_" + str(tile_id) + ".tif")
    filename_vv = os.path.join(data_dir, "Rondonia_S1_VV_TileID_" + str(tile_id) + ".tif")
    # ---------
    # name of the file with datetimes
    filenames_vh = os.path.join(data_dir, "Rondonia_S1_VH_TileID_" + str(tile_id) + ".csv")
    filenames_vv = os.path.join(data_dir, "Rondonia_S1_VV_TileID_" + str(tile_id) + ".csv")
    # ---------
    # specify the output netCDF filename
    out_filename = "Rondonia_S1_TileID_" + str(tile_id) + ".nc"
    # ---------
    # specify the output directory:
    out_dir = r'/mnt/raid/milutin/upScaling/Rondonia/S1_all_tiles/S1_allDaysDataCubes/'
    # -----------------------------------------------------------------------
    # change the working dir:
    # -----------------------------------------------------------------------
    os.chdir(data_dir)
    os.getcwd()
    # -----------------------------------------------------------------------
    # import Hansen data (Global Forest Change - gfc) from GEE:
    # -----------------------------------------------------------------------
    da_gfc = xr.open_rasterio(gfc_filename)
    # check of the input disturbance is Hansen single-, or multi-band raster
    if LOSSYEAR_MAP:
        da_gfc.coords["band"] = ["lossyear"]
    else:
        da_gfc.coords["band"] = [label_1 for label_1 in da_gfc.descriptions]
    # select only the lossyer band
    da_loss = da_gfc.sel(band='lossyear', drop=True)

    # -----------------------------------------------------------------------
    # import S1-VH data:
    # -----------------------------------------------------------------------
    # read the multi-band raster
    da_vh = xr.open_rasterio(filename_vh).rename({'band': 'time'})
    # read the file with datetimes in the file names:
    df_filenames_vh = pd.read_csv(filenames_vh)
    # assign the time stamps:
    da_vh.coords['time'] = [pd.Timestamp(description.split("_")[4][0:8]).to_datetime64() for description in df_filenames_vh.Date.values]

    # either filter out S1B, or keep S1B and sort dates:
    if KEEP_ONLY_S1A:
        ind_S1A = [description.split("_")[0][0:3] == 'S1A' for description in df_filenames_vh.Date.values]
        da_vh_sorted = da_vh[ind_S1A]
    else:
        # sort the dataArray according to time coordinate:
        da_vh_sorted = da_vh.sortby(da_vh.coords.get('time'))

    # -----------------------------------------------------------------------
    # check if there are at all multiple dates in the data:
    # !!! NOTE !!! This is done only for the VH times.
    # The VV times are assumed identical to VH
    # -----------------------------------------------------------------------
    len_all_times = da_vh_sorted.time.values.shape[0]
    len_unique_times = np.unique(da_vh_sorted.time.values).shape[0]
    # if there are no unique dates, set the flag to false:
    if len_all_times == len_unique_times:
        MULTIPLE_DATES = False

    # -----------------------------------------------------------------------

    # find the images with identical dates and mosaic them (with max function)
    # !!! Note !!! this step is time consuming (ca. 30min), consider optimization
    if MULTIPLE_DATES:
        da_vh_grouped = da_vh_sorted.groupby('time').max(dim='time')
    else:
        da_vh_grouped = da_vh_sorted

    # da_vh_grouped.to_netcdf('da_vh_grouped.nc', mode='w')
    # da_vh_grouped = xr.open_dataarray('da_vh_grouped.nc')

    # convert to xarray data structure:
    ds_vh = da_vh_grouped.to_dataset(name="VH")
    # -----------------------------------------------------------------------
    # import S1-VV data:
    # -----------------------------------------------------------------------
    # read the multi-band raster
    da_vv = xr.open_rasterio(filename_vv).rename({'band': 'time'})
    # read the file with datetimes in the file names:
    df_filenames_vv = pd.read_csv(filenames_vv)
    # assign the time stamps:
    da_vv.coords['time'] = [pd.Timestamp(description.split("_")[4][0:8]).to_datetime64() for description in df_filenames_vv.Date.values]

    # either filter out S1B, or keep S1B and sort dates:
    if KEEP_ONLY_S1A:
        ind_S1A = [description.split("_")[0][0:3] == 'S1A' for description in df_filenames_vv.Date.values]
        da_vv_sorted = da_vv[ind_S1A]
    else:
        # sort the dataArray according to time coordinate:
        da_vv_sorted = da_vv.sortby(da_vv.coords.get('time'))

    # find the images with identical dates and mosaic them (with max function)
    # !!! Note !!! this step is time consuming (ca. 30min), consider optimization
    if MULTIPLE_DATES:
        da_vv_grouped = da_vv_sorted.groupby('time').max(dim='time')
    else:
        da_vv_grouped = da_vv_sorted

    # convert to xarray data structure:
    ds_vv = da_vv_grouped.to_dataset(name="VV")
    # -----------------------------------------------------------------------
    # combine the datasets:
    # -----------------------------------------------------------------------
    # as the time dimension is different, the two ds-s are aligned to the intersection of ds-s time stamps (join='inner')
    ds_align = xr.align(ds_vh, ds_vv, join='inner')
    # convert tuple of aligned ds-s ino a single ds:
    ds = xr.merge(ds_align)
    # -----------------------------------------------------------------------
    # calculate extra parameters if specified:
    # -----------------------------------------------------------------------
    if CALC_CR_RVI:
        # calculate the cross ratio (CR)
        ds = ds.assign(CR=ds.VH / ds.VV)
        # calculate the Radar Vegetation Index (RVI)
        ds = ds.assign(RVI=(4*pow(10, ds.VH/10.))/(pow(10, ds.VV/10)+pow(10, ds.VV/10)))

    # -----------------------------------------------------------------------
    # add Hansen lossyear raster as a coordinate in the ds:
    # -----------------------------------------------------------------------
    ds.coords['LossYear'] = (('y', 'x'), da_loss.values)
    # -----------------------------------------------------------------------
    # save the ds as netCDF file
    # -----------------------------------------------------------------------
    ds.to_netcdf(os.path.join(out_dir, out_filename), mode='w')
    return "The processing of Tile " + str(tile_id) + " completed!"


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--tileID", type=int, dest="tile_id",
                        help="The id of the tile that should be processed [mandatory]")
    args = parser.parse_args()
    #
    collectionMultiband2xarray_ds(args.tile_id)

