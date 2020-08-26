# Author:       Milutin Milenkovic
# Created:      23/06/2020
# Copyright:    ???
# Licence:      ???

"""
This python function applies a workflow to process S1 xarray Dataset to calculate the TS categories and TS features.
The inputs are (a) id number of the tile, and (b) the polarisation band. Other parameters are hardcoded.
The outputs are two netCDF file and one HDF file:
(1) first netCDF file is an xarray Dataset with the 6-day regular S1 data cube
(2) second netCDF file an xarray Dataset with TS categories and TS features derived only for the disturbed pixels.
(3) the HDF file is a pandas data frame with the features from the disturbed pixels which TS recovered totally
in the recovery period (the category 2 pixels). The HDF file should be input for the statistical analysis.
"""

# ------------------------------------------------
import sys
import os
import numpy as np
import xarray as xr
import argparse
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from random import randrange
import random
import pandas as pd
import time
from pylab import *

import physt
import seaborn as sns
from pandas.plotting import register_matplotlib_converters
# ------------------------------------------------
# add my tools:
sys.path.append(r'/home/milutin/PycharmProjects')
# ---
from tile128.xarray_and_ts_tools_mm import label_lossyear, \
    ts_feature2ds, plot_TS, plot_2D_map, \
    s1_ts_features_ver2, ts_categoryLabel2ds_ver2, \
    aggregate2reg6DayCube

# ------------------------------------------------
# do some setting from the imported modules
register_matplotlib_converters()
# ------------------------------------------------


def generalized_TS_categories_analysis(tile_id, POLARIZATION):
    """
    This python function applies a workflow to process S1 xarray Dataset to calculate the TS categories and TS features.
    The inputs are (a) id number of the tile, and (b) the polarisation band. Other parameters are hardcoded.
    The outputs are two netCDF file and one HDF file:
    (1) first netCDF file is an xarray Dataset with the 6-day regular S1 data cube
    (2) second netCDF file an xarray Dataset with TS categories and TS features derived only for the disturbed pixels.
    (3) the HDF file is a pandas data frame with the features from the disturbed pixels which TS recovered totally
    in the recovery period (the category 2 pixels). The HDF file should be input for the statistical analysis.
    """
    ####################################################################################################################
    # set variables:
    ####################################################################################################################
    # working directory with the input data:
    # data_dir = r'/mnt/nobackup/milutin/Para_upScaling/Rondonia/test_10_tiles/'
    data_dir = r'/mnt/raid/milutin/upScaling/Rondonia/S1_all_tiles/S1_allDaysDataCubes/'
    # ---------------------------------------------------------------------------
    # the name of the xarray DataSet stored as the netCDF file
    # this file is the output of the 'collectionMultiband2xarray_ds.py' script
    nc_filename = "Rondonia_S1_TileID_" + str(tile_id) + ".nc"
    # output netCDF name:
    nc_6Day_output_name = "Rondonia_ds_6Day_TileID_" + str(tile_id) + ".nc"
    nc_output_name = "Rondonia_dsDisturb_TileID_" + str(tile_id) + "_pol" + str(POLARIZATION) + ".nc"
    # output name for the panda dataFrame:
    df_output_name = "Rondonia_df_2_TileID_" + str(tile_id) + "_pol" + str(POLARIZATION) +".hdf"
    # ---------------------------------------------------------------------------
    # set output folder for 6-day netCDF files:
    out_folder_nc = r'/mnt/raid/milutin/upScaling/Rondonia/S1_all_tiles/S1_6DayDataCubes/'
    # set output folder for disturbed pixels (nc files for disturb pixels, and hdf files):
    out_DistrbFolder = r'/mnt/raid/milutin/upScaling/Rondonia/S1_all_tiles/S1_DisturbedPixelsOut/'
    # ---------------------------------------------------------------------------
    # set the start/end time of the analysis
    # e.g., min. 3 month before/after the time window to deal with the running mean
    start_time = '2016-10-01'
    end_time = '2020-06-01'
    # set the start/end time the history period
    start_history = '2017-01-01'
    end_history = '2018-01-01'
    # ---------------------------------------------------------------------------
    # set the 'lossyear' corresponding to the year of disturbance that is analysed (e.g., value 18 for year 2018)
    LOSS_YEAR = 18
    # set the number of observation used for running mean operator
    # !!! NOTE !!! take care about the time sampling, e.g. 15 is 90 days for the 6-day cube
    sample_num = 15
    # st the threshold [dB] for the std of the running mean TS in the history period
    # std larger than th_std indicates that disturbance happened already in the history period
    th_std = 0.6
    # the label of the disturbed TS that recovered within the monitoring time
    # this label should be set according to the labels introduced in the 'ts_categoryLabel2ds_ver2' function
    LABEL_TS_RECOVERED = 2
    # ---------------------------------------------------------------------------
    os.chdir(data_dir)
    os.getcwd()
    ####################################################################################################################
    # read and prepare the xarray DataSet
    ####################################################################################################################
    ds = xr.open_dataset(nc_filename)
    # select data in the specified time window
    ds = ds.sel(time=slice(start_time, end_time))
    ####################################################################################################################
    # aggregate and resample to 6-day data cube
    ####################################################################################################################
    ds = aggregate2reg6DayCube(ds)
    # save the 6-day cube
    ds.to_netcdf(os.path.join(out_folder_nc, nc_6Day_output_name), mode='w')
    ####################################################################################################################
    # calculate the TS features for all the disturbed pixels
    ####################################################################################################################
    # get the disturbed pixels:
    disturbed_pixels = np.argwhere(ds.LossYear.values == LOSS_YEAR)
    # ----------------------------------------------------------------------
    # prepare an empty list to store all the parameters
    dist_pix_par_all = list()
    # loop trough all the disturbed pixels and get the TS features
    for dist_pix in disturbed_pixels:
        ts_dist_pix = ds[POLARIZATION][:, dist_pix[0], dist_pix[1]]
        dist_pix_par_all.append(s1_ts_features_ver2(ts_dist_pix))

    ####################################################################################################################
    # assign the derived TS features to a new dataset xarray
    ####################################################################################################################
    # define a list with names of TS feature
    # !!! NOTE !!! the order of names follows the order of the output features in the 's1_ts_features_ver2' function
    featNames = ['exception_label', 'ref_mean', 'error_margin', 'num_of_segments', 'TS_end_flag', 'TS_end_flag_long',
                 'TS_end_mag', 'seg_id', 'seg_size', 'max_mag', 'max_mag_date', 't_pre', 't_post', 't_total',
                 'max_mag_org', 'max_mag_org_date', 't_mag_org', 'seg2_size']
    # prepare new xarray DataSet with only disturbed pixels:
    ds_disturb = ds.where(ds.LossYear.values == LOSS_YEAR).mean(dim='time')
    # update 'ds_disturb' by assigning the TS features as additional Data Variables in the xarray DataSet
    for featName in featNames:
        ts_feature2ds(ds_disturb, featName, dist_pix_par_all, disturbed_pixels)

    ####################################################################################################################
    # History period - Std of of the running mean xarray DataStructure
    ####################################################################################################################
    # get the running mean xarray DataSet for the History period,
    ds_disturb_runMean_history = ds.where(ds.LossYear.values == LOSS_YEAR)[POLARIZATION]\
        .rolling(time=sample_num, center=True)\
        .mean()\
        .sel(time=slice(start_history, end_history))
    # -------------
    # calculate the Std of running mean in the History
    # !!! NOTE !!! it reduces the DataSet to a single layer, thus xarray DataArray is the output.
    da_disturb_runMean_history_std = ds_disturb_runMean_history.std(dim='time')
    ####################################################################################################################
    # Categorize disturbed pixels according to its TS shape
    ####################################################################################################################
    # update 'ds_disturb' by assigning the TS category label as additional Coordinates in the xarray DataSet
    # newly assigned coordinate is named 'CategoryLabel'
    ts_categoryLabel2ds_ver2(ds_disturb, da_disturb_runMean_history_std, th_std)
    ####################################################################################################################
    # extra step to save the dataset as netCDF file to avoid above re-running
    ####################################################################################################################
    # save:
    # !!! NOTE !! it is necessary only once to save this data cube as it contains both VV and VH
    if POLARIZATION == 'VH':
        ds_disturb.to_netcdf(os.path.join(out_DistrbFolder, nc_output_name), mode='w')
    ####################################################################################################################
    # convert xarray Dataset into pandas DataFrame structure and select only the TS that fully recovered
    ####################################################################################################################
    df_11 = ds_disturb.to_dataframe()
    # select only TS with the class label 2 (i.e. the parameter LABEL_TS_RECOVERED):
    df_2 = df_11[df_11.CategoryLabel.eq(LABEL_TS_RECOVERED)]
    df_2.reset_index(inplace=True)
    ####################################################################################################################
    # extra step to save the panda dataFrame as netCDF file to avoid above re-running
    ####################################################################################################################
    # save
    df_2.to_hdf(os.path.join(out_DistrbFolder, df_output_name), key='df_2', mode='w')
    return "Tile " + str(tile_id) + " polarisation " + str(POLARIZATION) + "done!"


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--tileID", type=int, dest="tile_id",
                        help="The id of the tile that should be processed [mandatory]")
    parser.add_argument("-p", "--pol", dest="pol",
                        help="S1 polarisation band (e.g. VV, or VH) [mandatory]")
    args = parser.parse_args()
    #
    generalized_TS_categories_analysis(args.tile_id, args.pol)
