#!/home/return-mmilenkovic/miniconda3/envs/yeoda/bin

"""
    This program takes the name of a Equi7grid tile (T3-level) and calculates 
    time seres features as defined in auxilary_ts_tools_mm/features_from_S1_TS.
    
    Example:
    1. activate the yeoda env.
    2. run this line:
    python process_Equi7tile.py -t 'E078N066T3' -r 1 -c 1 -s 100 -o '/home/return-mmilenkovic/'

"""

import subprocess
import argparse
#
import numpy as np
import os, osr, glob
#
from datetime import datetime
import pandas as pd
import xarray as xr
import rioxarray
# import TUW packages
try:
    from yeoda.products.preprocessed import SIG0DataCube
except:
    sys.exit("Could not import yeoda/SIG0DataCube! Please check if the yeoda enviroment is properly set! \n\n")
    
try:
    from geopathfinder.naming_conventions.yeoda_naming import YeodaFilename
except:
    sys.exit("Could not import geopathfinder/YeodaFilename! Please check if the geopathfinder  is properly installed in the yeoda enviroment! \n\n")
# import my aux functions
try:
    from auxilary_ts_tools_mm import features_from_S1_TS, features_as_xrrray_ufunc
except:
    sys.exit("Could not import auxilary_ts_tools_mm! Please check if the code exsist in one of the pythonpath folders! \n\n")


#######################################################################################################

def process_Equi7tile(iTileName, iRowStep, iColStep, iChunkSize, oDir):
    # check the input:
    if not iTileName:
        print("ERROR: Name of the Equi7grid tile is Not Specified")
        print()
        parser.print_help()
        return
    
    if not iRowStep:
        print("ERROR: Row ID of the tille cunck is Not Specified")
        print()
        parser.print_help()
        return
    
    if not iColStep:
        print("ERROR: Column ID of the tille cunck is Not Specified")
        print()
        parser.print_help()
        return
    
    if not iChunkSize:
        print("ERROR: Cunck size is Not Specified")
        print()
        parser.print_help()
        return
    
    if not oDir:
        print("ERROR: Output directory is Not Specified")
        print()
        parser.print_help()
        return
    
    if not os.path.isdir(oDir):
        os.mkdir(oDir)
        
    #------------------------------------------------------------------------
    # Get the filepaths of the tiff fiels of the Equi7Tile
    #------------------------------------------------------------------------
    tile_dir1_path = r'/project/return/Share/EODC_SA020M/V01R01/' + iTileName
    tile_dir2_path = r'/project/return/Share/EODC_SA020M/V1M0R1/' + iTileName
    # specify other parameters:
    dimensions=['time', 'band', 'extra_field', 'sensor_field']
    #
    filepaths1 = glob.glob(os.path.join(tile_dir1_path,'*.tif'))
    filepaths2 = glob.glob(os.path.join(tile_dir2_path,'*.tif'))
    #------------------------------------------------------------------------
    # Prepare the datacube of the tile
    #------------------------------------------------------------------------
    # initatte:
    sig0_dc1 = SIG0DataCube(filepaths=filepaths1, dimensions=dimensions, filename_class=YeodaFilename, sres=20, continent='SA')
    sig0_dc2 = SIG0DataCube(filepaths=filepaths2, dimensions=dimensions, filename_class=YeodaFilename, sres=20, continent='SA')
    # filter by date:
    toi_start, toi_end = datetime(2017, 1, 1), datetime(2021, 1, 1)
    sig0_dc1 = sig0_dc1.filter_by_dimension([(toi_start, toi_end)], [(">=", "<")], name="time", inplace=True)
    sig0_dc2 = sig0_dc2.filter_by_dimension([(toi_start, toi_end)], [(">=", "<")], name="time", inplace=True)
    # select bands:
    sig0_vv_dc1 = sig0_dc1.filter_by_dimension('VV', name='band')
    sig0_vh_dc1 = sig0_dc1.filter_by_dimension('VH', name='band')
    #
    sig0_vv_dc2 = sig0_dc2.filter_by_dimension('VV', name='band')
    sig0_vh_dc2 = sig0_dc2.filter_by_dimension('VH', name='band')
    # merge and sort the datacubes:
    sig0_vv_dc = sig0_vv_dc1.unite(sig0_vv_dc2)
    sig0_vv_dc = sig0_vv_dc.sort_by_dimension('time', ascending=True)
    #
    sig0_vh_dc = sig0_vh_dc1.unite(sig0_vh_dc2)
    sig0_vh_dc = sig0_vh_dc.sort_by_dimension('time', ascending=True)
    #------------------------------------------------------------------------
    # Load the particular chunck of the datacube
    #------------------------------------------------------------------------
    start_row = iRowStep*iChunkSize
    start_col = iColStep*iChunkSize
    #
    sig0_vv_dc_chunk1 = sig0_vv_dc.load_by_pixels(start_row, start_col, row_size=iChunkSize, col_size=iChunkSize, dtype='xarray')
    sig0_vh_dc_chunk1 = sig0_vv_dc.load_by_pixels(start_row, start_col, row_size=iChunkSize, col_size=iChunkSize, dtype='xarray')
    # remane the xarray variable:
    sig0_vv_dc_chunk1 = sig0_vv_dc_chunk1.rename({'1':'sig0_vv'})
    sig0_vh_dc_chunk1 = sig0_vh_dc_chunk1.rename({'1':'sig0_vh'})
    # rescale the data in 2019 and 2020
    sig0_vv_dc_chunk1['sig0_vv'].loc[slice('2019-1-1','2021-1-1'), :, :] = sig0_vv_dc_chunk1.sel(time=slice('2019-1-1','2021-1-1')).apply(lambda x: np.round(x/10.,1)).sig0_vv.values
    sig0_vh_dc_chunk1['sig0_vh'].loc[slice('2019-1-1','2021-1-1'), :, :] = sig0_vh_dc_chunk1.sel(time=slice('2019-1-1','2021-1-1')).apply(lambda x: np.round(x/10.,1)).sig0_vh.values
    #------------------------------------------------------------------------
    # Apply the time-seres analysis per each x, y location in xarray
    #------------------------------------------------------------------------
    # prepare timestamps
    ts_time_stamps = sig0_vv_dc_chunk1['sig0_vv'][:,0, 0].time.values
    # get features for VV-Band per each pixel
    dist_out_vv = xr.apply_ufunc(features_as_xrrray_ufunc, 
                                 sig0_vv_dc_chunk1['sig0_vv'],
                                 ts_time_stamps,
                                 input_core_dims=[["time"], []],
                                 output_core_dims=[["features"]]
                                 )
    # get features for VH-Band per each pixel
    dist_out_vh = xr.apply_ufunc(features_as_xrrray_ufunc, 
                                 sig0_vh_dc_chunk1['sig0_vh'],
                                 ts_time_stamps,
                                 input_core_dims=[["time"], []],
                                 output_core_dims=[["features"]]
                                 )
    # convert output to dataset
    dist_out_vv_ds = dist_out_vv.to_dataset(dim='features')
    dist_out_vh_ds = dist_out_vh.to_dataset(dim='features')
    # give the proper names to the derived features:
    dist_out_vv_ds = dist_out_vv_ds.rename({0:'exception_label', 1:'ref_mean', 2:'error_margin',
                                            3:'num_of_segments', 4:'TS_end_flag', 5:'TS_end_flag_long', 6:'TS_end_mag',
                                            7:'seg_id', 8:'seg_size', 9:'max_mag', 10:'max_mag_date', 11:'t_pre', 12:'t_post', 13:'t_total',
                                            14:'max_mag_org', 15:'max_mag_org_date', 16:'t_mag_org',
                                            17:'seg2_size'})
    #
    dist_out_vh_ds = dist_out_vh_ds.rename({0:'exception_label', 1:'ref_mean', 2:'error_margin',
                                            3:'num_of_segments', 4:'TS_end_flag', 5:'TS_end_flag_long', 6:'TS_end_mag',
                                            7:'seg_id', 8:'seg_size', 9:'max_mag', 10:'max_mag_date', 11:'t_pre', 12:'t_post', 13:'t_total',
                                            14:'max_mag_org', 15:'max_mag_org_date', 16:'t_mag_org',
                                            17:'seg2_size'})
    #------------------------------------------------------------------------
    # Save the calculated features in a NetCDF file  
    #------------------------------------------------------------------------
    outName_VV = iTileName + '_' + str(iRowStep) + '_' + str(iColStep) + '_' + str(iChunkSize) + '_VV.nc'
    outName_VH = iTileName + '_' + str(iRowStep) + '_' + str(iColStep) + '_' + str(iChunkSize) + '_VH.nc'
    #
    dist_out_vv_ds.to_netcdf(os.path.join(oDir, outName_VV))
    dist_out_vh_ds.to_netcdf(os.path.join(oDir, outName_VH)) 

    
    
#######################################################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser( description=__doc__ , formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-t","--tilename", help="Name (str) of the Equi7grid tile to be processed [mandatory]", dest="iTileName")
    parser.add_argument("-r","--rowstep", type=int, help="Row ID (int) of the tille cunck to be processed [mandatory]", dest="iRowStep")
    parser.add_argument("-c","--colstep", type=int, help="Column ID (int) of the tille cunck to be processed [mandatory]", dest="iColStep")
    parser.add_argument("-s","--chunksize", type=int, help="Cunck size (int) in pixels [mandatory]", dest="iChunkSize")
    parser.add_argument("-o","--odir", help="Path to the output directory [mandatory]", dest="oDir")
    args = parser.parse_args()
    #-------------
    process_Equi7tile(args.iTileName, args.iRowStep, args.iColStep, args.iChunkSize, args.oDir)