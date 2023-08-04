#!/home/return-mmilenkovic/miniconda3/envs/yeoda/bin
#!/home/return-mmilenkovic/miniconda3/envs/yeoda/bin

"""
    This program takes the name of a Equi7grid tile (T3-level) and calculates 
    time seres features as defined in auxilary_ts_tools_mm/features_from_S1_TS.
    
    Example:
    1. activate the yeoda env.
    2. run this line:
    python process_Equi7tile.py -t 'E078N066T3' -r 1 -c 1 -s 100 -o '/path/to/your/output/folder/'

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

def process_Equi7tile(iOrbitPass, iTileName, iRowStep, iColStep, iChunkSize, oDir):
    # check the input:
    if not iOrbitPass:
        print("ERROR: Orbit Pass is is Not Specified")
        print()
        parser.print_help()
        return
    
    if not iOrbitPass.upper() in ['ASCENDING', 'DESCENDING', 'ALL']:
        print("ERROR: Orbit Pass is not corectly specified, either: ASCENDING, DESCENDING, or ALL")
        print()
        parser.print_help()
        return
    
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
        try:
             os.mkdir(oDir)
        except:
            print("The folder is alredy exsiting.")
       
    #------------------------------------------------------------------------
    # correct the row, col indexing to start from 0,0 instead of 1,1 
    #------------------------------------------------------------------------
    iRowStep = iRowStep - 1
    iColStep = iColStep - 1
    #------------------------------------------------------------------------
    # check if the chink is inside the AOI, i.e., in the tile, row, col list 
    #------------------------------------------------------------------------
    # path to the csv with list with tile, col, and row of chanks within AOI:
    my_csv_chanks = r'/home/return-mmilenkovic/mm_tools/use-case-return/List_tile_col_row_250chanks.txt.csv'
    # read it as pandas dataFrame:
    tcr_chanks = pd.read_csv(my_csv_chanks)
    # handle chanks outside the AOI, i.e. if the chank is not in the list:
    if not ((tcr_chanks['TILE'] == iTileName) & (tcr_chanks['COL'] == iColStep) & (tcr_chanks['ROW'] == iRowStep)).any():
        print("The chunk with row {myRow} and column {myCol} in tile {myTile} is outside the Amazon basin."
              .format(myRow=str(iRowStep), myCol=str(iColStep), myTile=iTileName))
        print()
        return
    
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
    # Prepare selections for 'ASCENDING' and 'DESCENDING' orbit passes: 
    #------------------------------------------------------------------------
    # get the unique list of ascending and descending orbits:
    desc_list = [aa for aa in sig0_vv_dc.inventory.extra_field.unique().tolist() if aa[0]=='D']
    asce_list = [aa for aa in sig0_vv_dc.inventory.extra_field.unique().tolist() if aa[0]=='A']
    # define list with orbit names
    if iOrbitPass.upper() == 'ALL':
        my_orbitNames = desc_list + asce_list
    elif iOrbitPass.upper() == 'DESCENDING':
        my_orbitNames = desc_list
    elif iOrbitPass.upper() == 'ASCENDING':
        my_orbitNames = asce_list

    # loop trough each orbit and create a list of filters
    sig0_vv_dc_list = []
    sig0_vh_dc_list = []
    for orbitName in my_orbitNames:
        #
        my_OrbitFilter_vv = sig0_vv_dc.filter_by_dimension(orbitName, name="extra_field")
        my_OrbitFilter_vv = my_OrbitFilter_vv.sort_by_dimension('time', ascending=True)
        sig0_vv_dc_list.append(my_OrbitFilter_vv)
        #
        my_OrbitFilter_vh = sig0_vh_dc.filter_by_dimension(orbitName, name="extra_field")
        my_OrbitFilter_vh = my_OrbitFilter_vh.sort_by_dimension('time', ascending=True)
        sig0_vh_dc_list.append(my_OrbitFilter_vh)
        
    #------------------------------------------------------------------------
    # Print all orbits
    #------------------------------------------------------------------------
    print('xxxxxxxx')
    print('Following orbits are listed: {}'.format(my_orbitNames))
    print('xxxxxxxx')
    #------------------------------------------------------------------------
    # Process orbit-wise
    #------------------------------------------------------------------------
    # set starting chunk row and column
    start_row = iRowStep*iChunkSize
    start_col = iColStep*iChunkSize
    # load data according to specified orbit pass:
    for orbitName, sig0_vv_dc, sig0_vh_dc in zip(my_orbitNames, sig0_vv_dc_list, sig0_vh_dc_list):
        # check if orbit has images: 
        if sig0_vv_dc.inventory.shape[0] == 0:
            print('xxxxxxxx')
            print('The orbit {} does not contain images'.format(orbitName))
            print('xxxxxxxx')
            print()
            continue
        #------------------------------------------------------------------------
        # Print all orbits
        #------------------------------------------------------------------------
        print()
        print('xxxxxxxx')
        print('Processing the orbit: {}'.format(orbitName))
        print('xxxxxxxx')
        #------------------------------------------------------------------------
        # Load the particular chunck of the datacube
        #------------------------------------------------------------------------
        sig0_vv_dc_chunk1 = sig0_vv_dc.load_by_pixels(start_row, start_col, row_size=iChunkSize, col_size=iChunkSize, dtype='xarray')
        sig0_vh_dc_chunk1 = sig0_vh_dc.load_by_pixels(start_row, start_col, row_size=iChunkSize, col_size=iChunkSize, dtype='xarray') 
        # rename the xarray variable:
        sig0_vv_dc_chunk1 = sig0_vv_dc_chunk1.rename({'1':'sig0_vv'})
        sig0_vh_dc_chunk1 = sig0_vh_dc_chunk1.rename({'1':'sig0_vh'})
        # rescale the data in 2019 and 2020
        sig0_vv_dc_chunk1['sig0_vv'].loc[slice('2019-1-1','2021-1-1'), :, :] = sig0_vv_dc_chunk1.sel(time=slice('2019-1-1','2021-1-1')).apply(lambda x: np.round(x/10.,1)).sig0_vv.values
        sig0_vh_dc_chunk1['sig0_vh'].loc[slice('2019-1-1','2021-1-1'), :, :] = sig0_vh_dc_chunk1.sel(time=slice('2019-1-1','2021-1-1')).apply(lambda x: np.round(x/10.,1)).sig0_vh.values
        #------------------------------------------------------------------------
        # deal with chunk size containing all nan values
        #------------------------------------------------------------------------
        if (sig0_vv_dc_chunk1['sig0_vv'].isnull().values).all() | (sig0_vh_dc_chunk1['sig0_vh'].isnull().values).all():
            print('The chunk with row {myRow} and column {myCol} in the orbit {myOrbName} contains only NaN values.'
                 .format(myRow=str(iRowStep), myCol=str(iColStep), myOrbName=orbitName))
            continue
        #------------------------------------------------------------------------
        # Apply the time-seres analysis per each x, y location in xarray
        #------------------------------------------------------------------------
        # prepare timestamps
        ts_time_stamps = sig0_vv_dc_chunk1['sig0_vv'][:,0, 0].time.values
        ts_time_stamps_vh = sig0_vh_dc_chunk1['sig0_vh'][:,0, 0].time.values
        #
        print('Lenght of vv timestamps: {}'.format(len(ts_time_stamps)))
        print('Lenght of vh timestamps: {}'.format(len(ts_time_stamps_vh)))
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
                                     ts_time_stamps_vh,
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
        # Calculate total 4-year statistics per time-series/pixel   
        #------------------------------------------------------------------------
        # for VV
        dist_out_vv_ds['count_total'] = sig0_vv_dc_chunk1['sig0_vv'].count('time')
        dist_out_vv_ds['min_total'] = sig0_vv_dc_chunk1['sig0_vv'].min('time')
        dist_out_vv_ds['max_total'] = sig0_vv_dc_chunk1['sig0_vv'].max('time')
        dist_out_vv_ds['mean_total'] = sig0_vv_dc_chunk1['sig0_vv'].mean('time')
        dist_out_vv_ds['median_total'] = sig0_vv_dc_chunk1['sig0_vv'].median('time')
        dist_out_vv_ds['std_total'] = sig0_vv_dc_chunk1['sig0_vv'].std('time')
        dist_out_vv_ds['q10_total'] = sig0_vv_dc_chunk1['sig0_vv'].quantile(0.1, dim='time').drop('quantile', dim=None)
        dist_out_vv_ds['q90_total'] = sig0_vv_dc_chunk1['sig0_vv'].quantile(0.9, dim='time').drop('quantile', dim=None)
        dist_out_vv_ds['q25_total'] = sig0_vv_dc_chunk1['sig0_vv'].quantile(0.25, dim='time').drop('quantile', dim=None)
        dist_out_vv_ds['q75_total'] = sig0_vv_dc_chunk1['sig0_vv'].quantile(0.75, dim='time').drop('quantile', dim=None)
        # get the absolute deviation/residual (from median)
        absMedRes_vv_da = xr.apply_ufunc(lambda ts:  abs(ts - np.nanmedian(ts)), 
                                         sig0_vv_dc_chunk1['sig0_vv'],
                                         input_core_dims=[["time"]],
                                         output_core_dims=[["time"]]
                                        )
        # get MAD (median absolute deviation)
        dist_out_vv_ds['MAD_total'] = absMedRes_vv_da.median('time')
        #------------------------------------------------------------
        # for VH
        dist_out_vh_ds['count_total'] = sig0_vh_dc_chunk1['sig0_vh'].count('time')
        dist_out_vh_ds['min_total'] = sig0_vh_dc_chunk1['sig0_vh'].min('time')
        dist_out_vh_ds['max_total'] = sig0_vh_dc_chunk1['sig0_vh'].max('time')
        dist_out_vh_ds['mean_total'] = sig0_vh_dc_chunk1['sig0_vh'].mean('time')
        dist_out_vh_ds['median_total'] = sig0_vh_dc_chunk1['sig0_vh'].median('time')
        dist_out_vh_ds['std_total'] = sig0_vh_dc_chunk1['sig0_vh'].std('time')
        dist_out_vh_ds['q10_total'] = sig0_vh_dc_chunk1['sig0_vh'].quantile(0.1, dim='time').drop('quantile', dim=None)
        dist_out_vh_ds['q90_total'] = sig0_vh_dc_chunk1['sig0_vh'].quantile(0.9, dim='time').drop('quantile', dim=None)
        dist_out_vh_ds['q25_total'] = sig0_vh_dc_chunk1['sig0_vh'].quantile(0.25, dim='time').drop('quantile', dim=None)
        dist_out_vh_ds['q75_total'] = sig0_vh_dc_chunk1['sig0_vh'].quantile(0.75, dim='time').drop('quantile', dim=None)
        # get the absolute deviation/residual (from median)
        absMedRes_vh_da = xr.apply_ufunc(lambda ts:  abs(ts - np.nanmedian(ts)), 
                                         sig0_vh_dc_chunk1['sig0_vh'],
                                         input_core_dims=[["time"]],
                                         output_core_dims=[["time"]]
                                        )
        # get MAD (median absolute deviation)
        dist_out_vh_ds['MAD_total'] = absMedRes_vh_da.median('time')
        #------------------------------------------------------------------------
        # Calculate anual statistics per time-series/pixel   
        #------------------------------------------------------------------------
        # get the xarray groupBy object with yearly data
        annual_vv = sig0_vv_dc_chunk1['sig0_vv'].groupby("time.year") 
        annual_vh = sig0_vh_dc_chunk1['sig0_vh'].groupby("time.year")
        #------------------------------------------------------------
        # for VV
        dist_out_vv_ds['count_2017'] = annual_vv.count('time').sel(year=2017, drop=True)
        dist_out_vv_ds['count_2018'] = annual_vv.count('time').sel(year=2018, drop=True)
        dist_out_vv_ds['count_2019'] = annual_vv.count('time').sel(year=2019, drop=True)
        dist_out_vv_ds['count_2020'] = annual_vv.count('time').sel(year=2020, drop=True)
        #
        dist_out_vv_ds['max_2017'] = annual_vv.max('time').sel(year=2017, drop=True)
        dist_out_vv_ds['max_2018'] = annual_vv.max('time').sel(year=2018, drop=True)
        dist_out_vv_ds['max_2019'] = annual_vv.max('time').sel(year=2019, drop=True)
        dist_out_vv_ds['max_2020'] = annual_vv.max('time').sel(year=2020, drop=True)
        #
        dist_out_vv_ds['min_2017'] = annual_vv.min('time').sel(year=2017, drop=True)
        dist_out_vv_ds['min_2018'] = annual_vv.min('time').sel(year=2018, drop=True)
        dist_out_vv_ds['min_2019'] = annual_vv.min('time').sel(year=2019, drop=True)
        dist_out_vv_ds['min_2020'] = annual_vv.min('time').sel(year=2020, drop=True)
        #
        dist_out_vv_ds['mean_2017'] = annual_vv.mean('time').sel(year=2017, drop=True)
        dist_out_vv_ds['mean_2018'] = annual_vv.mean('time').sel(year=2018, drop=True)
        dist_out_vv_ds['mean_2019'] = annual_vv.mean('time').sel(year=2019, drop=True)
        dist_out_vv_ds['mean_2020'] = annual_vv.mean('time').sel(year=2020, drop=True)
        #
        dist_out_vv_ds['median_2017'] = annual_vv.median('time').sel(year=2017, drop=True)
        dist_out_vv_ds['median_2018'] = annual_vv.median('time').sel(year=2018, drop=True)
        dist_out_vv_ds['median_2019'] = annual_vv.median('time').sel(year=2019, drop=True)
        dist_out_vv_ds['median_2020'] = annual_vv.median('time').sel(year=2020, drop=True)
        #
        dist_out_vv_ds['std_2017'] = annual_vv.std('time').sel(year=2017, drop=True)
        dist_out_vv_ds['std_2018'] = annual_vv.std('time').sel(year=2018, drop=True)
        dist_out_vv_ds['std_2019'] = annual_vv.std('time').sel(year=2019, drop=True)
        dist_out_vv_ds['std_2020'] = annual_vv.std('time').sel(year=2020, drop=True)
        #
        dist_out_vv_ds['q10_2017'] = annual_vv.quantile(0.1, 'time').sel(year=2017, drop=True).drop('quantile', dim=None)
        dist_out_vv_ds['q10_2018'] = annual_vv.quantile(0.1, 'time').sel(year=2018, drop=True).drop('quantile', dim=None)
        dist_out_vv_ds['q10_2019'] = annual_vv.quantile(0.1, 'time').sel(year=2019, drop=True).drop('quantile', dim=None)
        dist_out_vv_ds['q10_2020'] = annual_vv.quantile(0.1, 'time').sel(year=2020, drop=True).drop('quantile', dim=None)
        #
        dist_out_vv_ds['q90_2017'] = annual_vv.quantile(0.9, 'time').sel(year=2017, drop=True).drop('quantile', dim=None)
        dist_out_vv_ds['q90_2018'] = annual_vv.quantile(0.9, 'time').sel(year=2018, drop=True).drop('quantile', dim=None)
        dist_out_vv_ds['q90_2019'] = annual_vv.quantile(0.9, 'time').sel(year=2019, drop=True).drop('quantile', dim=None)
        dist_out_vv_ds['q90_2020'] = annual_vv.quantile(0.9, 'time').sel(year=2020, drop=True).drop('quantile', dim=None)
        #
        dist_out_vv_ds['q25_2017'] = annual_vv.quantile(0.25, 'time').sel(year=2017, drop=True).drop('quantile', dim=None)
        dist_out_vv_ds['q25_2018'] = annual_vv.quantile(0.25, 'time').sel(year=2018, drop=True).drop('quantile', dim=None)
        dist_out_vv_ds['q25_2019'] = annual_vv.quantile(0.25, 'time').sel(year=2019, drop=True).drop('quantile', dim=None)
        dist_out_vv_ds['q25_2020'] = annual_vv.quantile(0.25, 'time').sel(year=2020, drop=True).drop('quantile', dim=None)
        #
        dist_out_vv_ds['q75_2017'] = annual_vv.quantile(0.75, 'time').sel(year=2017, drop=True).drop('quantile', dim=None)
        dist_out_vv_ds['q75_2018'] = annual_vv.quantile(0.75, 'time').sel(year=2018, drop=True).drop('quantile', dim=None)
        dist_out_vv_ds['q75_2019'] = annual_vv.quantile(0.75, 'time').sel(year=2019, drop=True).drop('quantile', dim=None)
        dist_out_vv_ds['q75_2020'] = annual_vv.quantile(0.75, 'time').sel(year=2020, drop=True).drop('quantile', dim=None)
        #
        MAD_vv = annual_vv.map(lambda my_xr:  abs(my_xr - my_xr.median('time')).median('time'))
        #
        dist_out_vv_ds['MAD_2017'] = MAD_vv.sel(year=2017, drop=True)
        dist_out_vv_ds['MAD_2018'] = MAD_vv.sel(year=2018, drop=True)
        dist_out_vv_ds['MAD_2019'] = MAD_vv.sel(year=2019, drop=True)
        dist_out_vv_ds['MAD_2020'] = MAD_vv.sel(year=2020, drop=True)
        #------------------------------------------------------------
        # for VH
        dist_out_vh_ds['count_2017'] = annual_vh.count('time').sel(year=2017, drop=True)
        dist_out_vh_ds['count_2018'] = annual_vh.count('time').sel(year=2018, drop=True)
        dist_out_vh_ds['count_2019'] = annual_vh.count('time').sel(year=2019, drop=True)
        dist_out_vh_ds['count_2020'] = annual_vh.count('time').sel(year=2020, drop=True)
        #
        dist_out_vh_ds['max_2017'] = annual_vh.max('time').sel(year=2017, drop=True)
        dist_out_vh_ds['max_2018'] = annual_vh.max('time').sel(year=2018, drop=True)
        dist_out_vh_ds['max_2019'] = annual_vh.max('time').sel(year=2019, drop=True)
        dist_out_vh_ds['max_2020'] = annual_vh.max('time').sel(year=2020, drop=True)
        #
        dist_out_vh_ds['min_2017'] = annual_vh.min('time').sel(year=2017, drop=True)
        dist_out_vh_ds['min_2018'] = annual_vh.min('time').sel(year=2018, drop=True)
        dist_out_vh_ds['min_2019'] = annual_vh.min('time').sel(year=2019, drop=True)
        dist_out_vh_ds['min_2020'] = annual_vh.min('time').sel(year=2020, drop=True)
        #
        dist_out_vh_ds['mean_2017'] = annual_vh.mean('time').sel(year=2017, drop=True)
        dist_out_vh_ds['mean_2018'] = annual_vh.mean('time').sel(year=2018, drop=True)
        dist_out_vh_ds['mean_2019'] = annual_vh.mean('time').sel(year=2019, drop=True)
        dist_out_vh_ds['mean_2020'] = annual_vh.mean('time').sel(year=2020, drop=True)
        #
        dist_out_vh_ds['median_2017'] = annual_vh.median('time').sel(year=2017, drop=True)
        dist_out_vh_ds['median_2018'] = annual_vh.median('time').sel(year=2018, drop=True)
        dist_out_vh_ds['median_2019'] = annual_vh.median('time').sel(year=2019, drop=True)
        dist_out_vh_ds['median_2020'] = annual_vh.median('time').sel(year=2020, drop=True)
        #
        dist_out_vh_ds['std_2017'] = annual_vh.std('time').sel(year=2017, drop=True)
        dist_out_vh_ds['std_2018'] = annual_vh.std('time').sel(year=2018, drop=True)
        dist_out_vh_ds['std_2019'] = annual_vh.std('time').sel(year=2019, drop=True)
        dist_out_vh_ds['std_2020'] = annual_vh.std('time').sel(year=2020, drop=True)
        #
        dist_out_vh_ds['q10_2017'] = annual_vh.quantile(0.1, 'time').sel(year=2017, drop=True).drop('quantile', dim=None)
        dist_out_vh_ds['q10_2018'] = annual_vh.quantile(0.1, 'time').sel(year=2018, drop=True).drop('quantile', dim=None)
        dist_out_vh_ds['q10_2019'] = annual_vh.quantile(0.1, 'time').sel(year=2019, drop=True).drop('quantile', dim=None)
        dist_out_vh_ds['q10_2020'] = annual_vh.quantile(0.1, 'time').sel(year=2020, drop=True).drop('quantile', dim=None)
        #
        dist_out_vh_ds['q90_2017'] = annual_vh.quantile(0.9, 'time').sel(year=2017, drop=True).drop('quantile', dim=None)
        dist_out_vh_ds['q90_2018'] = annual_vh.quantile(0.9, 'time').sel(year=2018, drop=True).drop('quantile', dim=None)
        dist_out_vh_ds['q90_2019'] = annual_vh.quantile(0.9, 'time').sel(year=2019, drop=True).drop('quantile', dim=None)
        dist_out_vh_ds['q90_2020'] = annual_vh.quantile(0.9, 'time').sel(year=2020, drop=True).drop('quantile', dim=None)
        #
        dist_out_vh_ds['q25_2017'] = annual_vh.quantile(0.25, 'time').sel(year=2017, drop=True).drop('quantile', dim=None)
        dist_out_vh_ds['q25_2018'] = annual_vh.quantile(0.25, 'time').sel(year=2018, drop=True).drop('quantile', dim=None)
        dist_out_vh_ds['q25_2019'] = annual_vh.quantile(0.25, 'time').sel(year=2019, drop=True).drop('quantile', dim=None)
        dist_out_vh_ds['q25_2020'] = annual_vh.quantile(0.25, 'time').sel(year=2020, drop=True).drop('quantile', dim=None)
        #
        dist_out_vh_ds['q75_2017'] = annual_vh.quantile(0.75, 'time').sel(year=2017, drop=True).drop('quantile', dim=None)
        dist_out_vh_ds['q75_2018'] = annual_vh.quantile(0.75, 'time').sel(year=2018, drop=True).drop('quantile', dim=None)
        dist_out_vh_ds['q75_2019'] = annual_vh.quantile(0.75, 'time').sel(year=2019, drop=True).drop('quantile', dim=None)
        dist_out_vh_ds['q75_2020'] = annual_vh.quantile(0.75, 'time').sel(year=2020, drop=True).drop('quantile', dim=None)
        #
        MAD_vh = annual_vh.map(lambda my_xr:  abs(my_xr - my_xr.median('time')).median('time'))
        #
        dist_out_vh_ds['MAD_2017'] = MAD_vh.sel(year=2017, drop=True)
        dist_out_vh_ds['MAD_2018'] = MAD_vh.sel(year=2018, drop=True)
        dist_out_vh_ds['MAD_2019'] = MAD_vh.sel(year=2019, drop=True)
        dist_out_vh_ds['MAD_2020'] = MAD_vh.sel(year=2020, drop=True)
        #------------------------------------------------------------------------
        # Save the calculated features in a NetCDF file  
        #------------------------------------------------------------------------
        outName_VV = iTileName + '_' + str(iRowStep) + '_' + str(iColStep) + '_' + str(iChunkSize) + '_' + orbitName + '_VV.nc'
        outName_VH = iTileName + '_' + str(iRowStep) + '_' + str(iColStep) + '_' + str(iChunkSize) + '_' + orbitName + '_VH.nc'
        #
        dist_out_vv_ds.to_netcdf(os.path.join(oDir, outName_VV))
        dist_out_vh_ds.to_netcdf(os.path.join(oDir, outName_VH))
        
    print()
    return print("The processing finished for the chunk row {myRow} and column {myCol} in tile {myTile}."
                 .format(myRow=str(iRowStep), myCol=str(iColStep), myTile=iTileName))

    
    
#######################################################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser( description=__doc__ , formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-p","--orbitpass", help="Orbit pass name (str), either: Ascending, Descending, or All [mandatory]", dest="iOrbitPass")
    parser.add_argument("-t","--tilename", help="Name (str) of the Equi7grid tile to be processed [mandatory]", dest="iTileName")
    parser.add_argument("-r","--rowstep", type=int, help="Row ID (int) of the tille cunck to be processed [mandatory]", dest="iRowStep")
    parser.add_argument("-c","--colstep", type=int, help="Column ID (int) of the tille cunck to be processed [mandatory]", dest="iColStep")
    parser.add_argument("-s","--chunksize", type=int, help="Cunck size (int) in pixels [mandatory]", dest="iChunkSize")
    parser.add_argument("-o","--odir", help="Path to the output directory [mandatory]", dest="oDir")
    args = parser.parse_args()
    #-------------
    process_Equi7tile(args.iOrbitPass, args.iTileName, args.iRowStep, args.iColStep, args.iChunkSize, args.oDir)