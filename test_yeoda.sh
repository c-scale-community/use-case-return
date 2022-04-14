# login to the spyder node 

# active the yeoda env.
conda activate yeoda

# go to python:
python
import os, osr, glob
from datetime import datetime
import matplotlib.pyplot as plt
#
from yeoda.products.preprocessed import SIG0DataCube
from geopathfinder.naming_conventions.yeoda_naming import YeodaFilename

# specify parameters:
tile_dirpath = r'/project/return/Share/EODC_SA020M/V01R01/E078N066T3'
#tile_dirpath = r'/project/return/Share/EODC_SA020M/V01R01/E051N060T3'
dimensions=['time', 'band', 'extra_field', 'sensor_field']

filepaths = glob.glob(os.path.join(tile_dirpath,'*.tif'))

# read the data cube:
sig0_dc = SIG0DataCube(filepaths=filepaths, dimensions=dimensions, filename_class=YeodaFilename, sres=20, continent='SA')

# get info:
sig0_dc.inventory[dimensions]

# specify the start and end dates:
toi_start, toi_end = datetime(2020, 1, 1), datetime(2021, 1, 1)

# filter the datacube (temporal):
sig0_dc.filter_by_dimension([(toi_start, toi_end)], [(">", "<")], name="time", inplace=True)

# select bands:
sig0_vv_dc = sig0_dc.filter_by_dimension('VV', name='band')
sig0_vh_dc = sig0_dc.filter_by_dimension('VH', name='band')

# get a time series for a point of interest:
#poi = (-8.345586, -78.70486)
poi = (-2.70919, -54.95610)
sref = osr.SpatialReference()
sref.ImportFromEPSG(4326)
sig0_vv_ts = sig0_vv_dc.load_by_coords(*poi, sref=sref, dtype='numpy')
sig0_vh_ts = sig0_vh_dc.load_by_coords(*poi, sref=sref, dtype='numpy')

# plot a timeseries:
plt.scatter(sig0_vh_dc['time'].to_numpy(), sig0_vh_ts[:, 0, 0]/100)
plt.show()


# get a list of orbits in tile:
sig0_dc.inventory[dimensions].extra_field.unique() 

# select a single orbit: 
sig0_vv_desc_dc = sig0_vv_dc.filter_by_dimension('D039', name='extra_field')
sig0_vv_s1b_ts = sig0_vv_desc_dc.load_by_coords(*poi, sref=sref, dtype='numpy')

# plot a timeseries:
plt.scatter(sig0_vv_desc_dc['time'].to_numpy(), sig0_vv_s1b_ts[:, 0, 0]/100)
plt.show()


