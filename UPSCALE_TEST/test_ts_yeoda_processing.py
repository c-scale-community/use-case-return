#
# based on:
# https://github.com/c-scale-community/use-case-return/blob/master/test_ts_yeoda_processing.ipynb
#

import os, osr, glob
import numpy as np
# ------
# Choose a backend for matplotlib that doesn't display to the user, e.g. 'Agg' (note plots are overwritten automatically upon re-running)
import matplotlib
matplotlib.use('Agg')
# ------
import matplotlib.pyplot as plt
#%matplotlib inline # jupyter notebook specific (i.e. magic function %matplotlib inline to enable the inline plotting, where the plots/graphs will be displayed just below the cell where your plotting commands are written).
from datetime import datetime
import pandas as pd

# import TUW packages
from yeoda.products.preprocessed import SIG0DataCube
from geopathfinder.naming_conventions.yeoda_naming import YeodaFilename
#
#%load_ext autoreload # jupyter notebook specific
#%autoreload 2 # jupyter notebook specific
#%reload_ext autoreload # jupyter notebook specific

# # MM AUX functions
# import my aux functions
from auxilary_ts_tools_mm import features_from_yeodaTS, plot_TS, features_from_S1_TS

# Specify the folder with a S1 datacube (a 300x300 m2 Equi7Tile)
tile_dir1_path = r'/project/return/Share/EODC_SA020M/V01R01/E078N066T3'
tile_dir2_path = r'/project/return/Share/EODC_SA020M/V1M0R1/E078N066T3'
#
# specify other parameters:
dimensions=['time', 'band', 'extra_field', 'sensor_field']
filepaths1 = glob.glob(os.path.join(tile_dir1_path,'*.tif'))
filepaths2 = glob.glob(os.path.join(tile_dir2_path,'*.tif'))

# Read the datacube:
sig0_dc1 = SIG0DataCube(filepaths=filepaths1, dimensions=dimensions, filename_class=YeodaFilename, sres=20, continent='SA')
sig0_dc2 = SIG0DataCube(filepaths=filepaths2, dimensions=dimensions, filename_class=YeodaFilename, sres=20, continent='SA')
#
# get info:
sig0_dc1.inventory[dimensions].head(5)

# Filter by date:
toi_start, toi_end = datetime(2017, 1, 1), datetime(2021, 1, 1)
sig0_dc1 = sig0_dc1.filter_by_dimension([(toi_start, toi_end)], [(">=", "<")], name="time", inplace=True)
sig0_dc2 = sig0_dc2.filter_by_dimension([(toi_start, toi_end)], [(">=", "<")], name="time", inplace=True)

# Select bands:
sig0_vv_dc1 = sig0_dc1.filter_by_dimension('VV', name='band')
sig0_vh_dc1 = sig0_dc1.filter_by_dimension('VH', name='band')
#
sig0_vv_dc2 = sig0_dc2.filter_by_dimension('VV', name='band')
sig0_vh_dc2 = sig0_dc2.filter_by_dimension('VH', name='band')

# Merge and sort the datacubes:
sig0_vv_dc = sig0_vv_dc1.unite(sig0_vv_dc2)
sig0_vv_dc = sig0_vv_dc.sort_by_dimension('time', ascending=True)
#
sig0_vh_dc = sig0_vh_dc1.unite(sig0_vh_dc2)
sig0_vh_dc = sig0_vh_dc.sort_by_dimension('time', ascending=True)
#
sig0_vv_dc.inventory

# Get a time series for a point of interest:
poi = (-3.48472, -54.82820)
sref = osr.SpatialReference()
sref.ImportFromEPSG(4326)
#
sig0_vv_ts = sig0_vv_dc.load_by_coords(*poi, sref=sref, dtype='dataframe')
sig0_vh_ts = sig0_vh_dc.load_by_coords(*poi, sref=sref, dtype='dataframe')

# Remane the variable and drop na as well as indexing:
sig0_vv_ts = sig0_vv_ts.dropna().droplevel(['x', 'y']).rename(columns={'1' : 'vv'})
sig0_vh_ts = sig0_vh_ts.dropna().droplevel(['x', 'y']).rename(columns={'1' : 'vh'})
# plot:
my_xticks = pd.date_range(datetime(2017,1,1), datetime(2021,1,1), freq='YS')
ax2 = sig0_vv_ts.plot(style='r-', grid=True)
sig0_vh_ts.plot(style='b-', grid=True, ax=ax2)
plt.savefig('/home/return-roonk/S1MM/test/plots/fig1.png')
plt.close()

# Rescale the data in 2019 and 2020
sig0_vv_ts.loc['2019-1-1':'2021-1-1'] = sig0_vv_ts.loc['2019-1-1':'2021-1-1'].div(10.).round(1)
sig0_vh_ts.loc['2019-1-1':'2021-1-1'] = sig0_vh_ts.loc['2019-1-1':'2021-1-1'].div(10.).round(1)
# plot to check:
my_xticks = pd.date_range(datetime(2017,1,1), datetime(2021,1,1), freq='YS')
ax1 = sig0_vv_ts.plot(style='ro-', xticks=my_xticks, grid=True, figsize=(14,4), legend=True, xlabel='Time', ylabel='Bacscatter Intensity [dB]')
sig0_vh_ts.plot(style='bo-', ax=ax1, xticks=my_xticks, grid=True, figsize=(14,4), legend=True, xlabel='Time', ylabel='Bacscatter Intensity [dB]')
plt.savefig('/home/return-roonk/S1MM/test/plots/fig2.png')
plt.close()

# See the difference in time to decde on the minum sample interval:
np.unique(sig0_vh_ts.index.round('D').to_series().diff().values.astype('timedelta64[D]'))

# Round the time and resample to 6 day TS:
sig0_vh_ts.index = sig0_vh_ts.index.round('D')
#
sig0_vh_ts_6d = sig0_vh_ts.resample('6D').interpolate(method='linear')
sig0_vh_ts_6d.plot(style='bo-',grid=True, figsize=(14,4), legend=True, xlabel='Time', ylabel='Bacscatter Intensity [dB]')
plt.savefig('/home/return-roonk/S1MM/test/plots/fig3.png')
plt.close()

# # MM AUX functions
# # Plot and get the time seres features:
# plot_TS(sig0_vh_ts_6d)
# #
# myFeatures = features_from_yeodaTS(sig0_vh_ts_6d)
# #
# myFeatures

# # Test simple plot
# print("create simple plot")
# # Data for plotting
# t = np.arange(0.0, 2.0, 0.01)
# s = 1 + np.sin(2 * np.pi * t)
# fig, ax = plt.subplots()
# ax.plot(t, s)
# ax.set(xlabel='x-value', ylabel='y-value', title='Simple plot')
# ax.grid()
# fig.savefig("simple.png")
# plt.show()

print("done")
