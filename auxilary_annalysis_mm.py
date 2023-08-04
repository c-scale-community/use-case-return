# Author:       Milutin Milenkovic
# Created:      10/02/2023
# Copyright:    ???
# Licence:      ???


def multivarFit_s1_forAge(my_xr_ds):
    # modules
    import xarray as xr
    from sklearn.linear_model import LinearRegression
    # select secondary forest and convert to pandas data frame
    my_df = my_xr_ds.where((my_xr_ds['forAge2017'] > 0) & (my_xr_ds['forAge2017'] < 50)).to_dataframe()
    # drop unneccery data:
    my_df.drop(columns='spatial_ref', inplace=True)
    my_df.dropna(subset=['forAge2017'], inplace=True)
    # group by year:
    df_med = my_df.groupby(by='forAge2017').mean().reset_index()
    # prepare data for the MULTIVARIATE fit
    X = df_med.drop(['forAge2017'], axis = 1)
    y = df_med.reset_index()[['forAge2017']]
    # fit the model
    lm = LinearRegression()
    model = lm.fit(X,y)
    return model.score(X,y)



def plot_featSpaceDist(c_min, c_max, r_min, r_max):
    import os
    import xarray as xr
    import rioxarray
    from rioxarray.merge import merge_datasets, merge_arrays
    import numpy as np
    import geopandas as gpd
    import pandas as pd
    #
    from sklearn.preprocessing import MinMaxScaler
    #
    import matplotlib.pyplot as plt
    # ----------------------------------------
    myTile = 'E078N066T3' # Para
    myOrbit = 'D039' # Para
    myPol = 'VH'
    # ----------------------------------------
    # set the folter wih the output files:
    data_folder = r'/project/return/Share/mm/S1_SA_TEST_UPSCALE/TILE_WISE/AOI_PA/' +  myTile + '_smallChanks'
    # Equi7 SA wkt:
    PROJ = 'PROJCS["Azimuthal_Equidistant",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.017453292519943295]],PROJECTION["Azimuthal_Equidistant"],PARAMETER["false_easting",7257179.23559],PARAMETER["false_northing",5592024.44605],PARAMETER["central_meridian",-60.5],PARAMETER["latitude_of_origin",-14.0],UNIT["Meter",1.0]]'
    # ----------------------------------------
    # read multiple ckunks:
    ChunkDatasetList = []
    # set the rwos and columns:
    #myCols, myRows = np.meshgrid(np.arange(34, 38), np.arange(42, 46)) # Rondonia
    myCols, myRows = np.meshgrid(np.arange(c_min, c_max), np.arange(r_min, r_max)) # Para
    # read each one
    for myRow, myCol in zip(myRows.flatten(), myCols.flatten()):
        # set filename:
        fileName = myTile + '_' + str(myRow) + '_' + str(myCol) + '_250_' + myOrbit + '_'+ myPol + '.nc'
        if os.path.exists(os.path.join(data_folder, fileName)):
            aux_ds = rioxarray.open_rasterio(os.path.join(data_folder, fileName))
            aux_ds.rio.write_crs(PROJ, inplace=True)
            ChunkDatasetList.append(aux_ds)
        # set crs:
           
    # ----------------------------------------
    # mearge chanks spatialy:
    my_out = merge_datasets(ChunkDatasetList)
    # ----------------------------------------
    # Read forest age and LULC data:
    # fileName2 = r'myforAgeLCLU_2017_2020_equi7_E066N060T3_20m.tif' # Rondonia
    fileName2 = r'myforAgeLCLU_2017_2020_equi7_E078N066T3_20m.tif' # Para
    forAgeLULC = rioxarray.open_rasterio(fileName2, band_as_variable=True)
    # convert to dataset:
    forAgeLULC_ds = forAgeLULC.to_dataset('band').rename({1: 'forAge2017', 2: 'lulc2018', 3: 'lulc2019', 4: 'forAge2020', 5: 'ECO_ID'})
    # ----------------------------------------
    # Correct coordinates (temporary soluton!!!)
    forAgeLULC_ds = forAgeLULC_ds.assign_coords(x=forAgeLULC_ds.x - 10 )
    forAgeLULC_ds = forAgeLULC_ds.assign_coords(y=forAgeLULC_ds.y + 10 )
    # ----------------------------------------
    # Join datasets and clip forAgeLCLU to my_out coordinates:
    all_ds = xr.merge([my_out.isel(band=0, drop=True), forAgeLULC_ds], join='inner') 
    # ############################# 
    # Add annual statistics
    # ############################# 
    # Add the combined annual statistics:
    all_ds['qDiff_2017'] = all_ds['q90_2017'] - all_ds['q10_2017']
    all_ds['dRange_2017'] = all_ds['max_2017'] - all_ds['min_2017']
    all_ds['iqRange_2017'] = all_ds['q75_2017'] - all_ds['q25_2017']
    # set the relevant statistics:
    myStatnames = ['mean_2017', 'median_2017', 'std_2017', 'MAD_2017', 'dRange_2017', 'qDiff_2017', 'iqRange_2017', 'forAge2017']
    # mearge the arrays into new dataset:
    my_ds = all_ds[myStatnames]
    # select secondary and OLD forest and convert to pandas data frame
    my_df = my_ds.where((my_ds['forAge2017'] > 0) | (my_ds['forAge2017'] == 50)).to_dataframe()
    # drop unneccery data:
    my_df.drop(columns='spatial_ref', inplace=True)
    my_df.dropna(subset=['forAge2017'], inplace=True)
    # ############################# 
    # Split Old forest in negative age classes:
    # #############################
    # get secondary forest pixels
    my_df_secFor = my_ds.where((my_ds['forAge2017'] > 0) & (my_ds['forAge2017'] < 50)).to_dataframe()
    my_df_secFor.drop(columns='spatial_ref', inplace=True)
    my_df_secFor.dropna(subset=['forAge2017'], inplace=True)
    # get old forest pixels:
    my_df_oldFor = my_ds.where(my_ds['forAge2017'] == 50).to_dataframe()
    my_df_oldFor.drop(columns='spatial_ref', inplace=True)
    my_df_oldFor.dropna(subset=['forAge2017'], inplace=True)
    # #############################
    # Feature space conversion:
    # #############################
    my_refData = my_df_oldFor.reset_index()
    my_refData.drop(columns=['y', 'x', 'forAge2017'], inplace=True)
    # get the mean of the refernce data
    myRef_mean = my_refData.mean()
    myRef_mean['forAge2017'] = 0
    # get the mean of the secFor data:
    my_df_secFor = my_df_secFor.reset_index()
    secFor_mean = my_df_secFor.groupby(by='forAge2017').mean()
    my_secForData = secFor_mean.drop(columns=['y', 'x'])
    # append secFor data with ref data
    my_Data = my_secForData.reset_index().append(myRef_mean, ignore_index=True)
    # define ref vector:
    ref = my_Data[my_Data.forAge2017==0].drop(columns=['forAge2017']).values
    # calcilate distance per each row, i.e., forest age
    euqlidianDist = my_Data.drop(columns=['forAge2017']).apply(lambda row: np.linalg.norm(row-ref), axis=1, raw=True)
    # assign distance to df
    my_Data['Dist'] = euqlidianDist
    # #############################
    # Scale the data
    # #############################
    # make a copy of the data witout the exsisting distance:
    my_Data_v2 = my_Data.drop(columns=['Dist', 'forAge2017']).copy()
    # initiate a model to normalize the data: 
    my_norm = MinMaxScaler().fit(my_Data_v2)
    # get the normalized data as nd array:
    aux_norm = my_norm.transform(my_Data_v2)
    # make a df of the normalized data:
    my_Data_scaled = pd.DataFrame(aux_norm, columns=my_Data_v2.columns)
    # calculate the distance from the normalized data:
    norDist = my_Data_scaled.apply(lambda row: np.linalg.norm(row-ref), axis=1, raw=True)
    # assign back the droped colimns and normalized distance:
    my_Data_scaled['forAge2017'] = my_Data['forAge2017'].values
    my_Data_scaled['Dist'] = my_Data['Dist'].values
    my_Data_scaled['norDist'] = norDist
    # #############################
    # plotting
    # #############################
    fig1, ax1 = plt.subplots(1, 1)
    my_Data_scaled[my_Data_scaled.forAge2017>0].plot(x='forAge2017', y='Dist', style='bo-', ax=ax1)
    # -----
    fig2, ax2 = plt.subplots(1, 1)
    my_Data_scaled[my_Data_scaled.forAge2017>0].plot(x='forAge2017', y='norDist', style='bo-', ax=ax2)
    return