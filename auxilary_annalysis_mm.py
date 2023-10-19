# Author:       Milutin Milenkovic
# Created:      10/02/2023
# Copyright:    ???
# Licence:      ???



def getSlopeInt_qDiff_forAge(myTile, myOrbit, myPol, NumChunks, ChunkStart_col, ChunkStart_row, ):
    """
    
    Calculate the slope and intercepth of the linear regression model fitted trough 
    the annual percentile diference (p90 - p10) mean-aggregated over the secoundary forest age.
    
    :param myTile: The Equi7 tile name to be processed, e.g. 'E078N066T3'
    :param myOrbit: The decending orbit to be processed, e.g. 'D039'
    :param myPol: the polarisation chanel (as a string), e.g. 'VV', or 'VH'
    :param NumChunks: the number of chunks to consider for the calculation of slope and intercept 
    :param ChunkStart_col: the column of the starting chank
    :param ChunkStart_row: the row of the starting chank
    :return: a list including linear regresion parameters, performance stats, and nuber of secoundary pixels
    
    """
    # modules
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
    from sklearn.linear_model import LinearRegression
    from sklearn.metrics import mean_squared_error
    # ------------------
    # set parameters:
    # ------------------
    #
    # specify variables to read from the NetCDF files:
    myVariable_list = ['q90_2017', 'q10_2017']
    # specify the number of chunks to consider for the calculation of slope and intercept:
    NumChunks = 4
    # specify the starting chank:
    ChunkStart = 1
    # -------
    # set  path to the netCDF files:
    data_folder = r'/project/return/Share/mm/S1_SA_TEST_UPSCALE/TILE_WISE/AOI_PA/' +  myTile
    # set path to the Forest Age and LULC file:
    fileName2 = r'/home/return-mmilenkovic/mm_tools/use-case-return/input/myforAgeLCLU_2017_2020_equi7_E078N066T3_20m.tif'
    # set the Equi7 projection WKT string:
    PROJ='PROJCS["Azimuthal_Equidistant",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.017453292519943295]],PROJECTION["Azimuthal_Equidistant"],PARAMETER["false_easting",7257179.23559],PARAMETER["false_northing",5592024.44605],PARAMETER["central_meridian",-60.5],PARAMETER["latitude_of_origin",-14.0],UNIT["Meter",1.0]]'
    # -----------------------------------
    # read and mearge ckunks (nc-files):
    # -----------------------------------
    ChunkDatasetList = []
    # set the rwos and columns:
    myCols, myRows = np.meshgrid(np.arange(ChunkStart_col, ChunkStart_col + NumChunks), 
                                 np.arange(ChunkStart_row, ChunkStart_row + NumChunks))
    # read each one into a list
    for myRow, myCol in zip(myRows.flatten(), myCols.flatten()):
        # set filename:
        fileName = myTile + '_' + str(myRow) + '_' + str(myCol) + '_250_' + myOrbit + '_'+ myPol + '.nc'
        folder_name = str(myRow + 1) + '_' + str(myCol + 1)
        if os.path.exists(os.path.join(data_folder, folder_name, fileName)):
            aux_ds = rioxarray.open_rasterio(os.path.join(data_folder, folder_name, fileName), variable=myVariable_list)
            # set crs:
            aux_ds.rio.write_crs(PROJ, inplace=True)
            ChunkDatasetList.append(aux_ds)
    
    # mearge chanks spatialy:
    my_out = merge_datasets(ChunkDatasetList)   
    # -----------------------------------
    # Read Forest Age and LULC data:
    # -----------------------------------
    forAgeLULC = rioxarray.open_rasterio(fileName2, band_as_variable=True)
    # convert to dataset:
    forAgeLULC_ds = forAgeLULC.to_dataset('band').rename({1: 'forAge2017', 2: 'lulc2018', 3: 'lulc2019', 4: 'forAge2020', 5: 'ECO_ID'})
    # -----------------------------------
    # Correct corrdinates
    # NOTE: !!! REVISIT THIS AGAIN !!!
    # -----------------------------------
    forAgeLULC_ds = forAgeLULC_ds.assign_coords(x=forAgeLULC_ds.x - 10 )
    forAgeLULC_ds = forAgeLULC_ds.assign_coords(y=forAgeLULC_ds.y + 10 )
    # ---------------------------------------------------
    # Join two datasets and calculate the qDiff
    # ----------------------------------------------------
    all_ds = xr.merge([my_out.isel(band=0, drop=True), forAgeLULC_ds], join='inner')
    # annual statistics:
    all_ds['qDiff_2017'] = all_ds['q90_2017'] - all_ds['q10_2017']
    # ----------------------------------------------------
    # Convert only Forest pixels to the panda dataframe
    # ----------------------------------------------------
    # select secondary and OLD forest and convert to pandas data frame
    my_df = all_ds.where((all_ds['forAge2017'] > 0) | (all_ds['forAge2017'] == 50)).to_dataframe()
    # drop unneccery data:
    my_df.drop(columns='spatial_ref', inplace=True)
    my_df.dropna(subset=['forAge2017'], inplace=True)
    # ----------------------------------------------------
    # Split Old forest and Secoundary forest
    # ----------------------------------------------------
    # get secondary forest pixels
    my_df_secFor = all_ds.where((all_ds['forAge2017'] > 0) & (all_ds['forAge2017'] < 50)).to_dataframe()
    my_df_secFor.drop(columns='spatial_ref', inplace=True)
    my_df_secFor.dropna(subset=['forAge2017'], inplace=True)
    # get old forest pixels:
    my_df_oldFor = all_ds.where(all_ds['forAge2017'] == 50).to_dataframe()
    my_df_oldFor.drop(columns='spatial_ref', inplace=True)
    my_df_oldFor.dropna(subset=['forAge2017'], inplace=True)
    # ----------------------------------------------------
    # Get the number of avaliable secoundary forest pixels
    # ----------------------------------------------------
    # get the number of secFor pixels per age
    secFor_count = my_df_secFor.groupby(by='forAge2017').count()
    # get the total numbert of secoundary forest pixels:
    NumTotalSecFor = secFor_count.q90_2017.values.sum()
    # get the average num of secFor pixels in 5 count-largest ages:
    myPixNu_largest = np.int32(np.floor(np.mean(secFor_count.q90_2017.nlargest(5).values)))
    # get the average num of secFor pixels in 5 count-smallest ages:
    myPixNum_smallest = np.int32(np.floor(np.mean(secFor_count.q90_2017.nsmallest(5).values)))
    # get the average num of secFor pixels in the first 5 years:
    #myPixNu_first5 = np.int32(np.floor(np.mean(secFor_count[:5].q90_2017.values)))
    # get the average num of secFor pixels in the first 5 years:
    #myPixNu_last5 =np.int32(np.floor(np.mean(secFor_count.tail(5).q90_2017.values)))
    #
    # --------------------------------------------------------
    # Aggregate the S1 qDiff per Forest Age class using mean:
    # --------------------------------------------------------
    # get the mean of the secFor data:
    my_df_secFor = my_df_secFor.reset_index()
    secFor_mean = my_df_secFor.groupby(by='forAge2017').mean()
    my_secForData = secFor_mean.drop(columns=['y', 'x'])
    # ----------------------------------------------------
    # Fit the linear regression model:
    # ----------------------------------------------------
    # prepare data for the MULTIVARIATE fit
    X = my_secForData.reset_index()[['forAge2017']]
    y = my_secForData.reset_index()[['qDiff_2017']]
    # fit the model
    lm = LinearRegression()
    model = lm.fit(X,y)
    # ----------------------------------------------------
    # Get the model parameters and performance stats
    # ----------------------------------------------------
    mySlope = model.coef_[0]
    myInt = model.intercept_[0]
    # the coef. of determination (R2)
    myR2 = model.score(X,y)
    # the RMSE
    myRMSE = mean_squared_error(y, model.predict(X), squared=False)
    #
    return [mySlope[0], myInt, myR2, myRMSE, NumTotalSecFor, myPixNu_largest, myPixNum_smallest]


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