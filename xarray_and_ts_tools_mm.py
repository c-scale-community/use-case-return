# Author:       Milutin Milenkovic
# Created:      23/06/2020
# Copyright:    ???
# Licence:      ???

def label_lossyear(ds):
    """ Label polygons in Hansen GEE dataset for the lossyear-s 2017, 2018, and 2019.
    The input is an xarray Data Set derived using either:
     (a) collectionMultiband2xarray_ds.py, or
     (b) collection2xarray_ds.py """
    # --------------------
    from skimage import measure
    # --------------------
    if 'LossYear' not in ds.coords:
        raise Exception('The coordinate "LossYear" not defined in the input xarray Data Set')

    # --------------------
    LossYear = 17
    # get the connected components (regions) and label them:
    mask_2017_labels = measure.label(1*(ds.LossYear.values == LossYear), background=0)
    # get the region properties (a python list with objects)
    regions_2017 = measure.regionprops(mask_2017_labels)
    regions_2017_sorted = sorted(regions_2017, key=lambda x: x.area, reverse=True)
    # --------------------
    LossYear = 18
    # get the connected components (regions) and label them:
    mask_2018_labels = measure.label(1*(ds.LossYear.values == LossYear), background=0)
    # get the region properties (a python list with objects)
    regions_2018 = measure.regionprops(mask_2018_labels)
    regions_2018_sorted = sorted(regions_2018, key=lambda x: x.area, reverse=True)
    # --------------------
    LossYear = 16
    # get the connected components (regions) and label them:
    mask_2016_labels = measure.label(1*(ds.LossYear.values == LossYear), background=0)
    # get the region properties (a python list with objects)
    regions_2016 = measure.regionprops(mask_2016_labels)
    regions_2016_sorted = sorted(regions_2016, key=lambda x: x.area, reverse=True)
    # -----------------------------------------------------------------------
    # add the region labels as new coordinate in dataset:
    # -----------------------------------------------------------------------
    ds.coords['RegionLabels_2017'] = (('y', 'x'), mask_2017_labels)
    ds.coords['RegionLabels_2018'] = (('y', 'x'), mask_2018_labels)
    ds.coords['RegionLabels_2016'] = (('y', 'x'), mask_2016_labels)
    # -----------------------------------------------------------------------
    return ds


def disturbance_parameters(ts_pix_1,ts_ref,sample_num,alpha1):
    """
    This function takes two time seres (TSs), a reference TS and a disturbed-pixel TS,
    and derives different disturbance parameters.
    ----------
    The flowing steps are performed:
    (a) calculate the rolling mean from input TSs
            - parameter #1: 'sample_num' defines the sample length of the rolling window
    (b) calculate the difference of the two mean TSs (ref. TS - single-pix TS)
    (c) calculate the statistics of the TS residuals and derive T-test statistics
            - parameter #2: 'alpha1' confidence probability
    (d) threshold the difference according to t-test margin and extract signal segments above the threshold
    (e) derive the disturbance parameters
    ----------
    The output is a list of the following values:
    - max_mag:  the maximum magnitude (dB, float) of the largest disturbance segment derived from the running mean TS
    - t_pre:    time difference (days, Int) between the maximum magnate and the start of the largest segment
    - t_post:   time difference (days, Int) between the maximum magnate and the end of the largest segment
    - max_mag_date:     absolute date (Y-M-D, numpy.datetime64) of 'max_mag'
    - max_mag_org:      maximum magnitude (dB, float) of the largest disturbance segment derived from the original TS
    - t_mag_org:        time difference (days, Int) between 'max_mag' and 'max_mag_org'
    - max_mag_org_date: absolute date (Y-M-D, numpy.datetime64) of 'max_mag_org'
    - seg2_size:        the size of the 2nd segment
    """
    # -------------------------------------------
    from scipy import stats
    import numpy as np
    from skimage import measure
    # -------------------------------------------
    # check if the single-pixel TS has valid values at all:
    if np.all(np.isnan(ts_pix_1.values)):
        return [1000, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
    # ------------------ (a) --------------------
    # get the rolling means
    ts_ref_mean = ts_ref.rolling(time=sample_num, center=True).mean()
    ts_pix_1_mean = ts_pix_1.rolling(time=sample_num, center=True).mean()
    # ------------------ (b) --------------------
    # get the difference between reference and single-pixel mean TS-s:
    ts_diff = ts_ref_mean - ts_pix_1_mean
    # ------------------ (c) --------------------
    # get the residual TSs:
    ts_ref_res = ts_ref - ts_ref_mean
    ts_pix_1_res = ts_pix_1 - ts_pix_1_mean
    # get the std of the residual TSs:
    ref_res_std = np.nanstd(ts_ref_res.values)
    pix_1_res_std = np.nanstd(ts_pix_1_res.values)
    # get the T-test statistics (critical value and error margin):
    deg_of_freedom = sample_num + sample_num - 2
    t_critical = stats.t.ppf(q=1 - alpha1 / 2, df=deg_of_freedom)
    error_margin = t_critical * np.sqrt((ref_res_std ** 2.) / sample_num + (pix_1_res_std ** 2.) / sample_num)
    # ------------------ (d) --------------------
    # dealing with nan-s values due to the running mean
    ts_diff_values = ts_diff.values
    ts_diff_values[np.isnan(ts_diff_values)] = 0
    diff_bool = ts_diff_values > error_margin
    # -----------------------------------------------------------
    # different checks for extra cases:
    # -------
    # check if the first value is already above margin (indicates a disturbance outside the time span):
    ind_1st_valid = int(sample_num/2) # note!!!! this works only for the even sample_num
    if abs(ts_diff_values[ind_1st_valid]) > error_margin:
        return [1002, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, error_margin]
    # -------
    # check if the first value is above margin and negative (indicates a TS above the ref_TS):
    if abs(ts_diff_values[ind_1st_valid]) > error_margin and ts_diff_values[ind_1st_valid] < 0:
        return [1003, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, error_margin]
    # -------
    # check if there is at least one element above the margin:
    if np.all(diff_bool == False):
        return [1001, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, error_margin]
    # -----------------------------------------------------------
    # identify time segments (consecutive magnitudes) in the thresholded difference :
    diff_bool_labels = measure.label(diff_bool)
    # as measure.regionprops only works for images, an extra auxiliary row is added here:
    diff_bool_labels_2D = np.array([diff_bool_labels, diff_bool_labels])
    # get the segments' properties
    diff_bool_regions = measure.regionprops(diff_bool_labels_2D)
    # get size of each segment as array
    seg_size = np.array([aa.area for aa in diff_bool_regions])
    # get the indices of the start and the end of the largest segment:
    seg_ind_start = diff_bool_regions[seg_size.argmax()].bbox[1]
    seg_ind_end = diff_bool_regions[seg_size.argmax()].bbox[3]
    # ------------------ (e) --------------------
    # get the maximum disturbance magnitude in the largest segment:
    ts_diff_seg1 = ts_diff.where(diff_bool_labels == diff_bool_regions[seg_size.argmax()].label)
    max_mag = np.nanmax(ts_diff_seg1.values)
    max_mag_ind = np.nanargmax(ts_diff_seg1.values)
    max_mag_date = ts_diff_seg1.coords['time'].values[max_mag_ind]
    # get the relative time from the maximum disturbance:
    d_time = (ts_diff.coords['time'].values - max_mag_date).astype('timedelta64[D]').astype(int)
    # get the disturbance parameters related to the start and the end of the largest segment:
    t_pre = d_time[seg_ind_start]
    t_post = d_time[seg_ind_end]
    # get the total number of segments:
    num_of_segments = len(seg_size)
    # size of the 2nd largest segment:
    if len(seg_size) > 1:
        seg2_size = seg_size[1]/2.
        # the division by 2 is due to the extra array row added above to trick 'measure.regionprops'
    else:
        seg2_size = 0
    # get parameters from the original TS (no running mean applied):
    diff_org = ts_ref_mean - ts_pix_1
    diff_org_seg1 = diff_org.where(diff_bool_labels == diff_bool_regions[seg_size.argmax()].label)
    # get the maximum magnitude
    max_mag_org = np.nanmax(diff_org_seg1.values)
    max_mag_org_ind = np.nanargmax(diff_org_seg1.values)
    max_mag_org_date = diff_org_seg1.coords['time'].values[max_mag_org_ind]
    t_mag_org = d_time[max_mag_org_ind]
    # -------------------------------------------
    return [max_mag, t_pre, t_post, max_mag_org, t_mag_org, seg2_size, max_mag_date, max_mag_org_date, error_margin]


def get_xy_from_onclick(event):
    """
    This function is used to digitize manually a certain number of pixels from a python plot.
    The number of pixels to digitize is hardcoded below with the parameter 'click_num'.
    The python plot should be active and the following code should be run:
    my_coords = []
    cid = fig.canvas.mpl_connect('button_press_event', get_xy_from_onclick)
    """
    x = event.x
    xdata = event.xdata
    y = event.y
    ydata = event.ydata
    #
    print('a pixel x = %d, y = %d  was selected' % (xdata,ydata))
    #
    global my_coords
    my_coords.append([xdata, ydata])
    # Disconnect after specified number of clicks
    click_num = 2
    if len(my_coords) == click_num:
        fig.canvas.mpl_disconnect(cid)
    #
    return my_coords


def plot_pixel_ts_loc(my_coor, ds, LOSS_YEAR, POLARIZATION):
    """
    The function plots a single-pixel TS, reference TS and the spatial location of the pixel. This plot is used
    as an easy way to visually inspect the shape of the TS.

    :param my_coor: a list with x and y coordinates of the pixel of interest (PoI)
    :param ds: an xarray Dataset derived from the Sentinel-1 collection
    :param LOSS_YEAR: the losyear value from the Hansen data
    :param POLARIZATION: the polarisation chanel (as a string), e.g. 'VV', or 'VH'
    :return: a Python plot with the PoI TS, reference TS, and spatial location of the pixel
    """
    import numpy as np
    from matplotlib import pyplot as plt
    from matplotlib import cm
    import pandas as pd
    # define some plotting parameters
    start_date = np.datetime64('2016-01-01')
    end_date = np.datetime64('2020-01-01')
    datelist = pd.date_range(start=start_date, end=end_date, freq='YS').tolist()
    #
    gfc_cmap = cm.get_cmap('viridis', 20).reversed()
    gfc_cmap.colors[0][-1] = 0
    # -------------
    # select the ts:
    ts_pix_1 = ds[POLARIZATION].sel(dict(x=my_coor[0], y=my_coor[1]), method='nearest')
    ts_ref = ds[POLARIZATION].where(ds.LossYear.values == 0).mean(dim=['x','y'])
    # -------------
    fig, (ax_ts, ax_map) = plt.subplots(2,figsize=(8, 9), gridspec_kw={"height_ratios": (.2, .8)})
    ts_pix_1.plot(color='grey',ax=ax_ts)
    ts_ref.plot(color='green',ax=ax_ts)
    ax_ts.set_xticks(datelist)
    ax_ts.set_axisbelow(True)
    ax_ts.set_yticks(np.arange(-20, -9, 1))
    ax_ts.minorticks_on()
    ax_ts.grid(which='major', linestyle='-', linewidth=1)
    ax_ts.grid(which='minor', linestyle=':', linewidth=0.5)
    ax_ts.legend(['disturbed pixels TS','non-disturbed pixels TS'],loc="lower left")
    ax_ts.set_title(r'Location: x={}, y={} '.format(my_coor[0],my_coor[1]))
    plt.tight_layout()
    #
    #ds.LossYear.plot(ax=ax_map, cmap=gfc_cmap, levels=np.arange(20), cbar_kwargs={"ticks": np.arange(20)}, add_colorbar=False)
    ds.LossYear.where(ds.LossYear == LOSS_YEAR).plot(ax=ax_map, cmap=gfc_cmap, levels=np.arange(20), add_colorbar=False)
    plt.scatter(my_coor[0],my_coor[1], marker='o', s=20, color='red', axes=ax_map)
    ax_map.set_aspect('equal', 'box')
    ax_map.grid(linestyle=':', linewidth=0.5)
    aux1_x = my_coor[0]
    ax_map.set_xlim(aux1_x - 1500, aux1_x + 1500)
    aux1_y = my_coor[1]
    ax_map.set_ylim(aux1_y - 1500, aux1_y + 1500)
    plt.tight_layout()


def plot2_pixel_ts_loc(my_coor, ts_pix_1, ts_ref, ds, alpha1, sample_num, error_margin1, dist_pix_ind, LOSS_YEAR):
    """
    The function plots a single-pixel TS, reference TS and the spatial location of the pixel. This plot is used
    as an easy way to visually inspect the shape of the TS.

    The difference from 'plot_pixel_ts_loc' is that here the mean TSs and its difference TS are also plotted.

    :param my_coor: a list with x and y coordinates of the pixel of interest (PoI)
    :param ts_pix_1: a time series of the pixel of interest (PoI) given as xarray Dataset
    :param ts_ref: the refernce time series derived as xarray Dataset
    :param ds: an xarray Dataset derived from the Sentinel-1 collection
    :param alpha1: the confidence probability for the two-sided t-test of the difference between two means
    :param sample_num: the number of samples used to derive the running mean time series
    :param error_margin1: the error margin associated to the t-test critical value
    :param dist_pix_ind: the index of a randomly-selected pixel from the list of the all disturbed pixels
    :param LOSS_YEAR: the losyear value from the Hansen data
    :return: a Python plot with the PoI TS, reference TS, and spatial location of the pixel
    """
    import numpy as np
    from matplotlib import pyplot as plt
    from matplotlib import cm
    import pandas as pd
    from scipy import stats
    # get the running means, difference signals, and confidence intervals for plotting
    ts_ref_mean = ts_ref.rolling(time=sample_num, center=True).mean()
    ts_ref_res = ts_ref - ts_ref_mean
    ref_res_std = np.nanstd(ts_ref_res.values)
    ts_ref_mean_confInts = stats.t.interval(1 - alpha1, sample_num - 1, loc=ts_ref_mean.values,
                                            scale=np.sqrt(ref_res_std ** 2))
    #
    ts_pix_1_mean = ts_pix_1.rolling(time=sample_num, center=True).mean()
    ts_pix_1_res = ts_pix_1 - ts_pix_1_mean
    pix_1_res_std = np.nanstd(ts_pix_1_res.values)
    ts_pix_1_mean_confInts = stats.t.interval(1 - alpha1, sample_num - 1, loc=ts_pix_1_mean.values,
                                              scale=np.sqrt((pix_1_res_std ** 2) / sample_num))
    #
    ts_diff = ts_ref_mean - ts_pix_1_mean
    #
    # ----------------------------------------------------------------------
    # plotting the profiles and the location of the pixel
    # ----------------------------------------------------------------------
    start_date = np.datetime64('2016-09-01')
    end_date = np.datetime64('2020-01-01')
    datelist = pd.date_range(start=start_date, end=end_date, freq='YS').tolist()
    #
    gfc_cmap = cm.get_cmap('viridis', 20).reversed()
    gfc_cmap.colors[0][-1] = 0
    # ---------
    fig, (ax_ts1, ax_dif, ax_map) = plt.subplots(3, figsize=(8, 9), gridspec_kw={"height_ratios": (.2, .2, .6)})
    #
    ax_ts1.fill_between(ts_pix_1.coords['time'].values, ts_pix_1_mean_confInts[0], ts_pix_1_mean_confInts[1],
                        facecolor='lightblue', alpha=0.7)
    ts_pix_1.plot(ls='-', color='grey', lw=0.75, ax=ax_ts1)
    ts_pix_1_mean.plot(ls='-', color='blue', ax=ax_ts1)
    ax_ts1.fill_between(ts_ref.coords['time'].values, ts_ref_mean_confInts[0], ts_ref_mean_confInts[1],
                        facecolor='lightgreen', alpha=0.7)
    ts_ref_mean.plot(ls='-', color='green', ax=ax_ts1)
    #
    ax_ts1.set_xticks(datelist)
    ax_ts1.set_axisbelow(True)
    ax_ts1.set_yticks(np.arange(-20, -9, 1))
    ax_ts1.minorticks_on()
    ax_ts1.grid(which='major', linestyle='-', linewidth=1)
    ax_ts1.grid(which='minor', linestyle=':', linewidth=0.5)
    ax_ts1.legend(['single-pixels TS', 'single-pixels mean TS', 'reference mean TS'], loc="lower left")
    ax_ts1.set_title(r'Disturbed pixel index: {}'.format(dist_pix_ind))
    ax_ts1.set_xlabel('')
    ax_ts1.set_xlim('2016-09-01', '2020-03-01')
    plt.tight_layout()
    # -------------
    ts_diff.plot(ls='-', color='red', lw=1, ax=ax_dif)
    ax_dif.axhline(y=error_margin1, lw=0.75, ls='--', color='black')
    ax_dif.axhline(y=-error_margin1, lw=0.75, ls='--', color='black')
    #
    ax_dif.set_title('')
    ax_dif.set_xticks(datelist)
    ax_dif.set_axisbelow(True)
    # ax_dif.set_yticks(np.arange(-5, 5, 1))
    ax_dif.set_yticks(np.arange(np.floor(np.nanmin(ts_diff.values)), np.ceil(np.nanmax(ts_diff.values)) + 1, 1))
    ax_dif.minorticks_on()
    ax_dif.grid(which='major', linestyle='-', linewidth=1)
    ax_dif.grid(which='minor', linestyle=':', linewidth=0.5)
    ax_dif.set_xlabel('')
    ax_dif.legend(['difference TS-s', '95 % confidence margin'], loc="upper left")
    ax_dif.set_xlim('2016-09-01', '2020-03-01')
    # -------------
    ds.LossYear.where(ds.LossYear == LOSS_YEAR).plot(ax=ax_map, cmap=gfc_cmap, levels=np.arange(20), add_colorbar=False)
    plt.scatter(my_coor[0], my_coor[1], marker='o', s=20, color='red', axes=ax_map)
    ax_map.set_aspect('equal', 'box')
    ax_map.grid(linestyle=':', linewidth=0.5)
    aux1_x = my_coor[0]
    ax_map.set_xlim(aux1_x - 1500, aux1_x + 1500)
    aux1_y = my_coor[1]
    ax_map.set_ylim(aux1_y - 1500, aux1_y + 1500)
    plt.tight_layout()

def s1_ts_features(ts_pix_1):
    """
    This function calculates a number of features from S1-TS.

    The main segment to be analysed here is the 1st segment (in case of TS with multiple disturbance segments).

    :param  ts_pix_1: TS given as the xarray Data Array.
    :return: The output is a list with 18 parameters derived from the input TS.
    """
    #
    # hardcoded parameters:
    alpha = 0.99
    sample_num = 8
    segment_threshold = 8
    seg_id = 0
    # -------------------------------------------
    # import modules
    from scipy import stats
    import numpy as np
    from skimage import measure
    import numpy as xr
    # -------------------------------------------
    # define the default output array
    # there are currently 18 output parameters
    # --
    exception_label = 0
    ref_mean = np.nan
    error_margin = np.nan
    num_of_segments = np.nan
    TS_end_flag = np.nan
    TS_end_flag_long = 2
    TS_end_mag = np.nan
    # seg_id
    seg_size = np.nan
    max_mag = np.nan
    max_mag_date = np.nan
    t_pre = np.nan
    t_post = np.nan
    t_total = np.nan
    max_mag_org = np.nan
    max_mag_org_date = np.nan
    t_mag_org = np.nan
    seg2_size = np.nan
    # the output array
    my_output = [exception_label, ref_mean, error_margin,
                 num_of_segments, TS_end_flag, TS_end_flag_long, TS_end_mag,
                 seg_id, seg_size, max_mag, max_mag_date, t_pre, t_post, t_total,
                 max_mag_org, max_mag_org_date, t_mag_org,
                 seg2_size]
    # -------------------------------------------
    # check if the single-pixel TS has valid values at all:
    if np.all(np.isnan(ts_pix_1.values)):
        exception_label = 1
        my_output[0] = exception_label
        return my_output
    # ------------------ (a) --------------------
    # get the running mean of single-pixel TS:
    ts_pix_1_mean = ts_pix_1.rolling(time=sample_num, center=True).mean()
    ts_pix_1_mean = ts_pix_1_mean.dropna(dim='time')
    # ------------------ (b) --------------------
    # get the reference and its features:
    ref_mean = ts_pix_1_mean.sel(time=slice('2017-01-01', '2018-01-01')).values.mean()
    # get the std og the running mean:
    ref_std = ts_pix_1_mean.sel(time=slice('2017-01-01', '2018-01-01')).values.std()
    # ------------------ (c) -------------------
    # get the confidence intervals:
    [mean_DOWN, mean_UP] = stats.norm.interval(alpha, loc=ref_mean, scale=ref_std)
    error_margin = np.abs(mean_DOWN - ref_mean)
    # ------------------ (d) --------------------
    # get the difference between reference and single-pixel mean TS-s:
    ts_diff = ref_mean - ts_pix_1_mean
    # label dates that exceeds the threshold
    diff_bool = np.abs(ts_diff.values) > error_margin
    # ignore the values before the history period, i.e. the values in 2016 are set to False:
    diff_bool[ts_diff.time.dt.year.values == 2016] = False
    # ignore the values in the history period, i.e. the values in 2017 are set to False:
    diff_bool[ts_diff.time.dt.year.values == 2017] = False
    # ------------------ (e) --------------------
    # check if there is at least one element above the margin:
    if np.all(diff_bool == False):
        exception_label = 2
        my_output[0] = exception_label
        return my_output

    # ------------------ (f) --------------------
    # identify segments and get their properties
    # --
    # label segments in TS (consecutive labeled dates from the step "d")
    diff_bool_labels = measure.label(diff_bool)
    # as measure.regionprops only works for images, an extra auxiliary row is added here:
    diff_bool_labels_2D = np.array([diff_bool_labels, diff_bool_labels])
    # ---
    # this is new part added latter to tackle the case when the first segment has negative magnitude:
    # ---
    # get the max magnitude of the first segment:
    diff_bool_regions = measure.regionprops(diff_bool_labels_2D)
    aux_ts_diff_seg1 = ts_diff.where(diff_bool_labels == diff_bool_regions[0].label)
    aux_max_mag_ind = np.nanargmax(np.abs(aux_ts_diff_seg1.values))
    aux_max_mag = aux_ts_diff_seg1.values[aux_max_mag_ind]
    # only if the max_mag is positive, then segments with small separation should be connected
    if aux_max_mag > 0:
        # ------------------ (g) --------------------
        # check and meagre small intervals between the segments (the inverse segments)
        #---
        # get the inverse fo the segments boolean, i.e. the inverse segments
        inverse_bool = (diff_bool_labels_2D == 0)
        # label the inverse segments
        inverse_bool_labels = measure.label(inverse_bool)
        # get properties of the inverse segments
        inverse_bool_regions = measure.regionprops(inverse_bool_labels)
        # form a list with the lengths of the inverse segments
        seg_separation = [aux1.area/2 for aux1 in inverse_bool_regions]
        # form a list with labels of the inverse segments that corresponds to the length list
        separation_label = [aux1.label for aux1 in inverse_bool_regions]
        # check which inverse segments should be merged with their neighboring segments
        sep_label_ind = np.array(seg_separation) < segment_threshold
        # the first inverse segment should be false because it is the history period, but it is set here again, just in case
        sep_label_ind[0] = False
        # the last inverse segment should be "false" when the TS ends with this segment
        if diff_bool[-1] == False:
            sep_label_ind[-1] = False
        # get the labels of inverse segments to be merged
        sep_labels_to_meagre = np.array(separation_label)[sep_label_ind]
        # get the boolean flags for elements in diff_bool to be changed
        ind_to_change = np.in1d(inverse_bool_labels[0, :], sep_labels_to_meagre)
        # change the diff_bool elements, i.e. the resulting diff_bool_new will have connected segments
        diff_bool_new = diff_bool
        diff_bool_new[ind_to_change] = True
        # ------------------ (h) --------------------
        # the segments with small separation are now connected.
        # Thus, new segments and their properties should be calculated, i.e. the step "f" again for diff_bool_new
        # --
        # label new segments
        diff_bool_labels = measure.label(diff_bool_new)
        #  extra auxiliary row is added:
        diff_bool_labels_2D_new = np.array([diff_bool_labels, diff_bool_labels])
        # get the segments' properties
        diff_bool_regions = measure.regionprops(diff_bool_labels_2D_new)
    # ------------------ (i) --------------------
    # get the global TS properties from running mean TS:
    # ---
    # total number of segments:
    num_of_segments = diff_bool_labels.max()
    # check if the last value of the TS is within or outside the error_margin
    # TS_end_flag: (-1) below the reference, (0) in reference bounds, and (1) above the reference
    # TS_end_flag_long is set to 1 only when the majority of the last quarter is outside error_margin
    # TS_end_mag is the magnitude of the last elment of TS, or
    # the mean magnitude of the last quarter for TS_end_flag_long=1
    TS_end_mag = ts_diff.values[-1]
    if diff_bool[-1] == False:
        TS_end_flag = 0
        # check if all the elements in the last 3 months are smaller than error_margin
        if np.all(diff_bool[-segment_threshold:] == False):
            TS_end_flag_long = 0
    elif (diff_bool[-1] == True) and (TS_end_mag > 0) :
        TS_end_flag = -1
        # check if the majority in the the last 3 months is larger than error_margin
        if stats.mode( diff_bool[-segment_threshold:])[0][0]:
            TS_end_flag_long = 1
            TS_end_mag = np.nanmean(ts_diff.values[-1-segment_threshold:])
    elif (diff_bool[-1] == True) and (TS_end_mag < 0):
        TS_end_flag = 1
        if stats.mode( diff_bool[-segment_threshold:])[0][0]:
            TS_end_flag_long = 1
            TS_end_mag = np.nanmean(ts_diff.values[-1 - segment_threshold:])
    # ------------------ (j) --------------------
    # get the properties of the first segment from the running mean TS:
    # ---
    # derive the parameters for the first segment
    seg_size = diff_bool_regions[seg_id].area/2
    # get the indices of the start and the end of the largest segment:
    seg_ind_start = diff_bool_regions[seg_id].bbox[1]
    seg_ind_end = diff_bool_regions[seg_id].bbox[3] - 1
    # --- note ---
    # -1 above is because bbox gives the first index outside the segment, and for the segments at the end of TS that
    # returns an index outside the length of the TS
    # ------------
    # get the maximum disturbance magnitude in the segment
    ts_diff_seg1 = ts_diff.where(diff_bool_labels == diff_bool_regions[seg_id].label)
    max_mag_ind = np.nanargmax(np.abs(ts_diff_seg1.values))
    max_mag = ts_diff_seg1.values[max_mag_ind]
    max_mag_date = ts_diff_seg1.coords['time'].values[max_mag_ind]
    # get the np array with relative times (in days) from the maximum disturbance:
    d_time = (ts_diff.coords['time'].values - max_mag_date).astype('timedelta64[D]').astype(int)
    # get the disturbance parameters related to the start and the end of the largest segment:
    t_pre = d_time[seg_ind_start]
    t_post = d_time[seg_ind_end]
    t_total = np.abs(t_pre) + np.abs(t_post)
    # ---
    # size of the 2nd segment:
    if num_of_segments > 1:
        seg2_size = diff_bool_regions[seg_id+1].area/2
        # the division by 2 is due to the extra array row added above to trick 'measure.regionprops'
    else:
        seg2_size = np.nan
    # ------------------ (k) --------------------
    # get parameters from the original TS (no running mean applied):
    # ---
    # exclude the values from the biginging and end of the TS according to the running mean TS
    # this is important because previous indices should math also the original TS.
    ts_pix_1_short = ts_pix_1.sel(time=ts_pix_1_mean.coords['time'].values)
    # get the difference from the original TS
    diff_org = ref_mean - ts_pix_1_short
    # get the values with the first segment
    diff_org_seg1 = diff_org.where(diff_bool_labels == diff_bool_regions[seg_id].label)
    # get the maximum magnitude
    max_mag_org_ind = np.nanargmax(np.abs(diff_org_seg1.values))
    max_mag_org = diff_org_seg1.values[max_mag_org_ind]
    max_mag_org_date = diff_org_seg1.coords['time'].values[max_mag_org_ind]
    # the time difference between original max magnitude and the running mean max magnitude
    t_mag_org = d_time[max_mag_org_ind]
    # -------------------------------------------
    # the output array
    my_output = [exception_label, ref_mean, error_margin,
                 num_of_segments, TS_end_flag, TS_end_flag_long, TS_end_mag,
                 seg_id, seg_size, max_mag, max_mag_date, t_pre, t_post, t_total,
                 max_mag_org, max_mag_org_date, t_mag_org,
                 seg2_size]
    # -------------------------------------------
    return my_output




def s1_ts_features_ver2(ts_pix_1):
    """
    This function calculates a number of features from S1-TS.

    The function is an improved version of the function 's1_ts_features'. The improvements includes more robust
    definitions of certain TS-, and segment- features (see e.g. the 't_total' variable). Furthermore, here,
    the main segment to be analysed is one with the maximum disturbance magnitude, and NOT the 1st segment.
    The latter is the case in 's1_ts_features'. Then, the disturbance with the signal higher than the reference
    is ignored (see the definition of the 'diff_bool' variable in the code). The definition of the 'TS_end_flag_long'
    is more strict here compared to 's1_ts_features'.

    :param  ts_pix_1: TS given as the xarray Data Array.
    :return: The output is a list with 18 parameters derived from the input TS.
    """
    #
    # hardcoded parameters:
    alpha = 0.99
    sample_num = 15
    segment_threshold = 15
    seg_id = 0
    start_history = '2017-01-01'
    end_history = '2018-01-01'
    start_recovery = '2020-01-01'
    end_recovery = '2020-04-01'
    # -------------------------------------------
    # import modules
    from scipy import stats
    import numpy as np
    from skimage import measure
    import numpy as xr
    # -------------------------------------------
    # define the default output list
    # there are currently 18 output parameters
    # --
    exception_label = 0
    ref_mean = np.nan
    error_margin = np.nan
    num_of_segments = np.nan
    TS_end_flag = np.nan
    TS_end_flag_long = 2
    TS_end_mag = np.nan
    # seg_id
    seg_size = np.nan
    max_mag = np.nan
    max_mag_date = np.nan
    t_pre = np.nan
    t_post = np.nan
    t_total = np.nan
    max_mag_org = np.nan
    max_mag_org_date = np.nan
    t_mag_org = np.nan
    seg2_size = np.nan
    # the output array
    my_output = [exception_label, ref_mean, error_margin,
                 num_of_segments, TS_end_flag, TS_end_flag_long, TS_end_mag,
                 seg_id, seg_size, max_mag, max_mag_date, t_pre, t_post, t_total,
                 max_mag_org, max_mag_org_date, t_mag_org,
                 seg2_size]
    # -------------------------------------------
    # check if the single-pixel TS has valid values at all:
    if np.all(np.isnan(ts_pix_1.values)):
        exception_label = 1
        my_output[0] = exception_label
        return my_output
    # ------------------ (a) --------------------
    # get the running mean of single-pixel TS:
    ts_pix_1_mean = ts_pix_1.rolling(time=sample_num, center=True).mean()
    # remove the NaN-s
    ts_pix_1_mean = ts_pix_1_mean.dropna(dim='time')
    # remove everything after the recovery end period date
    ts_pix_1_mean = ts_pix_1_mean.sel(time=slice(start_history, end_recovery))
    # ------------------ (b) --------------------
    # get the reference and its features:
    ref_mean = ts_pix_1_mean.sel(time=slice(start_history, end_history)).values.mean()
    # get the std og the running mean:
    ref_std = ts_pix_1_mean.sel(time=slice(start_history, end_history)).values.std()
    # ------------------ (c) -------------------
    # get the confidence intervals:
    [mean_DOWN, mean_UP] = stats.norm.interval(alpha, loc=ref_mean, scale=ref_std)
    error_margin = np.abs(mean_DOWN - ref_mean)
    # ------------------ (d) --------------------
    # get the difference between reference and single-pixel mean TS-s:
    ts_diff = ref_mean - ts_pix_1_mean
    # label dates that exceeds the threshold
    # !!! NOTE !!! the disturbance with the signal above the reference level is ignored
    diff_bool = ts_pix_1_mean < mean_DOWN
    # ignore the values before the history period, i.e. the values in 2016 are set to False:
    diff_bool[ts_diff.time.dt.year.values == 2016] = False
    # ignore the values in the history period, i.e. the values in 2017 are set to False:
    diff_bool[ts_diff.time.dt.year.values == 2017] = False
    # ------------------ (e) --------------------
    # check if there is at least one element above the margin:
    if np.all(diff_bool.values == False):
        exception_label = 2
        my_output[0] = exception_label
        return my_output

    # ------------------ (f) --------------------
    # identify segments and get their properties
    # --
    # label segments in TS (consecutive labeled dates from the step "d")
    diff_bool_labels = measure.label(diff_bool)
    # as measure.regionprops only works for images, an extra auxiliary row is added here:
    diff_bool_labels_2D = np.array([diff_bool_labels, diff_bool_labels])
    # get the regions, i.e. the segments of TS where disturbance happens:
    diff_bool_regions = measure.regionprops(diff_bool_labels_2D)
    # ---
    # this is new part added to tackle the case when there are more than one segment:
    # ---
    # total number of segments:
    num_of_segments = diff_bool_labels.max()
    #
    if num_of_segments > 1:
        # ------------------ (g) --------------------
        # check and meagre small intervals between the segments (the inverse segments)
        #---
        # get the inverse fo the segments boolean, i.e. the inverse segments
        inverse_bool = (diff_bool_labels_2D == 0)
        # label the inverse segments
        inverse_bool_labels = measure.label(inverse_bool)
        # get properties of the inverse segments
        inverse_bool_regions = measure.regionprops(inverse_bool_labels)
        # form a list with the lengths of the inverse segments
        seg_separation = [aux1.area/2 for aux1 in inverse_bool_regions]
        # form a list with labels of the inverse segments that corresponds to the length list
        separation_label = [aux1.label for aux1 in inverse_bool_regions]
        # check which inverse segments should be merged with their neighboring segments
        sep_label_ind = np.array(seg_separation) <= segment_threshold
        # the first inverse segment should be false because it is the history period,
        # but it is set here again, just in case
        sep_label_ind[0] = False
        # the last inverse segment should be "false" when the TS ends with this segment
        if diff_bool.values[-1] == False:
            sep_label_ind[-1] = False
        # get the labels of inverse segments to be merged
        sep_labels_to_meagre = np.array(separation_label)[sep_label_ind]
        # get the boolean flags for elements in diff_bool to be changed
        ind_to_change = np.in1d(inverse_bool_labels[0, :], sep_labels_to_meagre)
        # change the diff_bool elements, i.e. the resulting diff_bool_new will have connected segments
        diff_bool_new = diff_bool
        diff_bool_new[ind_to_change] = True
        # ------------------ (h) --------------------
        # the segments with small separation are now connected.
        # Thus, new segments and their properties should be calculated, i.e. the step "f" again for diff_bool_new
        # --
        # label new segments
        diff_bool_labels = measure.label(diff_bool_new)
        #  extra auxiliary row is added:
        diff_bool_labels_2D_new = np.array([diff_bool_labels, diff_bool_labels])
        # get the segments' properties
        diff_bool_regions = measure.regionprops(diff_bool_labels_2D_new)
        # update the total number of segments:
        num_of_segments = diff_bool_labels.max()
    # ------------------ (i) --------------------
    # get the global TS properties from running mean TS:
    # ---
    # check if the last value of the TS is within or outside the error_margin
    # TS_end_flag: (-1) below the reference, (0) in reference bounds, and (1) above the reference
    # TS_end_flag_long is set to 1 when one, or more observations in the Q4 is outside the error_margin
    # TS_end_mag is the magnitude of the last observation of TS, or
    # the mean magnitude of the last quarter for TS_end_flag_long=1
    TS_end_mag = ts_diff.values[-1]
    if diff_bool.values[-1] == True:
        TS_end_flag = np.int(np.sign(TS_end_mag))
        TS_end_flag_long = 1
    elif diff_bool.values[-1] == False:
        TS_end_flag = 0
        # check if at least one element in the recovery period is outside the error margin:
        if np.any(diff_bool.sel(time=slice(start_recovery, end_recovery)).values):
            TS_end_flag_long = 1
        else:
            TS_end_flag_long = 0
    # --------------------------------------------
    # new part added to ensure in case of multiple segments that we analyse the segment with the largest disturbance
    # --------------------------------------------
    if num_of_segments > 1:
        # update the segment id value:
        seg_id = diff_bool_labels[np.argmax(ts_diff.values)] - 1
    # ------------------ (j) --------------------
    # get the properties of the selected segment from the running mean TS:
    # ---
    # derive the parameters for the selected segment
    seg_size = diff_bool_regions[seg_id].area/2
    # get the indices of the start and the end of the selected segment:
    seg_ind_start = diff_bool_regions[seg_id].bbox[1]
    seg_ind_end = diff_bool_regions[seg_id].bbox[3] - 1
    # --- note ---
    # -1 above is because bbox gives the first index outside the segment, and for the segments at the end of TS that
    # returns an index outside the length of the TS
    # ------------
    # get the maximum disturbance magnitude in the segment
    ts_diff_seg1 = ts_diff.where(diff_bool_labels == diff_bool_regions[seg_id].label)
    max_mag_ind = np.nanargmax(np.abs(ts_diff_seg1.values))
    max_mag = ts_diff_seg1.values[max_mag_ind]
    max_mag_date = ts_diff_seg1.coords['time'].values[max_mag_ind]
    # get the np array with relative times (in days) from the maximum disturbance:
    d_time = (ts_diff.coords['time'].values - max_mag_date).astype('timedelta64[D]').astype(int)
    # get the disturbance parameters related to the start and the end of the largest segment:
    t_pre = d_time[seg_ind_start]
    t_post = d_time[seg_ind_end]
    # !!! NOTE !!! this part is new. total time now considers multiple segments
    # in that case, the total time is between tthe start of the first segment and the end of the last segment
    if num_of_segments > 1:
        first_seg_start_ind = diff_bool_regions[0].bbox[1]
        last_seg_end_ind = diff_bool_regions[-1].bbox[3] - 1
        t_total = (ts_diff.coords['time'].values[last_seg_end_ind] -
                   ts_diff.coords['time'].values[first_seg_start_ind])\
            .astype('timedelta64[D]').astype(int)
    else:
        t_total = np.abs(t_pre) + np.abs(t_post)
    # ---
    # size of another segment (either one after the max mag segment, or the first one):
    if num_of_segments > 1:
        if seg_id == 0:
            seg2_size = diff_bool_regions[seg_id+1].area/2
            # the division by 2 is due to the extra array row added above to trick 'measure.regionprops'
        else:
            seg2_size = diff_bool_regions[0].area / 2
    else:
        seg2_size = np.nan
    # ------------------ (k) --------------------
    # get parameters from the original TS (no running mean applied):
    # ---
    # exclude the values from the beginning and end of the TS according to the running mean TS
    # this is important because previous indices should math also the original TS.
    ts_pix_1_short = ts_pix_1.sel(time=ts_pix_1_mean.coords['time'].values)
    # get the difference from the original TS
    diff_org = ref_mean - ts_pix_1_short
    # get the values with the first segment
    diff_org_seg1 = diff_org.where(diff_bool_labels == diff_bool_regions[seg_id].label)
    # get the maximum magnitude
    max_mag_org_ind = np.nanargmax(np.abs(diff_org_seg1.values))
    max_mag_org = diff_org_seg1.values[max_mag_org_ind]
    max_mag_org_date = diff_org_seg1.coords['time'].values[max_mag_org_ind]
    # the time difference between original max magnitude and the running mean max magnitude
    t_mag_org = d_time[max_mag_org_ind]
    # -------------------------------------------
    # the output array
    my_output = [exception_label, ref_mean, error_margin,
                 num_of_segments, TS_end_flag, TS_end_flag_long, TS_end_mag,
                 seg_id, seg_size, max_mag, max_mag_date, t_pre, t_post, t_total,
                 max_mag_org, max_mag_org_date, t_mag_org,
                 seg2_size]
    # -------------------------------------------
    return my_output


def plot_TS(ts_pix_1, POLARIZATION):
    #---------------------
    import numpy as np
    from matplotlib import pyplot as plt
    from matplotlib import cm
    import pandas as pd
    from scipy import stats
    # ------------------
    # define the history time:
    start_history = '2017-01-01'
    end_history = '2018-01-01'
    # define some plotting parameters
    start_plot = start_history
    end_plot = '2020-04-01'
    # specify the sample number (3 months = 90 days / 6-day = 15):
    sample_num = 15
    # ------------------
    # calculate end plot time
    # ------------------
    ts_pix_1_mean = ts_pix_1.rolling(time=sample_num, center=True).mean()
    # plotting the TSs, the difference TS and the spatial location of the pixel
    my_coor = [ts_pix_1.coords['x'].values.flatten()[0], ts_pix_1.coords['y'].values.flatten()[0]]
    # get the reference:
    ref_mean = ts_pix_1_mean.sel(time=slice(start_history, end_history)).values.mean()
    # get the std og the running mean:
    ref_std = ts_pix_1_mean.sel(time=slice(start_history, end_history)).values.std()
    # get the confidence intervals:
    [mean_DOWN, mean_UP] = stats.norm.interval(0.99, loc=ref_mean, scale=ref_std)
    total_num = ts_pix_1.count().values
    # ---------------------------------------------
    # plotting
    # ---------------------------------------------
    datelist = pd.date_range(start=np.datetime64(start_plot), end=np.datetime64(end_plot), freq='YS').tolist()
    #
    gfc_cmap = cm.get_cmap('viridis', 20).reversed()
    gfc_cmap.colors[0][-1] = 0
    # -------------
    fig, ax_ts = plt.subplots(1, figsize=(8, 3))
    ts_pix_1.plot(color='grey', ax=ax_ts)
    ts_pix_1_mean.plot(ls='-', color='blue', ax=ax_ts)
    ax_ts.fill_between(ts_pix_1.coords['time'].values, mean_DOWN*np.ones(total_num), mean_UP*np.ones(total_num),
                            facecolor='lightgreen', alpha=0.7)
    ax_ts.axhline(y=ref_mean, lw=2, ls='-', color='green')
    ax_ts.set_xticks(datelist)
    ax_ts.set_axisbelow(True)
    if POLARIZATION == 'VH':
        ax_ts.set_yticks(np.arange(-20, -9, 1))
    else:
        ax_ts.set_yticks(np.arange(-14, -3, 1))
    ax_ts.minorticks_on()
    ax_ts.grid(which='major', linestyle='-', linewidth=1)
    ax_ts.grid(which='minor', linestyle=':', linewidth=0.5)
    ax_ts.legend(['original TS', 'running mean TS', 'reference line', '99 % confidence'], loc="lower left")
    ax_ts.set_title(r'{} Pol., Location [m]: x={}, y={} '.format(POLARIZATION,
                                                                         np.int(my_coor[0]),
                                                                         np.int(my_coor[1])))
    ax_ts.set_xlim((start_plot, end_plot))
    plt.tight_layout()


def plot_2D_map(my_coor, ds, ds_var_name):
    # ---------------------
    import numpy as np
    from matplotlib import pyplot as plt
    from matplotlib import cm
    import pandas as pd
    from scipy import stats
    # ------------------
    if ds_var_name == 'LossYear':
        # ---------------
        # define some plotting parameters
        gfc_cmap = cm.get_cmap('viridis', 20).reversed()
        gfc_cmap.colors[0][-1] = 0
    else:
        gfc_cmap = cm.get_cmap('viridis', 20)
    # -------------
    fig, ax_map = plt.subplots(1, figsize=(6, 6))
    if ds_var_name == 'LossYear':
        ds.get(ds_var_name).where(ds.LossYear == 18).plot(ax=ax_map, cmap=gfc_cmap, levels=np.arange(20), add_colorbar=False)
    else:
        ds.get(ds_var_name).where(ds.LossYear == 18).plot(ax=ax_map, cmap=gfc_cmap, levels=np.arange(20),
                                                          vmax=np.nanquantile(ds.get(ds_var_name).values,0.98),
                                                          vmin=np.nanquantile(ds.get(ds_var_name).values,0.02),
                                                          )
    plt.scatter(my_coor[0], my_coor[1], marker='o', s=20, color='red', axes=ax_map)
    ax_map.set_aspect('equal', 'box')
    ax_map.grid(linestyle=':', linewidth=0.5)
    aux1_x = my_coor[0]
    ax_map.set_xlim(aux1_x - 1500, aux1_x + 1500)
    aux1_y = my_coor[1]
    ax_map.set_ylim(aux1_y - 1500, aux1_y + 1500)
    plt.tight_layout()

def ts_feature2ds(ds_disturb, featName, dist_pix_par_all, disturbed_pixels):
    """
    The function assigns the feature values as new variable in the xarray Dataset
    """
    # -----
    import pandas as pd
    import numpy as np
    # -----
    featNames = ['exception_label', 'ref_mean', 'error_margin', 'num_of_segments', 'TS_end_flag',
                 'TS_end_flag_long',
                 'TS_end_mag', 'seg_id', 'seg_size', 'max_mag', 'max_mag_date', 't_pre', 't_post', 't_total',
                 'max_mag_org', 'max_mag_org_date', 't_mag_org', 'seg2_size']
    # -----
    if featNames.index(featName) == 10 or featNames.index(featName) == 15:
        myFeature_year = [pd.to_datetime(aux1[featNames.index(featName)]).year for aux1 in dist_pix_par_all]
        myFeature_DOY = [pd.to_datetime(aux1[featNames.index(featName)]).dayofyear for aux1 in dist_pix_par_all]
        # transform variables as 2D-arrays:
        myFeature_year_array = np.nan * np.ones([ds_disturb.dims['y'], ds_disturb.dims['x']])
        myFeature_year_array[disturbed_pixels[:, 0], disturbed_pixels[:, 1]] = myFeature_year
        #
        myFeature_DOY_array = np.nan * np.ones([ds_disturb.dims['y'], ds_disturb.dims['x']])
        myFeature_DOY_array[disturbed_pixels[:, 0], disturbed_pixels[:, 1]] = myFeature_DOY
        # assign new varibles:
        ds_disturb[featName + '_YEAR'] = (['y', 'x'], myFeature_year_array)
        ds_disturb[featName + '_DOY'] = (['y', 'x'], myFeature_DOY_array)
    else:
        myFeature = [aux1[featNames.index(featName)] for aux1 in dist_pix_par_all]
        # transform variables as 2D-arrays:
        myFeature_array = np.nan * np.ones([ds_disturb.dims['y'], ds_disturb.dims['x']])
        myFeature_array[disturbed_pixels[:, 0], disturbed_pixels[:, 1]] = myFeature
        # assign new varibles:
        ds_disturb[featName] = (['y', 'x'], myFeature_array)


def ts_categoryLabel2ds(ds_disturb, ds_Distrb_std, th_std):
    """
    The function assigns a category label (8 TS categories in total) for each disturbed pixel
    as a new coordinate in the input xarray Dataset
    Labels:
    0 - not disturbed in 2018 according to the Hansen 'lossyear' layer
    1 - TS with all nan values (usually border pixels)
    2 - TS with single disturb. segment and fully recovered in the Q4 of 2019
    3 - disturbance in the history period (2017)
    4 - no change compared to the history period
    5 - disturbance with negative magnitude (disturbance backscatter larger than the reference)
    6 - disturbance with positive magnitude, but the signal never recovered in the Q4 of 2019
    7 - TS with positive disturbance, but with multiple segments
    -----
    Note: only labels in 'ds_disturb.pixel_labels_v4' contains all the labels
    """
    # -------------------
    import numpy as np
    # ------ round 1 ----------------------------------------------------------
    # pixels with all nan values in TS
    nanTS_mask = 1 * (ds_disturb.exception_label.values == 1)
    # typical TS (exceed the error margin)
    ChangeTS_mask = 2 * (ds_disturb.exception_label.values == 0)
    # Label the no-change to large and small Std
    ds_Distrb_std_noChange = ds_Distrb_std.where(ds_disturb.exception_label == 2)
    # label pixels with TS disturbed in history period:
    largeStdTS_mask = 3 * (ds_Distrb_std_noChange.values > th_std)
    # no-change TS (does not exceed the error margin)
    noChangeTS_mask = 4 * (ds_Distrb_std_noChange.values <= th_std)
    # sum up all the masks
    pix_labels_array = nanTS_mask + ChangeTS_mask + largeStdTS_mask + noChangeTS_mask
    # assign the mask to the dataset:
    ds_disturb.coords['pixel_labels'] = (('y', 'x'), pix_labels_array)
    # ------ round 2 ----------------------------------------------------------
    # ----------------------------------------------------------------------
    #  Disturbance Magnitude - classify pixels with positive and negative values
    # ----------------------------------------------------------------------
    ds_typical_TS = ds_disturb.where(ds_disturb.pixel_labels == 2)
    # split the class 2 (typical disturbance) to pixels with negative magnitude:
    neg_mag_mask = 5 * (ds_typical_TS.max_mag.values < 0)
    # split the class 2 (typical disturbance) to pixels with positive magnitude:
    pos_mag_mask = 2 * (ds_typical_TS.max_mag.values > 0)
    # ----
    # calculate new labels
    # ----
    # set the class 2 to zero as this class is now splited in two more:
    previous_labels = np.zeros(pix_labels_array.shape) + pix_labels_array
    previous_labels[pix_labels_array == 2] = 0
    # sum up all the masks:
    pix_labels_array_v2 = previous_labels + neg_mag_mask + pos_mag_mask
    # assign the mask to the dataset:
    ds_disturb.coords['pixel_labels_v2'] = (('y', 'x'), pix_labels_array_v2)
    # ------ round 3 ----------------------------------------------------------
    # ----------------------------------------------------------------------
    #  Multiple segments -  further classification
    # ----------------------------------------------------------------------
    ds_typical_TS_v2 = ds_disturb.where(ds_disturb.pixel_labels_v2 == 2)
    # ---------
    # split the class 2 (typical disturbance) to pixels with multiple segments
    multiSegmnts_mask = 6 * (ds_typical_TS_v2.num_of_segments.values > 1)
    # split the class 2 (typical disturbance) to pixels with single disturbance segment:
    singleSegmnt_mask = 2 * (ds_typical_TS_v2.num_of_segments.values == 1)
    # ----
    # calculate new labels
    # ----
    # set the class 2 to zero as this class is now splited in two more:
    previous_labels = np.zeros(pix_labels_array.shape) + pix_labels_array_v2
    previous_labels[pix_labels_array_v2 == 2] = 0
    # sum up all the masks:
    pix_labels_array_v3 = previous_labels + multiSegmnts_mask + singleSegmnt_mask
    # assign the mask to the dataset:
    ds_disturb.coords['pixel_labels_v3'] = (('y', 'x'), pix_labels_array_v3)
    # ------ round 4 ----------------------------------------------------------
    # ----------------------------------------------------------------------
    #  Positive Disturbance Magnitude
    # ----------------------------------------------------------------------
    ds_typical_TS_v3 = ds_disturb.where(ds_disturb.pixel_labels_v3 == 2)
    # ---------
    #  Classify pixels with no-recovery, and typical recovery
    # ---------
    # split the class 2 (typical single-segment disturbance) to pixels with positive disturbance that never recover
    posMag_noRecovery_mask = 7 * (ds_typical_TS_v3.TS_end_flag_long.values > 0)
    # split the class 2 (typical disturbance) to pixels with positive magnitude that recover in the last quoter:
    posMag_yesRecovery_mask = 2 * (ds_typical_TS_v3.TS_end_flag_long.values == 0)
    # ----
    # calculate new labels
    # ----
    # set the class 2 to zero as this class is now splited in two more:
    previous_labels = np.zeros(pix_labels_array.shape) + pix_labels_array_v3
    previous_labels[pix_labels_array_v3 == 2] = 0
    # sum up all the masks:
    pix_labels_array_v4 = previous_labels + posMag_noRecovery_mask + posMag_yesRecovery_mask
    # assign the mask to the dataset:
    ds_disturb.coords['pixel_labels_v4'] = (('y', 'x'), pix_labels_array_v4)


def ts_categoryLabel2ds_ver2(ds_disturb, ds_Distrb_std, th_std):
    """
    The function assigns a category label (6 TS categories in total) for each disturbed pixel
    as a new coordinate with the name 'CategoryLabel' in the input xarray Dataset.

    The function is an modified version of the function 'ts_categoryLabel2ds', and the TS features used here
    should be calculated using the function 's1_ts_features_ver2'.

    !!! NOTE !!! old (deprecated) names for the new coordinate were: 'pixel_labels_v4', or 'pixel_labels'. Those
    names should be replaced with the actual coordinate name, i.e. 'CategoryLabel'.

    Output category labels:

    0 - not disturbed in 2018 according to the Hansen 'lossyear' layer

    1 - TS with all nan values (usually border pixels)

    2 - Disturbed and fully recovered by the Q4 of 2019 period

    3 - Disturbed and not recovered by the Q4 of 2019 period

    4 - no change because of a disturbance in the history period (2017)

    5 - no change compared to the history period

    :param ds_disturb: the input xarray dataset with the TS features calculated with the 's1_ts_features_ver2' function,
     and then assigned using the 'ts_feature2ds' function
    :param ds_Distrb_std: xarray data array with the Std. calculated for the running mean TS in the history period
    (dim. 'time' not exist !!!)
    :param th_std: the Std. value to threshold 'ds_Distrb_std'
    :return: the input xarray 'ds_disturb' with a new coordinate 'CategoryLabel' that contains the labels
    """
    # -------------------
    import numpy as np
    # ------ round 1 ----------------------------------------------------------
    # pixels with all nan values in TS
    nanTS_mask = 1*(ds_disturb.exception_label.values == 1)
    # typical TS (exceed the error margin)
    ChangeTS_mask = 2*(ds_disturb.exception_label.values == 0)
    # Label the no-change to large and small Std
    ds_Distrb_std_noChange = ds_Distrb_std.where(ds_disturb.exception_label == 2)
    # label pixels with TS disturbed in history period:
    largeStdTS_mask = 3*(ds_Distrb_std_noChange.values > th_std)
    # no-change TS (does not exceed the error margin)
    noChangeTS_mask = 4*(ds_Distrb_std_noChange.values <= th_std)
    # sum up all the masks
    pix_labels_array = nanTS_mask + ChangeTS_mask + largeStdTS_mask + noChangeTS_mask
    # assign the mask to the dataset:
    ds_disturb.coords['pixel_labels'] = (('y', 'x'), pix_labels_array)
    # ds_disturb.coords['pixel_labels'] = (('y', 'x'), pix_labels_array) # deprecated name 'pixel_labels'
    # ------ round 2 ----------------------------------------------------------
    ds_typical_TS_v3 = ds_disturb.where(ds_disturb.pixel_labels == 2)
    # ---------
    #  Classify pixels with no-recovery, and typical recovery
    # ---------
    # split the class 2 to pixels that never recovered
    posMag_noRecovery_mask = 5*(ds_typical_TS_v3.TS_end_flag_long.values > 0)
    # split the class 2 to pixels that fully recovered by the Q4:
    posMag_yesRecovery_mask = 2*(ds_typical_TS_v3.TS_end_flag_long.values == 0)
    # ----
    # calculate new labels
    # ----
    # set the class 2 to zero as this class is now split in two more:
    previous_labels = np.zeros(pix_labels_array.shape) + pix_labels_array
    previous_labels[pix_labels_array == 2] = 0
    # sum up all the masks:
    pix_labels_array_v4 = previous_labels + posMag_noRecovery_mask + posMag_yesRecovery_mask
    # assign the mask to the dataset:
    ds_disturb.coords['CategoryLabel'] = (('y', 'x'), pix_labels_array_v4)
    # ds_disturb.coords['pixel_labels_v4'] = (('y', 'x'), pix_labels_array_v4) # deprecated name 'pixel_labels_v4'

def aggregate2reg6DayCube(ds):
    """
    Aggregate the input data structure to regular 6-day data structure.
    New times are defined as 6-day multiples to the time of the first image that has more than 98 % of valid pixels.
    The aggregation function is hardcoded as the mean of all images in the 6-day interval.
    """
    import numpy as np
    import pandas as pd
    #
    # get the number of valid pixels per image:
    ds_count = ds.VH.groupby('time').count(dim=['x', 'y'])
    # select only images that cover 98 % of tile or more:
    ds_fullCover = ds_count.where(ds_count > 0.98*ds_count.values.max(), drop=True)
    # get the start time of the bin with the first full cover image:
    FirstFullCover_BinTime = ds_fullCover.time.values[0] - np.timedelta64(3, 'D')
    # if the ds times start earlier than the time of the first full cover image, then correct the startingBinTime
    # otherwise startingBinTime is identical to the start time of the first bin with full cover image
    if FirstFullCover_BinTime >= ds.time.values[0]:
        startingBinTime = pd.date_range(FirstFullCover_BinTime, ds.time.values[0], freq='-6D')[-1]
    else:
        startingBinTime = FirstFullCover_BinTime
    # specify all the start-end times for the bins and their labels for the resampling and aggregation
    my_bins = pd.date_range(startingBinTime, ds.time.values[-1], freq='6D')
    my_labels = pd.date_range(startingBinTime + np.timedelta64(3, 'D'), my_bins[-1], freq='6D')
    # resample to the cube with 6 days time step (aggregated with the mean)
    ds2 = ds.groupby_bins('time', bins=my_bins, labels=my_labels, right=True, include_lowest=True).mean()
    ds2 = ds2.rename({'time_bins': 'time'})
    # interpolate linearly the NaN values
    ds3 = ds2.interpolate_na(dim='time', method='linear')
    # !!! NOTE !! the linear na interpolation may introduce NaN-s in last or first image when those contain NaN pixels
    # Thus, check if there are NaN-s in the last image, and if so, set those NaN-s with the values of the nearest image
    nanPix_NumLast = np.isnan(ds3.VH[-1, :, :].values).flatten().sum()
    if nanPix_NumLast > 0:
        ImageLast = ds3.isel(time=-1)
        Images2ndLast = ds3.isel(time=-2)
        newImageLast = ImageLast.where(~ImageLast.VH.isnull(), Images2ndLast.where(ImageLast.isnull()))
        ds3.VH[-1, :, :] = newImageLast.VH
        ds3.VV[-1, :, :] = newImageLast.VV
    # do the same for the first image:
    nanPix_NumFirst = np.isnan(ds3.VH[0, :, :].values).flatten().sum()
    if nanPix_NumFirst > 0:
        ImageFirst = ds3.isel(time=0)
        Images2nd = ds3.isel(time=1)
        newImageFirst = ImageFirst.where(~ImageFirst.VH.isnull(), Images2nd.where(ImageFirst.isnull()))
        ds3.VH[0, :, :] = newImageFirst.VH
        ds3.VV[0, :, :] = newImageFirst.VV
    return ds3
