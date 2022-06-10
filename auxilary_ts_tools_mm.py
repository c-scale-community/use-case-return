# Author:       Milutin Milenkovic
# Created:      13/05/2022
# Copyright:    ???
# Licence:      ???


def features_as_xrrray_ufunc(my_ndArray, ts_time_stamps):
    # reduce the 3d array (x, y, time) to a2d array (x, y) applying an user function:
    import numpy as np
    my_2d_array = np.apply_along_axis(wrapper_4_np1dArray, axis=-1, arr=my_ndArray, ts_time_stamps=ts_time_stamps)
    return my_2d_array

def wrapper_4_np1dArray(ts_aaray, ts_time_stamps):
    #  make a panda seres from the input 
    import pandas as pd
    ts_pix_1 = pd.Series(data=ts_aaray, index=ts_time_stamps).dropna()
    # deal with the partly identical index
    oidx = ts_pix_1.index
    nidx = pd.date_range(oidx.min(), oidx.max(), freq='6D')
    # interpolate into 6d time seres:
    #my_6D_ts = ts_pix_1.resample('6D').interpolate(method='linear')
    my_6D_ts = ts_pix_1.reindex(oidx.union(nidx)).interpolate(method='linear').round(1)
    # calculate the time seres features:
    myFeatures = features_from_S1_TS(my_6D_ts)
    return myFeatures


def features_from_S1_TS(ts_pix_1):
    """
    This function calculates a number of features from a Sentinel-1 time seres (TS).

    The function analyses the disturbances in the regular (6 day) time seres (TS). The disturbances are detected when 
    the 3-mont mean (hardcoded) TS exceeds its 99% confidence intervals. The confidence intervals were drawn relatively 
    from the first year of the TS treated as the history period. The last year of TS is treated as the recovery period.
    
    The TS analysis is performed segment-wise and the disturbance features are derived from: 
    (a) the segment with the largest disturbance magnitude, 
    (b) from the folowing or the first segment, and 
    (c) from the original (non-smoothed) TS
    
    
    !!! NOTE !!! ------------
    - The total disturbance time (the 't_total' variable) is here mesured between the start of the first segment and end of
    the last segment, provided that a TS has multiple segments. 
    - The disturbance with the signal higher than the reference is labeled with the class 3 and is not processed further.
    - The start and end of TS are here hardcoded and should be adopted directly in the code

    :param  ts_pix_1: TS given as the pandas Series.
    :return: The output is a list with 17 parameters derived from the input TS.
    """
    #
    # hardcoded parameters:
    alpha = 0.99
    sample_num = 15
    segment_threshold = 15
    seg_id = 0
    #
    start_history = '2017-01-01'
    end_history = '2018-01-01'
    start_recovery = '2020-01-01'
    end_recovery = '2021-01-01'
    # specify what is the time interval of the inpunt seres in days
    ts_dTime = 6
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
    # check if the TS has valid values at all:
    if np.all(np.isnan(ts_pix_1.values)):
        exception_label = 1
        my_output[0] = exception_label
        return my_output
    # ------------------ (a) --------------------
    # get the running mean of the TS:
    ts_pix_1_mean = ts_pix_1.rolling(sample_num, center=True).mean()
    # remove the NaN-s
    ts_pix_1_mean = ts_pix_1_mean.dropna()
    # remove everything after the recovery end period date
    ts_pix_1_mean = ts_pix_1_mean.loc[slice(start_history, end_recovery)]
    # ------------------ (b) --------------------
    # get the reference and its features:
    ref_mean = ts_pix_1_mean.loc[slice(start_history, end_history)].values.mean()
    # get the std og the running mean:
    ref_std = ts_pix_1_mean.loc[slice(start_history, end_history)].values.std()
    # ------------------ (c) -------------------
    # get the confidence intervals:
    [mean_DOWN, mean_UP] = stats.norm.interval(alpha, loc=ref_mean, scale=ref_std)
    error_margin = np.abs(mean_DOWN - ref_mean)
    # ------------------ (d) --------------------
    # get the difference between reference and the mean TS-s:
    ts_diff = ref_mean - ts_pix_1_mean
    # label dates that exceeds the threshold
    # !!! NOTE !!! the disturbance with the signal above the reference level is ignored
    diff_bool = ts_pix_1_mean < mean_DOWN
    # ignore the values before the history period, i.e. the values in 2017 or earler are set to False:
    diff_bool[ts_diff.index.year.values <= 2017] = False
    # ignore the values in the history period, i.e. the values in 2018 are set to False:
    diff_bool[ts_diff.index.year.values == 2018] = False
    # ------------------ (e) --------------------
    # check if the whole mean TS is above the DOWN margin, i.e., no "typical" disturbance:
    if np.all(diff_bool.values == False):
        # further check if any magnitude is above the UP bounds:
        diff_bool_UP = ts_pix_1_mean > mean_UP
        if np.any(diff_bool_UP.values == True):
            # label this TS as the "atypical disturbance"
            exception_label = 3
        else:
            # label this TS as "NO disturbance"
            exception_label = 2
            
        # update the output   
        my_output[0] = exception_label
        my_output[1] = ref_mean
        my_output[2] = error_margin
        return my_output
    # ------------------ (f) --------------------
    # identify segments and get their properties
    # --
    # label segments in TS (consecutive labeled dates from the step "d")
    diff_bool_labels = measure.label(diff_bool.values)
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
        # get the inverse for the segments boolean, i.e. the inverse segments
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
        diff_bool_labels = measure.label(diff_bool_new.values)
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
    # TS_end_flag_long is set to 1 when one, or more observations in the recovery perod is outside the error_margin
    # TS_end_mag is the magnitude of the last observation of TS, or
    # the mean magnitude of the last quarter for TS_end_flag_long=1
    TS_end_mag = ts_diff.values[-1]
    if diff_bool.values[-1] == True:
        TS_end_flag = np.int(np.sign(TS_end_mag))
        TS_end_flag_long = 1
    elif diff_bool.values[-1] == False:
        TS_end_flag = 0
        # check if at least one element in the recovery period is outside the error margin:
        if np.any(diff_bool.loc[slice(start_recovery, end_recovery)].values):
            TS_end_flag_long = 1
        else:
            TS_end_flag_long = 0
    # --------------------------------------------
    # new part added to ensure in case of multiple segments that we analyse the segment with the largest disturbance
    # --------------------------------------------
    if num_of_segments > 1:
        # update the segment id value:
        seg_id = diff_bool_labels[np.argmax(ts_diff.values)] - 1
        #seg_id = seg_id[0]
    # ------------------ (j) --------------------
    # get the properties of the selected segment from the running mean TS:
    # ---
    # derive the parameters for the selected segment
    # and multipy by the time interval to get the segment size in days
    seg_size = (diff_bool_regions[seg_id].area/2)*ts_dTime
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
    max_mag_date = ts_diff_seg1.index.values[max_mag_ind]
    # get the np array with relative times (in days) from the maximum disturbance:
    d_time = (ts_diff.index.values - max_mag_date).astype('timedelta64[D]').astype(int)
    # get the disturbance parameters related to the start and the end of the largest segment:
    t_pre = d_time[seg_ind_start]
    t_post = d_time[seg_ind_end]
    # !!! NOTE !!! this part is new. total time now considers multiple segments
    # in that case, the total time is between the start of the first segment and the end of the last segment
    if num_of_segments > 1:
        first_seg_start_ind = diff_bool_regions[0].bbox[1]
        last_seg_end_ind = diff_bool_regions[-1].bbox[3] - 1
        t_total = (ts_diff.index.values[last_seg_end_ind] -
                   ts_diff.index.values[first_seg_start_ind])\
            .astype('timedelta64[D]').astype(int)
    else:
        # add one segment that refers to max mag as this segment is neither pre- or post-disturbance
        t_total = np.abs(t_pre) + np.abs(t_post) + ts_dTime
    # ---
    # size of another segment (either one after the max mag segment, or the first one):
    if num_of_segments > 1:
        if seg_id == 0:
            # and multipy by the time interval to get the segment size in days
            seg2_size = (diff_bool_regions[seg_id+1].area/2)*ts_dTime
            # the division by 2 is due to the extra array row added above to trick 'measure.regionprops'
        else:
            # and multipy by the time interval to get the segment size in days
            seg2_size = (diff_bool_regions[0].area/2)*ts_dTime
    else:
        seg2_size = np.nan
    # ------------------ (k) --------------------
    # get parameters from the original TS (no running mean applied):
    # ---
    # exclude the values from the beginning and end of the TS according to the running mean TS
    # this is important because previous indices should math also the original TS.
    ts_pix_1_short = ts_pix_1.loc[ts_pix_1_mean.index.values]
    # get the difference from the original TS
    diff_org = ref_mean - ts_pix_1_short
    # get the values with the first segment
    diff_org_seg1 = diff_org.where(diff_bool_labels == diff_bool_regions[seg_id].label)
    # get the maximum magnitude
    max_mag_org_ind = np.nanargmax(np.abs(diff_org_seg1.values))
    max_mag_org = diff_org_seg1.values[max_mag_org_ind]
    max_mag_org_date = diff_org_seg1.index.values[max_mag_org_ind]
    # the time difference (in days) between original max magnitude and the running mean max magnitude
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


def features_from_yeodaTS(ts_pix_1):
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
    end_recovery = '2021-01-01'
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
    ts_pix_1_mean = ts_pix_1.rolling(sample_num, center=True).mean()
    # remove the NaN-s
    ts_pix_1_mean = ts_pix_1_mean.dropna()
    # remove everything after the recovery end period date
    ts_pix_1_mean = ts_pix_1_mean.loc[slice(start_history, end_recovery)]
    # ------------------ (b) --------------------
    # get the reference and its features:
    ref_mean = ts_pix_1_mean.loc[slice(start_history, end_history)].values.mean()
    # get the std og the running mean:
    ref_std = ts_pix_1_mean.loc[slice(start_history, end_history)].values.std()
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
    diff_bool[ts_diff.index.year.values == 2016] = False
    # ignore the values in the history period, i.e. the values in 2017 are set to False:
    diff_bool[ts_diff.index.year.values == 2017] = False
    # ------------------ (e) --------------------
    # check if there is at least one element above the margin:
    if np.all(diff_bool.values == False):
        exception_label = 2
        my_output[0] = exception_label
        my_output[2] = error_margin
        return my_output

    # ------------------ (f) --------------------
    # identify segments and get their properties
    # --
    # label segments in TS (consecutive labeled dates from the step "d")
    diff_bool_labels = measure.label(diff_bool.values)
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
        diff_bool_labels = measure.label(diff_bool_new.values)
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
        seg_id = seg_id[0]
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
    max_mag_date = ts_diff_seg1.index.values[max_mag_ind]
    # get the np array with relative times (in days) from the maximum disturbance:
    d_time = (ts_diff.index.values - max_mag_date).astype('timedelta64[D]').astype(int)
    # get the disturbance parameters related to the start and the end of the largest segment:
    t_pre = d_time[seg_ind_start]
    t_post = d_time[seg_ind_end]
    # !!! NOTE !!! this part is new. total time now considers multiple segments
    # in that case, the total time is between the start of the first segment and the end of the last segment
    if num_of_segments > 1:
        first_seg_start_ind = diff_bool_regions[0].bbox[1]
        last_seg_end_ind = diff_bool_regions[-1].bbox[3] - 1
        t_total = (ts_diff.index.values[last_seg_end_ind] -
                   ts_diff.index.values[first_seg_start_ind])\
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
    ts_pix_1_short = ts_pix_1.loc[ts_pix_1_mean.index.values]
    # get the difference from the original TS
    diff_org = ref_mean - ts_pix_1_short
    # get the values with the first segment
    diff_org_seg1 = diff_org.where(diff_bool_labels == diff_bool_regions[seg_id].label)
    # get the maximum magnitude
    max_mag_org_ind = np.nanargmax(np.abs(diff_org_seg1.values))
    max_mag_org = diff_org_seg1.values[max_mag_org_ind]
    max_mag_org_date = diff_org_seg1.index.values[max_mag_org_ind]
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



def plot_TS_fromPandaSeres(ts_pix_1):
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
    end_plot = '2021-01-01'
    # specify the sample number (3 months = 90 days / 6-day = 15):
    sample_num = 15
    # ------------------
    # calculate end plot time
    # ------------------
    ts_pix_1_mean = ts_pix_1.rolling(sample_num, center=True).mean()
    # plotting the TSs, the difference TS and the spatial location of the pixel
    #my_coor = [ts_pix_1.coords['x'].values.flatten()[0], ts_pix_1.coords['y'].values.flatten()[0]]
    # get the reference:
    ref_mean = ts_pix_1_mean.loc[slice(start_history, end_history)].dropna().values.mean()
    # get the std og the running mean:
    ref_std = ts_pix_1_mean.loc[slice(start_history, end_history)].dropna().values.std()
    # get the confidence intervals:
    [mean_DOWN, mean_UP] = stats.norm.interval(0.99, loc=ref_mean, scale=ref_std)
    total_num = ts_pix_1.count()
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
    ax_ts.fill_between(ts_pix_1.index.values, mean_DOWN*np.ones(total_num), mean_UP*np.ones(total_num),
                            facecolor='lightgreen', alpha=0.7)
    ax_ts.axhline(y=ref_mean, lw=2, ls='-', color='green')
    ax_ts.set_xticks(datelist)
    ax_ts.set_axisbelow(True)
    if ts_pix_1.name.split('_')[-1].upper() == 'VH':
        ax_ts.set_yticks(np.arange(-22, -9, 2))
    else:
        ax_ts.set_yticks(np.arange(-14, -3, 2))
    ax_ts.minorticks_on()
    ax_ts.grid(which='major', linestyle='-', linewidth=1)
    ax_ts.grid(which='minor', linestyle=':', linewidth=0.5)
    ax_ts.legend(['original TS', 'running mean TS', '99 % confidence', 'reference line'], loc="lower left")
    ax_ts.set_title(r'{} Pol. '.format(ts_pix_1.name.split('_')[-1].upper()))
    ax_ts.set_xlim((start_plot, end_plot))
    plt.tight_layout()
    
    

def plot_TS(ts_pix_1):
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
    end_plot = '2021-01-01'
    # specify the sample number (3 months = 90 days / 6-day = 15):
    sample_num = 15
    # ------------------
    # calculate end plot time
    # ------------------
    ts_pix_1_mean = ts_pix_1.rolling(sample_num, center=True).mean()
    # plotting the TSs, the difference TS and the spatial location of the pixel
    #my_coor = [ts_pix_1.coords['x'].values.flatten()[0], ts_pix_1.coords['y'].values.flatten()[0]]
    # get the reference:
    ref_mean = ts_pix_1_mean.loc[slice(start_history, end_history)].dropna().values.mean()
    # get the std og the running mean:
    ref_std = ts_pix_1_mean.loc[slice(start_history, end_history)].dropna().values.std()
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
    ax_ts.fill_between(ts_pix_1.index.values, mean_DOWN*np.ones(total_num), mean_UP*np.ones(total_num),
                            facecolor='lightgreen', alpha=0.7)
    ax_ts.axhline(y=ref_mean, lw=2, ls='-', color='green')
    ax_ts.set_xticks(datelist)
    ax_ts.set_axisbelow(True)
    if ts_pix_1.columns[0].upper() == 'VH':
        ax_ts.set_yticks(np.arange(-22, -9, 2))
    else:
        ax_ts.set_yticks(np.arange(-14, -3, 2))
    ax_ts.minorticks_on()
    ax_ts.grid(which='major', linestyle='-', linewidth=1)
    ax_ts.grid(which='minor', linestyle=':', linewidth=0.5)
    ax_ts.legend(['original TS', 'running mean TS', '99 % confidence', 'reference line'], loc="lower left")
    ax_ts.set_title(r'{} Pol. '.format(ts_pix_1.columns[0].upper()))
    ax_ts.set_xlim((start_plot, end_plot))
    plt.tight_layout()





