"""
Set of utilities written by Claire Bedbrook for counting peaks in jellyfish
activity trace
"""

import matplotlib.pyplot as plt
import skimage.io
import os
import numpy as np
import scipy.signal
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import jb_utils as jb
import cb_utils as cb

###############################################################################
###############################################################################
   
   
def peak_counter(df, jelly, thresh, fps, lam_big, lam_small, data_start, data_end, n):
    """
    Input:
        df: dataframe with each column for one jellyfish and each row for a 
                single frame's measured pixel intentisty 
        jelly: jellyfish number
        thresh: threshold for normalized activtiy trace to i.d. peaks
        fps: frames per second
        lam_big: filtering parameter for basline
        lam_small: filtering parameter for trace smoothing
        data_start: first frame (time = 0)
        data_end: last frame of jump (time = 20 min)
        n: jump number
    
    Returns:
        peak_dist_lst: distance (# frames) between peaks
        peak_dist_t_lst: distance (time) between peaks, i.e. IPI
        peak_numb_lst: # number of peaks within time period (20 min)
        peak_inds: array of frame number with pulse peak
    """
    
    # Arrays for relevant peak counting values
    peak_dist_lst = [] # distance (# frames) between peaks
    peak_dist_t_lst = [] # distance (time) between peaks, i.e. IPI
    peak_numb_lst = [] # number of peaks
    
    # check if the file has data for analysis                  
    if np.isnan(df[jelly][1]) == True:
        peak_dist_lst.append(0)
        peak_numb_lst.append(0)
        print('no data, cannot count peaks')
    
    # if data, then continue   
    else:                
        # Scale by activity trace by 1000 for filtering purposes
        trace = 1000 * df[jelly][data_start:data_end].values
        
        # Array of frames in the file
        frames = np.array(range(len(df[jelly][data_start:data_end])))
        
        # Convert frames to time in seconds 
        time = frames * 1.0 / (fps) # time in s
        
        # Find the baseline by smoothing the trace  
        baseline = jb.nw_kernel_smooth(frames, frames, trace, jb.epan_kernel, lam_big)
        
        # Find the max based on the delta from baseline to max hight of the pulse
        # from a section of the trace
        delta_max1 = np.max(trace[0:1000] - baseline[0:1000])
        delta_max2 = np.max(trace[2000:3000] - baseline[2000:3000])
        delta_max3 = np.max(trace[4000:5000] - baseline[4000:5000])
        delta_max4 = np.max(trace[6000:7000] - baseline[6000:7000])
        delta_max5 = np.max(trace[9000:10000] - baseline[9000:10000])
        delta_max6 = np.max(trace[11000:12000] - baseline[11000:12000])
        delta_max7 = np.max(trace[14000:15000] - baseline[14000:15000])
        delta_max8 = np.max(trace[16000:17000] - baseline[16000:17000])
        delta_max9 = np.max(trace - baseline)
        
        # average the max measurements 
        delta_max = (delta_max1 + delta_max2 + delta_max3 + delta_max4 \
                    + delta_max5 + delta_max6 + delta_max7 + delta_max8 \
                    + delta_max9)/9
                
        # Use the delta between baseline and top of the pulse to normalize the 
        # trace:
        trace_norm = (trace - baseline) / ((baseline + delta_max) - baseline)
        
        # Find max and min of the trace to find peaks 
        max_inds = scipy.signal.argrelmax(trace_norm, order=5)    
                
        # Find places where two contiguous points are equal- i.e. looking
        # for the flat peaks
        plataeu_peak= []
        for i in range(len(trace)-1):
            if (trace[i+1] - trace[i] == 0.0) & (trace_norm[i] > thresh):
                plataeu_peak.append(i)
                max_inds = np.column_stack((max_inds, i))
        
        # sort full list of max_inds with appended plataeu_peak 's
        max_inds = np.sort(max_inds[0])
                                
        # Find peak_inds: go through list of all max_inds and exclude max 
        # within four frames of each other
        peak_inds = []
        for i in range(len(max_inds) - 1):
            if (trace_norm[max_inds][i] > thresh) & \
                    (max_inds[i+1] != max_inds[i]) & \
                    (max_inds[i+1] != max_inds[i] + 1) & \
                    (max_inds[i+1] != max_inds[i] + 2) & \
                    (max_inds[i+1] != max_inds[i] + 3) & \
                    (max_inds[i+1] != max_inds[i] + 4):
                peak_inds.append(max_inds[i])
        
        # Plot the normalized trace with max min
        plt.figure('Hour_' + str(n))
        plt.plot(time, trace_norm, alpha=0.5)
        plt.plot(time[peak_inds], trace_norm[peak_inds], '.r')
        plt.xlabel('Time [s]', fontsize=16)
        plt.ylabel('Normalized Intensity', fontsize=16)
        
        # Compute the total number of peaks
        peak_numb_lst.append(len(peak_inds))  
        
        # Compute distance between peaks (i.e. IPI)
        for i in range(len(peak_inds)):
            if i < len(peak_inds)-1:
                peak_dists = peak_inds[i + 1] -  peak_inds[i]
                peak_dists_t = time[peak_inds][i]
                peak_dist_lst.append(peak_dists / (15.0))  # for 15 fps
                peak_dist_t_lst.append(peak_dists_t)
    plt.show()
    return [peak_dist_lst, peak_dist_t_lst, peak_numb_lst, peak_inds] 


def df_maker_IPI_time_stamp(p_lst, n):
    """
    make dataframe for single jellyfish of inter-peak distance with time-stamp
    (i.e. IPI times)
    Input:
        p_lst: array of frame number with pulse peak for each jump (hour)
        n: item, jump (hour) 
    Output:
        df: dataframe for each jellyfish for each hour with time of each pulse
            and inter-pulse interval (i.e. IPI times)
    """
    columns = len(p_lst)
    i_test = []
    
    for i in range(columns):
        i_test.append(len(p_lst[i]))
    
    # Find the jump with the most IPI (rows) this will be the length of the df
    ind = np.max(i_test)
    
    # Add NaN to each column until length of all columns are equale
    for col in range(columns):
        if len(p_lst[col]) < ind:
            for w in range(ind - len(p_lst[col])):
                p_lst[col].append('NaN')
                
    df = pd.DataFrame(index=range(ind))
    
    n = n + 1
    df['time_J' + str(n-1)] = p_lst[2 * n - 2]
    df['IPI_J' + str(n-1)] = p_lst[2 * n-1]   
    return df

def df_maker_peak_time(p_lst, n):
    """
    make dataframe for each jellyfish/each hour with the frame number of each pulse
    
    Input:
        p_lst: array of frame number with pulse peak for each jump (hour)
        n: item, jump (hour) 
    Output:
        df: dataframe for single jellyfish of for single jellyfish of IPI times
            each column is jump (hour)
    """
    # Prep df for peak time 
    columns = len(p_lst)
    i_test = []
    for i in range(columns):
        i_test.append(len(p_lst[i]))
    
    ind = np.max(i_test)
    
    for col in range(columns):
        if len(p_lst[col]) < ind:
            for w in range(ind - len(p_lst[col])):
                p_lst[col].append('NaN')
                
    df = pd.DataFrame(index=range(ind))
       
    df['Peak_frame_J' + str(n)] = p_lst[n]
        
    return df  