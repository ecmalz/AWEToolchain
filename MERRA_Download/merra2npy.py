"""
Accesses the downloaded files from NASA and creates a data npy. file.
1. >>merra_download.py<<  : Download wind Data
2. >>merra2npy.py<<       : read files and save as .npy in order to have the needed format for optimisation
3. >>clustering.py<<      : cluster .npy files to get different clusters
this code is executed by merra_download


Python Version 2.7 (not system version)
- Author: Elena Malz, Chalmers 2018, elenama@chalmers.se
"""

from netCDF4 import Dataset
import numpy as np
import sys
import math
import os
from scipy.stats import circmean #http://scipy.github.io/devdocs/generated/scipy.stats.circmean.html
import matplotlib.pylab as plt

# directory = '/Users/elenama/Documents/Git/Pamplona16/'
# prefix = '42.8_-1.6_20160101_20161231'

# start       = [2016, 01, 01]    # year, month, day
# end         = [2016, 12, 31]    # year, month, day
# levels      = np.arange(62, 73) # max 72




def visualize(directory, prefix):
    print '#### VISUALIZE ####'
    path = directory + prefix
    winds = np.load(path+'_winds.npy')
    features = np.load(path+'_features.npy')
    heights = np.load(path+'_heights.npy')
    size_w    =  winds.shape[0]
    size_h    = heights.shape[0]
    print("winds shape", winds.shape)
    print("heights shape", heights.shape)

    plt.figure()
    for level in np.arange(size_h):
        plt.plot(heights[level])
    plt.title( 'Pressure level heights ')
    plt.xlabel('time step (3-hour intervals)')
    plt.ylabel('pressure level height [m]')
    plt.grid('on')
    plt.savefig(directory+'heights.pdf')

    plt.figure()
    for level in np.arange(size_w):
        wind_magnitudes = np.linalg.norm(winds[level,:,:],axis=1)
        plt.plot(wind_magnitudes)
    plt.title( ' Wind absolute values')
    plt.xlabel('time step (3-hour intervals)')
    plt.ylabel('pressure level wind magnitude [m/s]')
    plt.grid('on')
    plt.savefig(directory+'winds.pdf')

    plt.figure()
    for level in np.arange(size_w-3,size_w):
        plt.plot(heights[level])
    plt.title('3 lowest heights')
    plt.xlabel('time step (3-hour intervals)')
    plt.ylabel('pressure level height [m]')
    plt.grid('on')
    plt.savefig(directory+'heights_2m_10m_50m.pdf')

    plt.figure()
    for level in np.arange(size_h-3,size_h):
        wind_magnitudes = np.linalg.norm(winds[level,:,:],axis=1)
        plt.plot(wind_magnitudes)
    plt.title('lowest heights for')
    plt.xlabel('time step (3-hour intervals)')
    plt.ylabel('pressure level wind magnitude [m/s]')
    plt.grid('on')
    plt.savefig(directory+'winds_2m_10m_50m.pdf')
    print '#### figures saved to directory ####'


def compute_angular_distance(reference, destination):
    if destination<reference:
        if reference-destination<math.pi:    # Go right from reference to dest
            return reference-destination
        else:
            return -(2*math.pi-reference+destination)
    else:
        if destination-reference<math.pi:    # Go left from reference to dest
            return -(destination-reference)
        else:
            return (reference+2*math.pi-destination)


def from_wind_to_features(eastwind, northwind):
    absolute = math.hypot(eastwind,northwind)
    degree_rad = np.arctan2(northwind,eastwind) # assume eastwind is x direction
    #degrees = np.rad2deg(degree_rad)
    #if degrees<0:
    #    degrees = degrees + 360
    if degree_rad<0:
        degree_rad = degree_rad + 2*math.pi
    return (absolute,degree_rad)

fileList_high = []
fileList_low = []
dayStrings = []
dict_timeseries_heights = dict()
dict_timeseries_wind = dict()
dict_timeseries_features = dict()


def _init_(directory, start, end, levels):
    print '#### Save merra to npy .... ####'
    years = [start[0]]
    months = range(start[1],end[1]+1)
    days = range(start[2],end[2]+1)
    # get all time stamps available
    for year in years:
        for month in months:
            for day in days:
                year_str = str(year)
                month_str = str(month).zfill(2)
                day_str   = str(day).zfill(2)
                suffix = year_str+month_str+day_str
                filename_high = directory+"MERRA2_400.tavg3_3d_asm_Nv."+suffix+".nc4"
                if (os.path.exists(filename_high)):
                    fileList_high.append(filename_high)
                    dayStrings.append(suffix)

                    # We also require to have the corresponding file for lowest
                    filename_low = directory+"MERRA2_400.tavg1_2d_slv_Nx."+suffix+".nc4"
                    if not os.path.exists(filename_low):
                        print("WARNING: Can not find corresponding lower level", filename_low)
                    else:
                        fileList_low.append(filename_low)

    ndays = len(fileList_high)
    print("Found", len(fileList_high), "matching files in date order for high altitudes")
    print("Parsing data into time series")
    print(dayStrings)

    for dayString in dayStrings:
        nc_f = directory+'MERRA2_400.tavg3_3d_asm_Nv.'+dayString+'.nc4'
        nc_f_low = directory+'MERRA2_400.tavg1_2d_slv_Nx.'+dayString+'.nc4'
        nc_fid = Dataset(nc_f, 'r')
        nc_fid_low = Dataset(nc_f_low, 'r')

    latlonPairs = []
    for lat in nc_fid.variables['lat'][:]:
        for lon in nc_fid.variables['lon'][:]:
            latlonPairs.append([lat,lon])

    latlonPairs_numpy = []
    for lat in range(len(nc_fid.variables['lat'][:])):
        for lon in range(len(nc_fid.variables['lon'][:])):
            latlonPairs_numpy.append([lat,lon])
    nlev = len(levels)

    # prepare the arrays filled with zeros
    for (index_lat,index_lon) in latlonPairs:
        timeseries_heights = np.zeros((nlev+3,8*ndays))
        timeseries_heights[nlev] = 50*np.ones((8*ndays)) # 50 metres (fixed)
        timeseries_heights[nlev+1] = 10*np.ones((8*ndays)) # 10 metres (fixed)
        timeseries_heights[nlev+2] = 2*np.ones((8*ndays))  # 2  metres (fixed)
        timeseries_wind = np.zeros((nlev+3,8*ndays), dtype='2float32')
        timeseries_features = np.zeros((nlev+3,8*ndays), dtype='2float32')

        dict_timeseries_heights[index_lat,index_lon] = timeseries_heights
        dict_timeseries_wind[index_lat,index_lon] = timeseries_wind
        dict_timeseries_features[index_lat,index_lon] = timeseries_features

    # Running through days and save the data in the dicts
    # Creating the features directly from the wind speed.
    day = 0
    for dayString in dayStrings:
        nc_f = directory+'MERRA2_400.tavg3_3d_asm_Nv.'+dayString+'.nc4'
        nc_f_low = directory+'MERRA2_400.tavg1_2d_slv_Nx.'+dayString+'.nc4'

        nc_fid = Dataset(nc_f, 'r')
        nc_fid_low = Dataset(nc_f_low, 'r')

        print(nc_f, "day", day + 1)

        for (index_lat,index_lon), (numpy_lat,numpy_lon) in  zip(latlonPairs, latlonPairs_numpy):
            for time in range(8):       # 8 3-hour intervals
                                # get the low levels
                timeindex = time*3
                level = len(levels) # At 50m: Call this "Level 73"
                eastwinds = nc_fid_low.variables['U50M'][timeindex:timeindex+3,numpy_lat,numpy_lon]
                northwinds = nc_fid_low.variables['V50M'][timeindex:timeindex+3,numpy_lat,numpy_lon]
                eastwind_mean = np.mean(eastwinds)
                northwind_mean = np.mean(northwinds)
                timeseries_wind[level,day*8+time] = [eastwind_mean, northwind_mean]
                features = from_wind_to_features(eastwind_mean,northwind_mean)
                timeseries_features[level,day*8+time] = [features[0],features[1]]

                level = len(levels)+1 # At 10m: Call this "Level 74"
                eastwinds = nc_fid_low.variables['U10M'][timeindex:timeindex+3,numpy_lat,numpy_lon]
                northwinds = nc_fid_low.variables['V10M'][timeindex:timeindex+3,numpy_lat,numpy_lon]
                eastwind_mean = np.mean(eastwinds)
                northwind_mean = np.mean(northwinds)
                timeseries_wind[level,day*8+time] = [eastwind_mean, northwind_mean]
                features = from_wind_to_features(eastwind_mean,northwind_mean)
                timeseries_features[level,day*8+time] = [features[0],features[1]]

                level = len(levels)+2 # At 2m: Call this "Level 75"
                eastwinds = nc_fid_low.variables['U2M'][timeindex:timeindex+3,numpy_lat,numpy_lon]
                northwinds = nc_fid_low.variables['V2M'][timeindex:timeindex+3,numpy_lat,numpy_lon]
                eastwind_mean = np.mean(eastwinds)
                northwind_mean = np.mean(northwinds)
                timeseries_wind[level,day*8+time] = [eastwind_mean, northwind_mean]
                features = from_wind_to_features(eastwind_mean,northwind_mean)
                timeseries_features[level,day*8+time] = [features[0],features[1]]

                for level in range(nlev):     # 72 pressure levels
                    timeseries_heights[level,day*8+time] =  nc_fid.variables['H'][time,level,numpy_lat,numpy_lon] -  nc_fid.variables['PHIS'][time,numpy_lat,numpy_lon]/9.81 # height above sealevel - height of surface/9.81 = actual height (source: NASA email)
                    eastwind = nc_fid.variables['U'][time,level,numpy_lat,numpy_lon]                           # eastwind
                    northwind = nc_fid.variables['V'][time,level,numpy_lat,numpy_lon]                          # northwind
                    timeseries_wind[level,day*8+time] = [eastwind,northwind]
                    features = from_wind_to_features(eastwind,northwind)
                    timeseries_features[level,day*8+time] = [features[0],features[1]]

                circular_mean = circmean(timeseries_features[:,day*8+time,1])
                for level in range(len(timeseries_features)):
                    timeseries_features[level,day*8+time,1] = compute_angular_distance(circular_mean,timeseries_features[level,day*8+time,1])

            dict_timeseries_heights[index_lat,index_lon][:,day*8:day*8+8] = timeseries_heights[:,day*8:day*8+8]
            dict_timeseries_wind[index_lat,index_lon][:,day*8:day*8+8] = timeseries_wind[:,day*8:day*8+8]
            dict_timeseries_features[index_lat,index_lon][:,day*8:day*8+8]= timeseries_features[:,day*8:day*8+8]

        day = day + 1

    # Save to npy:
    start = dayStrings[0]
    stop = dayStrings[-1]



    for (index_lat,index_lon) in latlonPairs:
        prefix = str(format(index_lat,'.2f'))+'_'+str(format(index_lon,'.2f'))+'_'+start+'_'+stop


        np.save(directory+prefix+'_winds',dict_timeseries_wind[index_lat,index_lon][4:,:,:])
        np.save(directory+prefix+'_heights',dict_timeseries_heights[index_lat,index_lon][4:,:])
        np.save(directory+prefix+'_features',dict_timeseries_features[index_lat,index_lon][4:,:,:])
    return prefix





# # --- init Save to npy
# prefix = _init_(directory,start, end, levels)
#
# # --- Plot the result
# # prefix = '57.5_11.875_20160101_20161231'
# visualize(directory, prefix)

# -------------------------

# # --- CHECK GOTHENBURG .. something doesn't fit
# loc = '/Users/elenama/Documents/Git/Data/'
# winds = np.load(loc+'gothenburg_20160101_20161231_wind_lowest10.npy')
# features = np.load(loc+'gothenburg_20160101_20161231_features_lowest10.npy')
# mywinds = np.load('/Users/elenama/Documents/Git/Goteborg16/57.5_11.875_20160101_20160103_winds.npy')
# myfeatures = np.load('/Users/elenama/Documents/Git/Goteborg16/57.5_11.875_20160101_20160103_features.npy')
# print winds[:,1]
# print mywinds[:,1]
# print features[:,1]
# print myfeatures[:,1]

#
# # --- CHECK VIMMERBY .. hm features could be wrong since i generated them
# myheights = np.load('/Users/elenama/Documents/Git/Vimmerby12/57.5_15.625_20120101_20120103_heights.npy')
# mywinds = np.load('/Users/elenama/Documents/Git/Vimmerby12/57.5_15.625_20120101_20120103_winds.npy')
# myfeatures = np.load('/Users/elenama/Documents/Git/Vimmerby12/57.5_15.625_20120101_20120103_features.npy')
#
# loc = '/Users/elenama/Documents/Git/Data/2016GotVimmerby/'
# winds = np.load(loc+'ts_vimmerby_20120101_20121231_wind_lowest15.npy')
# heights = np.load(loc+'ts_vimmerby_20120101_20121231_heights_lowest15.npy')
# features = np.load(loc+'ts_vimmerby_20120101_20121231_features_lowest15.npy')
# print(winds[:,1])
# print mywinds[1:,1]
# print features[:,1]
# print myfeatures[1:,1]
