"""
Automated Download MERRA wind data. Generates a foler with desired name and downloads for the set time frame and coordinates the NW winddata for high and low altitudes.
Setup as : "https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+cURL+And+Wget"
High altitudes: https://disc.sci.gsfc.nasa.gov/datasets/M2T3NVASM_V5.12.4/summary
low altitudes: https://disc.sci.gsfc.nasa.gov/datasets/M2T1NXSLV_V5.12.4/summary?%2F
1. >>merra_download.py<<  : Download wind Data
2. >>merra2npy.py<<       : read files and save as .npy in order to have the needed format for optimisation
3. >>clustering.py<<      : cluster .npy files to get different clusters
(steps 2 and 3 can be called directly from this file!)
Python Version 2.7 (not system version)
- Author: Elena Malz, Chalmers 2017, elenama@chalmers.se


"""

from _downloadwind import *
import math
from scipy.stats import circmean #http://scipy.github.io/devdocs/generated/scipy.stats.circmean.html


cities = {'pamplona': [42.812526,-1.645775], 'hamburg': [53.551085,9.993682], 'esbjerg': [55.476466,8.459405],
         'stockholm': [59.329323,18.068581], 'gothenburg': [57.708870,11.974560], 'paris': [48.856614,2.352222],
          'vimmerby':[57.665556,15.854722], 'redsea':[21.732021,36.625196], 'zermatt':[46, 7.75], 'laPalma':[28.75,-17.99]}

# --- specify coordinates, dates, variables needed
city = 'pamplona'
coordinates = [cities[city][1],cities[city][1], cities[city][0], cities[city][0] ]  # Content: [minlongitude, maxlongitude, minlatitude , maxlatitude]
# coordinates = [ -180, 180,-90,90]  # Content: [minlongitude, maxlongitude, minlatitude , maxlatitude]

directory = '/Users/elenama/Documents/Git/' + city+'16/'  # -- Mkdir, Chdir and create text files in that directory where we also want to download the stuff


# coordinates = [11, 12, 57, 58]  # Content: [minlongitude, maxlongitude, minlatitude , maxlatitude]

start       = [2016, 01, 01]    # year, month, day
end         = [2016, 12, 31]    # year, month, day
# levels      = np.array([63, 70, 72]) # max 72 np.arange(62, 73)
levels      = np.arange(62, 73)

variables   = ['u', 'v', 'h', 'phis' ]
# u = eastwind, v = west wind, pl = pressure, h = height above sealevel, phis = height of surface/9.81 ...

url = create_MerraURL(coordinates, start, end, levels, variables)               # Get MERRA URL from coordinates
url_lowheights = create_lowheights_MerraURL(coordinates, start, end, variables)

if not os.path.exists(directory):
    os.makedirs(directory)
os.chdir(directory)                               # Change to that directory
name = 'wind'                                     # Give the textfile a name
curlfile = create_txtfile(url,directory,name) # Create textfile and get its name
curlfile_lowheights = create_txtfile(url_lowheights,directory,'lowheights')
print '###############################'
print '####   START DOWNLOAD    ##### '
print '###############################'
DoDownload(curlfile)
DoDownload(curlfile_lowheights)


# -- just check:
import glob
os.chdir(directory)
try:
    nc_f = glob.glob('*.nc4')[0]
    nc_fid = Dataset(nc_f, 'r')
    print nc_fid.variables['lat'][:]
    print 'download successful'
except: print '\n\n File cannot be opened, probably encrypted nc4 due to failed MERRA download'

for i in glob.glob('*.nc4'):
    nc_fid = Dataset(i, 'r')
    print nc_fid.Comment


import merra2npy
import clustering
#  init Save to npy

prefix = merra2npy._init_(directory,start, end, levels)
# Plot the result
merra2npy.visualize(directory, prefix)
#  cluster data
clustering.cluster(directory,prefix)
# Save cluster wind speed pairs
clustering.getClusterwindpairs(directory,prefix)
