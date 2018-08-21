
"""
Automated Download MERRA wind data
Setup as : "https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+cURL+And+Wget"
Python Version 2.7 (not system version)
- Author: Elena Malz, Chalmers 2017, elenama@chalmers.se
"""

# run with with anaconda python
from cookielib import CookieJar
from urllib import urlencode
import urllib2
from netCDF4 import Dataset
import requests
from BeautifulSoup import BeautifulSoup
import numpy as np
import sys
import os



def create_MerraURL(coordinates, start, end, levels, variables):
    """ Data from : https://disc.sci.gsfc.nasa.gov/datasets/M2T3NVASM_V5.12.4/summary """
    # --- generate the URL for downloading textfiles including the data url's
    textlevels = str()
    for l in levels:
        textlevels = textlevels+'&levels='+str(l)
    textvar = str()
    for v in variables:
        textvar = textvar + '&variable='+v

    url = ('https://goldsmr5.gesdisc.eosdis.nasa.gov/cgi-bin/OTF/HTTP_DpFileDownloadMERRA2.pl?DATASET=MERRA_DP&FCP_DIR=/ftp/private/tmp/&APPLICATION=SUBSET_MERRA2&FILTER=SUBSET_MERRA2&SUB_'
    'LONMIN='+str(coordinates[0])+'&SUB_LONMAX='+str(coordinates[1])+'&SUB_LATMAX='+str(coordinates[3])+'&SUB_LATMIN='+str(coordinates[2])+'&OUTPUT_FORMAT=nc4&LOOKUPID_List=M2T3NVASM&'
    'STARTYR='+str(start[0]).zfill(2)+'&STARTMON='+str(start[1]).zfill(2)+'&STARTDAY='+str(start[2]).zfill(2)+'&ENDYR='+str(end[0]).zfill(2)+'&ENDMON='+str(end[1]).zfill(2)+'&ENDDAY='+str(end[2]).zfill(2)+'&PREGRID=NONE'
    +textlevels+'&'+textvar+'&')
    return url

def create_lowheights_MerraURL(coordinates, start, end, variables):
    """ Data from : https://disc.sci.gsfc.nasa.gov/datasets/M2T1NXSLV_V5.12.4/summary?%2F """
    # --- generate the URL for downloading textfiles including the data url's
    variables = ['u10m', 'u2m', 'u50m', 'v10m', 'v2m', 'v50m']
    textvar = str()
    for v in variables:
        textvar = textvar + '&variable='+v
    url =('https://goldsmr4.gesdisc.eosdis.nasa.gov/cgi-bin/OTF/HTTP_DpFileDownloadMERRA2.pl?DATASET=MERRA_DP&FCP_DIR=/ftp/private/tmp/&APPLICATION=SUBSET_MERRA2&FILTER=SUBSET_MERRA2&SUB_'
    'LONMIN='+str(coordinates[0])+'&SUB_LONMAX='+str(coordinates[1])+'&SUB_LATMAX='+str(coordinates[3])+'&SUB_LATMIN='+str(coordinates[2])+'&OUTPUT_FORMAT=nc4&LOOKUPID_List=M2T1NXSLV&'
    'STARTYR='+str(start[0]).zfill(2)+'&STARTMON='+str(start[1]).zfill(2)+'&STARTDAY='+str(start[2]).zfill(2)+'&ENDYR='+str(end[0]).zfill(2)+'&ENDMON='+str(end[1]).zfill(2)+'&ENDDAY='+str(end[2]).zfill(2)
    +textvar+'&')
    return url

def create_textfileURL(page,):
    start_link = page.find("Save the")
    if start_link == -1:
        return None, 0
    startURL =  page.find("a href")
    start_quote = page.find('"', start_link)
    end_quote = page.find('"', start_quote + 1)
    urlend = page[start_quote + 1: end_quote]
    if page.find("goldsmr4") != -1:
        main =  'https://goldsmr4.gesdisc.eosdis.nasa.gov/cgi-bin/OTF/'
    else:
        main =  'https://goldsmr5.gesdisc.eosdis.nasa.gov/cgi-bin/OTF/'
    url = main + urlend
    return url

def create_txtfile(url, directory, name):
#     """open NASA merra data URL and create a text file listing the download URL's"""
    response    = requests.get(url)
    page        = str(BeautifulSoup(response.content)) # parse html
    txtfileURL  = create_textfileURL(page)     # Get the link for downloading the textfiles containing the URL's pointing at the data files
    txtfile     = requests.get(txtfileURL)  # Download the text file
    file = open(directory+name+'.txt', 'w')
    file.write(txtfile.text)
    file.close()

    # Reopen the text files and add an output to the download URL
    f = open(directory+name+'.txt', 'r')
    f2 = open(directory+name+'_curl.txt','w')
    for l in f:
        #print 'download nc4 file url : ', l
        fname = "MERRA2"+l.split('nc4')[0].split('MERRA2')[2]+"nc4"
        f2.write('url = ' + l)
        f2.write('output = '+ fname+'\n')
    f.close()
    f2.close()
    return name+'_curl.txt'


def DoDownload(downloadfile):
    # mystring = "curl -L -b ~/.urs_cookies -c ~/.urs_cookies -K ./" + downloadfile + ' --ipv4' #retry 3 --retry-delay 3 --retry-max-time 3
    mystring = "curl -n -c ~/.urs_cookies -b ~/.urs_cookies  -LJ  -K ./" + downloadfile  + ' --ipv4' #retry 3 --retry-delay 3 --retry-max-time 3
    print(mystring)
    os.system(mystring)
