#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 16:52:39 2022

@author: vsassard
"""

import sys
import obspy
import os
import time
import numpy as np
import pandas as pd
from mpi4py import MPI
from seisgo.utils import split_datetimestr
from seisgo import downloaders,noise
#########################################################
################ PARAMETER SECTION ######################
#########################################################

################# DOWNLOAD 4 CHANNELS FOR THE OBS #########################
tt0=time.time()

# paths and filenames
rootpath = "AACSE_region" # roothpath for the project
direc  = os.path.join(rootpath,'Raw')                   # where to store the downloaded data
#if not os.path.isdir(direc): os.mkdir(direc)
down_list  = os.path.join(direc,'station_off.txt')
# CSV file for station location info

# download parameters
source='IRIS'                                 # client/data center. see https://docs.obspy.org/packages/obspy.clients.fdsn.html for a list
max_tries = 3                                                  #maximum number of tries when downloading, in case the server returns errors.
use_down_list = False                                                # download stations from a pre-compiled list or not
flag      = True                                               # print progress when running the script; recommend to use it at the beginning
samp_freq = 5                                                  # targeted sampling rate at X samples per seconds
rmresp   = True                                                # instrument response removal
rmresp_out = 'VEL'
pressure_chan = ['HDH']				#Added by Xiaotao Yang. This is needed when downloading some special channels, e.g., pressure data. VEL output for these channels.
respdir   = os.path.join(rootpath,'resp')                       # directory where resp files are located (required if rm_resp is neither 'no' nor 'inv')
freqmin   = 0.002                                                # pre filtering frequency bandwidth
freqmax   = 0.5*samp_freq
# note this cannot exceed Nquist freq

# targeted region/station information: only needed when use_down_list is False
lamin,lamax,lomin,lomax= 50.792,60.011,-165.2052,-146.9971                # regional box: min lat, max lat, min lon, max lon
chan_list = ["HH1","HH2","HHZ","BHN","BHE","BHZ","HDH"]
net_list  = ["XO"]                   # network list
sta_list  = ["*"]                                               # station (using a station list is way either compared to specifying stations one by one)
start_date = "2018_01_01_0_0_0"                               # start date of download
end_date   = "2020_07_01_0_0_0"                               # end date of download
inc_hours  = 12                                                 # length of data for each request (in hour)
maxseischan = 3                                                  # the maximum number of seismic channels, excluding pressure channels for OBS stations.
ncomp      = maxseischan #len(chan_list)

# get rough estimate of memory needs to ensure it now below up in noise cross-correlations
cc_len    = 2*3600                                                # basic unit of data length for fft (s)
step      = 1*3600                                                 # overlapping between each cc_len (s)
MAX_MEM   = 5.0                                                 # maximum memory allowed per core in GB


##################################################
# we expect no parameters need to be changed below
# assemble parameters used for pre-processing waveforms in downloading
prepro_para = {'rmresp':rmresp,'rmresp_out':rmresp_out,'respdir':respdir,'freqmin':freqmin,'freqmax':freqmax,\
                'samp_freq':samp_freq}

downlist_kwargs = {"source":source, 'net_list':net_list, "sta_list":sta_list, "chan_list":chan_list, \
                    "starttime":start_date, "endtime":end_date, "maxseischan":maxseischan, "lamin":lamin, "lamax":lamax, \
                    "lomin":lomin, "lomax":lomax, "pressure_chan":pressure_chan, "fname":down_list}

########################################################
#################DOWNLOAD SECTION#######################
########################################################
#--------MPI---------
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if rank==0:
    if flag:
        print('station.list selected [%s] for data from %s to %s with %sh interval'%(use_down_list,start_date,end_date,inc_hours))

    if not os.path.isdir(direc):os.makedirs(direc)
    if use_down_list:
        stalist=pd.read_csv(down_list)
    else:
        stalist=downloaders.get_sta_list(**downlist_kwargs) # saves station list to "down_list" file, PANDA DATAFRAME
        # save parameters for future reference
        metadata = os.path.join(direc,'download_off_info.txt')
        fout = open(metadata,'w')
        fout.write(str({**prepro_para,**downlist_kwargs,'inc_hours':inc_hours,'ncomp':ncomp}));fout.close()
        
        
        
        
        
        
################# DOWNLOAD VERTICAL CHANNEL ONLY FOR THE ONSHORE STATIONS ##################
tt0=time.time()

# paths and filenames
rootpath = "AACSE_region" # roothpath for the project
direc  = os.path.join(rootpath,'Raw')                   # where to store the downloaded data
#if not os.path.isdir(direc): os.mkdir(direc)
down_list  = os.path.join(direc,'station_on.txt')
# CSV file for station location info

# download parameters
source='IRIS'                                 # client/data center. see https://docs.obspy.org/packages/obspy.clients.fdsn.html for a list
max_tries = 3                                                  #maximum number of tries when downloading, in case the server returns errors.
use_down_list = False                                                # download stations from a pre-compiled list or not
flag      = True                                               # print progress when running the script; recommend to use it at the beginning
samp_freq = 5                                                  # targeted sampling rate at X samples per seconds
rmresp   = True                                                # instrument response removal
rmresp_out = 'VEL'
pressure_chan = ['HDH']				#Added by Xiaotao Yang. This is needed when downloading some special channels, e.g., pressure data. VEL output for these channels.
respdir   = os.path.join(rootpath,'resp')                       # directory where resp files are located (required if rm_resp is neither 'no' nor 'inv')
freqmin   = 0.002                                                # pre filtering frequency bandwidth
freqmax   = 0.5*samp_freq
# note this cannot exceed Nquist freq

# targeted region/station information: only needed when use_down_list is False
lamin,lamax,lomin,lomax= 50.792,60.011,-165.2052,-146.9971                # regional box: min lat, max lat, min lon, max lon
chan_list = ["HHZ","BHZ"]
net_list  = ["6J","AK","AT","AV","GM","II","TA"]                   # network list
sta_list  = ["*"]                                               # station (using a station list is way either compared to specifying stations one by one)
start_date = "2018_01_01_0_0_0"                               # start date of download
end_date   = "2020_07_01_0_0_0"                               # end date of download
inc_hours  = 12                                                 # length of data for each request (in hour)
maxseischan = 1                                                  # the maximum number of seismic channels, excluding pressure channels for OBS stations.
ncomp      = maxseischan #len(chan_list)

# get rough estimate of memory needs to ensure it now below up in noise cross-correlations
cc_len    = 2*3600                                                # basic unit of data length for fft (s)
step      = 1*3600                                                 # overlapping between each cc_len (s)
MAX_MEM   = 5.0                                                 # maximum memory allowed per core in GB


##################################################
# we expect no parameters need to be changed below
# assemble parameters used for pre-processing waveforms in downloading
prepro_para = {'rmresp':rmresp,'rmresp_out':rmresp_out,'respdir':respdir,'freqmin':freqmin,'freqmax':freqmax,\
                'samp_freq':samp_freq}

downlist_kwargs = {"source":source, 'net_list':net_list, "sta_list":sta_list, "chan_list":chan_list, \
                    "starttime":start_date, "endtime":end_date, "maxseischan":maxseischan, "lamin":lamin, "lamax":lamax, \
                    "lomin":lomin, "lomax":lomax, "pressure_chan":pressure_chan, "fname":down_list}

########################################################
#################DOWNLOAD SECTION#######################
########################################################
#--------MPI---------
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if rank==0:
    if flag:
        print('station.list selected [%s] for data from %s to %s with %sh interval'%(use_down_list,start_date,end_date,inc_hours))

    if not os.path.isdir(direc):os.makedirs(direc)
    if use_down_list:
        stalist=pd.read_csv(down_list)
    else:
        stalist=downloaders.get_sta_list(**downlist_kwargs) # saves station list to "down_list" file, PANDA DATAFRAME
            # save parameters for future reference
        metadata = os.path.join(direc,'download_on_info.txt')
        fout = open(metadata,'w')
        fout.write(str({**prepro_para,**downlist_kwargs,'inc_hours':inc_hours,'ncomp':ncomp}));fout.close()
        
        
        
        
############ MERGE FILES CREATED #####################
