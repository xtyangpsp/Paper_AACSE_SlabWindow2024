#!/usr/bin/env python
# coding: utf-8

# # Tilt and Compliance Corrections for OBS Data: Continuous
# ### Xiaotao Yang @ Harvard University
# This notebook contains examples of compliance corrections using local data on the disk. The functions for tilt and compliance corrections are in module seispy.obsmaster.

# ## Step 0. Load needed packages.
# Some functions are imported from the utils.py and the obsmaster.py.

# In[ ]:


#import needed packages.
from seisgo import utils
from seisgo import obsmaster as obs
import sys
import time
import scipy
import obspy
import pyasdf
import datetime
from mpi4py import MPI
import os, glob
import numpy as np
# import pandas as pd
# import matplotlib.pyplot  as plt
# from obspy import UTCDateTime
from obspy.core import Stream, Trace

t0=time.time()
"""
1. Set global data path parameters.
"""
rootpath='/depot/xtyang/data/projects/vsassard'
rawdatadir = os.path.join(rootpath,'AACSE_region/Raw')
downloadexample=False #change to False or remove/comment this block if needed.

#directory to save the data after TC removal
tcdatadir = os.path.join(rootpath,'AACSE_region/TC_corrected')
cleantcdatadir=True #If True, the program will remove all *.h5 files under `tcdatadir` before running.

#parameters for orientation corrections for horizontal components
correct_obs_orient=True
obs_orient_file='orientation_AACSE_FINAL.csv'

#how to deal with stations with bad traces
drop_if_has_badtrace=True

"""
2. Tilt and compliance removal parameters
"""
window=7200
overlap=0.2
taper=0.08
qc_freq=[0.004,0.2]
plot_correction=False
normalizecorrectionplot=False
tc_subset=['ZP-H']
savetcpara=True                     #If True, the parameters are saved to a text file
                                    #in the `tcdatadir` directory.
requirePressure=True
for tcs in tc_subset:
    if 'P' in tcs: requirePressure=True; break
tcparaoutfile=os.path.join(tcdatadir,'tcparameters.txt')
OBS_cor_out = os.path.join(tcdatadir,'OBS_cor.txt')
OBS_uncor_out = os.path.join(tcdatadir,'OBS_uncor.txt')
land_out = os.path.join(tcdatadir,'land.txt')
#assemble all parameters into a dictionary.
tcpara={'window':window,'overlap':overlap,'taper':taper,'qc_freq':qc_freq,
        'tc_subset':tc_subset}

######################################################################
#### Normally, no changes are needed for the following processing ####
######################################################################
"""
3. Read local data and do correction
We use the wrapper function for tilt and compliance corrections. The data after noise removal
will be saved to the original file name BUT in different directory defined by `tcdatadir`.
"""
#-------- Set MPI parameters --------------------------------
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
if rank==0:
    ####################################
    #### Optional clean-up block ####
    ####################################
    if not os.path.isdir(rawdatadir):
        comm.barrier()
        raise IOError('Abort! Directory for raw data NOT found: '+rawdatadir)
        sys.exit()

    if not os.path.isdir(tcdatadir): os.mkdir(tcdatadir)
    dfilesTC0 = glob.glob(os.path.join(tcdatadir,'*.h5'))
    if cleantcdatadir and len(dfilesTC0)>0:
        print('Cleaning up TC removal directory before running ...')
        for df0 in dfilesTC0:os.remove(df0)
    ####################################
    ##### End of clean-up block #####
    ####################################
    print(tcpara)
    if savetcpara:
        fout = open(tcparaoutfile,'w')
        fout.write(str(tcpara));
        fout.close()
    dfiles0 = glob.glob(os.path.join(rawdatadir,'*.h5'))
    nfiles = len(dfiles0)
    splits0  = nfiles
    if nfiles < 1:
        raise IOError('Abort! no available seismic files in '+rawdatadir)
        sys.exit()
else:
    splits0,dfiles0 = [None for _ in range(2)]

# broadcast the variables
splits = comm.bcast(splits0,root=0)
dfiles  = comm.bcast(dfiles0,root=0)
#--------End of setting MPI parameters -----------------------
#make lists of corrected and uncorrected stations
sta_cor = []
sta_uncor = []
OBS_uncor = []

for ifile in range(rank,splits,size):
    #read obs orientation data.
    if correct_obs_orient:
        try:
            obs_orient_data=obs.get_orientations(obs_orient_file)
        except Exception as e:
            print(e)
            sys.exit()
    df=dfiles[ifile]
    print('Working on: '+df+' ... ['+str(ifile+1)+'/'+str(len(dfiles))+']')
    dfbase=os.path.split(df)[-1]
    df_tc=os.path.join(tcdatadir,dfbase)

    ds=pyasdf.ASDFDataSet(df,mpi=False,mode='r')
    netstalist = ds.waveforms.list()
    nsta = len(netstalist)

    tilt=[]
    sta_processed=[]
    for ista in netstalist:
        print('  station: '+ista)
        """
        Get the four-component data
        """
        try:
            inv = ds.waveforms[ista]['StationXML']
        except Exception as e:
            print('  No stationxml for %s in file %s'%(ista,df))
            inv = None

        all_tags = ds.waveforms[ista].get_waveform_tags()
        if len(all_tags) < 1 or len(all_tags) > 4: continue #empty waveform group.
        else: print(all_tags)

        tr1, tr2, trZ, trP=[None for _ in range(4)]
        #assign components by waveform tags.
        #This step may fail if the tags don't reflect the real channel information
        newtags=['-','-','-','-']
        for tg in all_tags:
            tr_temp = ds.waveforms[ista][tg][0]
            chan = tr_temp.stats.channel
            if chan[-1].lower() == 'h':trP=tr_temp;newtags[3]=tg
            elif chan[-1].lower() == '1' or chan[-1].lower() == 'e':tr1=tr_temp;newtags[0]=tg
            elif chan[-1].lower() == '2' or chan[-1].lower() == 'n':tr2=tr_temp;newtags[1]=tg
            elif chan[-1].lower() == 'z':trZ=tr_temp;newtags[2]=tg

        #sanity check.
        badtrace=False
        hasPressure=False
        if isinstance(trP,Trace):
            hasPressure=True
        if isinstance(trZ, Trace) and not isinstance(tr1, Trace) and not isinstance(tr2, Trace) and not isinstance(trP, Trace):
            print("Only vertical component, land station")
            newtags_tmp =[]
            sta_uncor.append(ista)
            outstream = Stream(traces = [trZ])
            newtags_tmp.append(utils.get_tracetag(trZ))
            utils.save2asdf(df_tc, outstream, newtags_tmp, sta_inv=inv)
        if isinstance(trZ, Trace) and isinstance(tr2, Trace) and not isinstance(tr1, Trace) or isinstance(trZ, Trace) and isinstance(tr1, Trace) and not isinstance(tr2, Trace) or isinstance(trZ, Trace) and isinstance(trP, Trace) and not isinstance(tr1, Trace) or isinstance(trZ, Trace) and isinstance(trP, Trace) and not isinstance(tr2, Trace) or isinstance(trZ, Trace) and isinstance(tr1, Trace) and not isinstance(trP, Trace) or isinstance(trZ, Trace) and isinstance(tr2, Trace) and not isinstance(trP, Trace):
            print("Problem, missing at least one channel, HH1, HH2, or HDH at least one of the others is present, save without TC removal")
            newtags_tmp=[]
            OBS_uncor.append(ista)
            outstream=Stream(traces=[trZ])
            #print(len(outstream))
            newtags_tmp.append(utils.get_tracetag(trZ))
            #print(len(newtags_tmp))
            utils.save2asdf(df_tc, outstream, newtags_tmp, sta_inv=inv)
        if not isinstance(tr1, Trace) and not isinstance(tr2, Trace) and not isinstance(trZ, Trace):
                print('  No seismic channels found. Drop the station: '+ista)
                continue
        if None not in [tr1,tr2,trZ]:
                if tr1.data.shape != tr2.data.shape or tr2.data.shape != trZ.data.shape or tr1.data.shape != trZ.data.shape:
                        badtrace = True
                        print(' Traces have different shapes, skipped ')
        for tr in [tr1, tr2, trZ]:
            if not isinstance(tr, Trace):
                print("  "+str(tr)+" is not a Trace object. "+ista)
                badtrace=True
                break
            elif np.sum(np.isnan(tr.data))>0:
                print('  NaN found in trace: '+str(tr)+". "+ista)
                badtrace=True
                break
            elif np.count_nonzero(tr.data) < 1:
                print('  All zeros in trace: '+str(tr)+". "+ista)
                badtrace=True
                break
        if badtrace:
            if not drop_if_has_badtrace:
                print("  Not enough good traces for TC removal! Save as is without processing!")
                outtrace=[]
                for tg in all_tags:
                    outtrace.append(ds.waveforms[ista][tg][0])
                utils.save2asdf(df_tc,Stream(traces=outtrace),all_tags,sta_inv=inv)
            else:
                print("  Encountered bad trace for "+ista+". Skipped!")
            continue
        elif requirePressure and not hasPressure: #if station doesn't have pressure channel, it might be an obs or a land station
            newtags_tmp=[]
            if isinstance(tr1, Trace) and isinstance(tr2, Trace) and correct_obs_orient and ista in obs_orient_data.keys():
                #correct horizontal orientations if in the obs_orient_data list.
                print("  Correcting horizontal orientations for: "+ista)
                trE,trN = obs.correct_orientations(tr1,tr2,obs_orient_data)
                #newtags_tmp.append(utils.get_tracetag(trE))
                #newtags_tmp.append(utils.get_tracetag(trN))
                #print(newtags_tmp)
            else: #save the station as is if it is not in the orientation database, assuming it is a land station.
                print("Missing one of the horizontsal components")
                #newtags_tmp.append(utils.get_tracetag(tr1))
                #newtags_tmp.append(utils.get_tracetag(tr2))
            newtags_tmp.append(utils.get_tracetag(trZ))
            outstream=Stream(traces=[trZ])
            print('  Saving '+ista+' without TC removal to: '+df_tc)
            utils.save2asdf(df_tc,outstream,newtags_tmp,sta_inv=inv)
            OBS_uncor.append(ista)
            continue

        """
        Call correction wrapper
        """
        try:
            spectra,transfunc,correct=obs.TCremoval_wrapper(
                tr1,tr2,trZ,trP,window=window,overlap=overlap,merge_taper=taper,
                qc_freq=qc_freq,qc_spectra=True,fig_spectra=False,
                save_spectrafig=False,fig_transfunc=False,correctlist=tc_subset)
            tilt.append(spectra['rotation'].tilt)
            sta_processed.append(ista)
            if plot_correction:
                obs.plotcorrection(trZ,correct,normalize=normalizecorrectionplot,freq=[0.005,0.1],
                                   size=(12,3),save=True,form='png')

            trZtc,tgtemp=obs.correctdict2stream(trZ,correct,tc_subset)
            if correct_obs_orient:
                print("  Correcting horizontal orientations for: "+ista)
                trE,trN = obs.correct_orientations(tr1,tr2,obs_orient_data)
                newtags=[]
                #newtags[0]=utils.get_tracetag(trE)
                #newtags[1]=utils.get_tracetag(trN)
                newtags.append(utils.get_tracetag(trZ))
                print(newtags)
                outstream=Stream(traces=trZtc[0])   #[trP,trN,trE,trZtc[0]])
                sta_cor.append(ista)
            else:
                newtags=[]
                newtags.append(utils.get_tracetag(trZ))
                outstream=Stream(traces=trZtc[0])
            """
            Save to ASDF file.
            """
            print('  saving to: '+df_tc)
            utils.save2asdf(df_tc,outstream,newtags,sta_inv=inv)
        except Exception as e:
            print(' Error in calling TCremoval procedures. Save uncorrected trace.')
            print(df+' : '+ista+' : '+str(e))
            outstream=Stream(traces=trZ)
            #print(len(outstream))
            nametags=[]
            nametags.append(utils.get_tracetag(trZ))
            #print(len(nametags))
            utils.save2asdf(df_tc,outstream,nametags,sta_inv=inv)
            continue
    #save auxiliary data to file.
    if len(tilt) > 0:
        print('  saving auxiliary data to: '+df_tc)
        tcpara_temp=tcpara
        tcpara_temp['tilt_stations']=sta_processed
        utils.save2asdf(df_tc,np.array(tilt),None,group='auxiliary',para={'data_type':'tcremoval',
                                                    'data_path':'tiltdir',
                                                    'parameters':tcpara_temp})

###############################################
comm.barrier()
if rank == 0:
    tend=time.time() - t0
    print('*************************************')
    print('<<< Finished all files in %7.1f seconds, or %6.2f hours for %d files >>>' %(tend,tend/3600,len(dfiles)))
    sys.exit()


print("OBS corrected for: ", sta_cor)
print("OBS not corrected for: ", OBS_uncor)
print("Land stations: ", sta_uncor)
print("Process ended successfully")

fout = open(OBS_cor_out,'w')
fout.write('\n'.join(sta_cor));
fout.close()

fout = open(OBS_uncor_out,'w')
fout.write('\n'.join(OBS_uncor));
fout.close()

fout = open(land_out,'w')
fout.write('\n'.join(sta_uncor));
fout.close()
