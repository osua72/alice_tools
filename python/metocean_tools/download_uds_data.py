#!/usr/bin/env python

# This works on Aotea3 or metocean docker' ops-base e.g.
# 
# docker run -it -v /data2:/data2 metocean/ops-base

import sys 
sys.path.append('/data2/simon/github/toolbox_simon/python_tools') 
import download_metocean_uds
import datetime

years = range(2000,2011)
months = range(1,13)
#output path
path_data  = '/data2/alice/0495_beachenergy_oilspill/' # location where queried data will be saved

# coverage same as MOANA backbone BL=[161.03 -51.9812] UR=[184.97 -31.0261]
frame_global = [161.03,184.97,-51.9812,-31.0261] # for  ww3_onr_2520 and  era5_wnd10m

frame_swan = [165,175,-48,-42] #for swan_nzra_nz-nzs BL=[165 -48] UR=[175 -42] hs tp dpm
# frame_nzra = [163.279,182.793,-48.297,-32.697] # for nzra1_nz  [163.279 -48.297] UR=[182.793 -32.697] 1hours
frame_nzra = [165,180.0,-47.5,-32.697] # for nzra1_nz  [163.279 -48.297] UR=[182.793 -32.697] 1hours need to adjust to ensure we get only good data

# hydro data is read straight from the ROMS moana backbone

for iy,yr in enumerate(years):
    for imth,mth in enumerate(months):
        tstart= datetime.datetime(yr,mth,1)
        if mth == 12:
            tend = datetime.datetime(yr+1,1,1)
        else:
            tend = datetime.datetime(yr,mth+1,1)
        
        # WINDS ###################################################################################################
        download_metocean_uds.download(uds_host = 'http://uds1.rag.metocean.co.nz:9191/uds',
                 datasets=['nzra1_nz'],
                 variables=['ugrd10m','vgrd10m'],
                 time_start=tstart,
                 time_end=tend,
                 timestep = 3.0,
                 boundary=frame_nzra,
                 datatype='hc')

        download_metocean_uds.download(uds_host = 'http://uds1.rag.metocean.co.nz:9191/uds',
                 datasets=['era5_wnd10m'],
                 variables=['ugrd10m','vgrd10m'],
                 time_start=tstart,
                 time_end=tend,
                 timestep = 3.0,
                 boundary=frame_global,
                 datatype='hc')

        # WAVES ####################################################################################################

        download_metocean_uds.download(uds_host = 'http://uds1.rag.metocean.co.nz:9191/uds',
                 datasets=['swan_nzra_nz-nzs'],
                 variables=['hs','tp','dpm'],
                 time_start=tstart,
                 time_end=tend,
                 timestep = 3.0,
                 boundary=frame_swan,
                 datatype='hc')

        download_metocean_uds.download(uds_host = 'http://uds1.rag.metocean.co.nz:9191/uds',
                 datasets=['ww3_onr_2520'],
                 variables=['hs','tp','dpm'],
                 time_start=tstart,
                 time_end=tend,
                 timestep = 3.0,
                 boundary=frame_global,
                 datatype='hc')