#!/usr/bin/env python
#
# Python Class to undertak all post-processing operation on OpenDrift netcdf output files
# 
# 
import os
import numpy as np
import argparse
import sys
import glob
import importlib
import gc
#sys.path.append('/home/agowardbrown/Documents/alice_tools/python')
# import download_metocean_uds
# import download_cmems
from datetime import datetime,timedelta
import opendrift
import logging
import xarray
import math
# for plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LogNorm
import shapefile as shp

from mpl_toolkits.basemap import Basemap
from osgeo import gdal

import datetime as dt

import utils
import utm

import cmocean 
import cmocean.cm as cmo

def calc_vol(dredgeordisp,OF_time,dredge,nb_parts,vol2rem,sites,sed,sed_frac):
    # python calc_vol [1 (dredging) / 2 (disposing)] [10/20/30(mins)] 
    # [TSHD,TSHD_OF,BHD] [# particles released per cycle] 
    # [volume of seidment to remove per site] [# sites] [sediment type] [sediment fraction]
    if dredge == 'TSHD' or dredge == 'TSHD_OF':
        dredge_rate = 117600 #m3/week
        dredge_capacity = 635.
        cycle_time = 3.5 # h
        loading_time = (cycle_time - 1.0)*3600 # loading time in seconds
        if dredge == 'TSHD_OF':
            overflow_time = OF_time*60 # overflow time in seconds
    elif dredge == 'BHD':
        dredge_rate = 50000 # m3/week
        dredge_capacity = 1860.
        cycle_time = 24 # h
        loading_time = (cycle_time - 1.0)*3600 # loading time in seconds

    if dredgeordisp == 1:
        vol2rem_site = vol2rem/sites
        print('Volume to remove per site (m3)', vol2rem_site)

        loading_time = (cycle_time - 1.0)*3600 # loading time in seconds
        print('Loading Time (s): ', loading_time)

        placement_duration = 10*60
        print('Placement Duration (s): ', placement_duration)

        if sed == 'sand':
            dry_weight = 1600 # kg/m3
        else: dry_weight = 500 # kg/m3
        print('Dry Weight (kg/m3)', dry_weight)

        weeks = vol2rem_site/dredge_rate
        print('Days per site: ', weeks * 7)

        cycles = (weeks*(7*24*3600))/loading_time
        print('# Cycles', cycles)

        total_loading_time = cycles*loading_time
        print('Total loading time: ', total_loading_time)

        Drate = dredge_rate*(7*24*3600) # m3/s
        Prate = vol2rem_site/total_loading_time # m3/s

        vol_per_cycle = np.round(loading_time*Prate) # m3
        print('Volume to dispose per cycle (m3)', vol_per_cycle) 

    elif dredgeordisp == 2:
        if sed == 'sand':
            dry_weight = 1600. # kg/m3
        else: dry_weight = 500. # kg/m3
        print('Dry Weight (kg/m3)', dry_weight)
        #    
        vol_per_cycle = dredge_capacity
        print('Volume to dispose per cycle (m3)', vol_per_cycle)    

    mass_per_cycle = vol_per_cycle*dry_weight*sed_frac # kg 
    print('mass of',sed,'removed per cycle: ',mass_per_cycle) # kg

    if dredgeordisp == 1:
        if dredge == 'TSHD_OF':
            pdh0 = 3/100 # 3 % of fines
            ss_mass_1 = mass_per_cycle*pdh0
            mass_hopper_1 = mass_per_cycle-ss_mass_1
            R0 = (overflow_time)/(loading_time)
            print('R0', R0)
            fsett = 0.25
            ftrap = 0.05
            mass_OF = R0*(1-fsett)*(1-ftrap)*mass_hopper_1 #check! this isn't a mass..
            pp0 = 20/100 # 20% of hopper mass OF = fines entrained
            ss_mass_2 = mass_OF*pp0 # OF plume kg
            ss_mass_3 = (1-pp0)*mass_OF # 80% in density current kg
            mass_sum = ss_mass_1+ss_mass_2+ss_mass_3
        if dredge == 'TSHD':
            pdh0 = 3/100 # 3% as fines
            ss_mass_1 = mass_per_cycle*pdh0
            mass_sum = ss_mass_1

        if dredge == 'BHD':
            pdh0 = 4/100 #
            ss_mass_1 = mass_per_cycle*pdh0 # bucket drip fraction
            mass_sum = ss_mass_1
        

    if dredgeordisp == 2:
        if dredge == 'TSHD' or dredge == 'TSHD_OF':
            pdh0 = (100.-73.)/100. #%
            ss_mass_1 = mass_per_cycle*pdh0

        if dredge == 'BHD':
            pdh0 = (100.-73.)/100. #% 23 % available for OD modelling
            ss_mass_1 = mass_per_cycle*pdh0 

        mass_sum = ss_mass_1
    print('Total mass of sed removed per cycle', mass_per_cycle)
    print('Fraction of sediment available for plume disp: ',pdh0,'%')
    print('Total mass released per cycle: ', mass_sum) 
    
    mass_per_part = mass_sum/nb_parts
    print('Mass per particle (kg): ', mass_per_part) 
    return mass_sum

def get_points(lon,lat,time,pdf_s,pdf_m,pdf_b,coords):
        # find where lon | lat == coords
        pdf_s_points = np.zeros([len(time),np.shape(coords)[1]])
        pdf_m_points = np.zeros([len(time),np.shape(coords)[1]])
        pdf_b_points = np.zeros([len(time),np.shape(coords)[1]])

        for i in range(0,np.shape(coords)[1]):
            print(i)
            lonidx = int([ti for ti, x in enumerate(np.abs(lon - coords[0][i])) if x == np.min(np.abs(lon - coords[0][i]))][0])
            latidx = int([ti for ti, x in enumerate(np.abs(lat - coords[1][i])) if x == np.min(np.abs(lat - coords[1][i]))][0])
            # get PDF
            pdf_s_points[:,i] = pdf_s[:,int(latidx),int(lonidx)]
            pdf_m_points[:,i] = pdf_m[:,int(latidx),int(lonidx)]
            pdf_b_points[:,i] = pdf_b[:,int(latidx),int(lonidx)]

        print(np.amax(pdf_b_points))

        return pdf_s_points,pdf_m_points,pdf_b_points

def calc_ssc(fdir,sed,nb_parts,method,dredge,key,switch):
    # days to allow settling to accumulate per area hardcoded for now:    
    # find all files which are for the sediment [key] of interest
    flist = [s for s in fdir if key in s]
    # get the fdir for saving later
    fdir,file = os.path.split(flist[0])
    print(fdir)
    print(flist)
    count = 0
    # read in daya and convert to kg/cell
    for file in flist:
        if key in file: # look for all files contianing the type of sediment 
            print(file)
            ds = xarray.open_dataset(file)
            #nb_parts = ds.attrs['nb_part_total']
            # calculate particle load for each sediment class / scenario
            mass_per_part = sed/nb_parts #kg to kg/particle (released per cycle)
            if key == 'sand':
                vol_per_part = mass_per_part / 1600.0 # kg/particle / kg / m3 = m3/particle
            elif key == 'silt' or key == 'clay':
                vol_per_part = mass_per_part / 500.0

            print('Mass per part: ', mass_per_part)
            print('Volume per part: ', vol_per_part)

            # get lon / lat / water depth
            if count == 0:
                lon = ds['lon']
                lat = ds['lat']
                h = ds['water_depth']
                ntime = ds['time'].shape[0]
                time = ds['time']
                if method == 'SSC':
                    pdf_active_s = ds['pdf_active[0_-3]'] * mass_per_part #kg/cell
                    pdf_active_m = ds['pdf_active[mid+1_mid-1]'] * mass_per_part     
                    pdf_active_b = ds['pdf_active[seafloor+2_seafloor]'] * mass_per_part     
                    nb_part_act = ds['nb_part_active'] * mass_per_part    
                elif method == 'DEP':
                    dt = (ds['time'][1].dt.minute-ds['time'][0].dt.minute)*60.0
                    last_step = ds['pdf_settled'][-1,:,:]
                    # [m3/cell] 
                    pdf_settled = last_step * vol_per_part 
                    #
                    nb_part_settled = ds['nb_part_settled'][-1] * vol_per_part                            
                count +=1
            else:
                h = np.mean([h,ds['water_depth']])
                ntime = ntime+(ds['time'].shape[0])
                time = np.concatenate([time,ds['time']])
                if method == 'SSC':
                    pdf_active_s = np.concatenate([pdf_active_s,ds['pdf_active[0_-3]'] * mass_per_part],axis=0)
                    pdf_active_m = np.concatenate([pdf_active_m,ds['pdf_active[mid+1_mid-1]'] * mass_per_part],axis=0)
                    pdf_active_b = np.concatenate([pdf_active_b,ds['pdf_active[seafloor+2_seafloor]'] * mass_per_part],axis=0)     
                    nb_part_act = np.concatenate([nb_part_act,ds['nb_part_active'] * mass_per_part])
                elif method == 'DEP':        
                    # gets the last timestep for the appropriate length of time it will take to dredge the area specified
                    last_step = ds['pdf_settled'][-1,:,:]
                    # [m3/cell] 
                    pdf_settled = np.dstack([pdf_settled,last_step * vol_per_part])
                    last_step = ds['nb_part_settled'][-1]
                    nb_part_settled = np.append(nb_part_settled,last_step * vol_per_part)

    # convert pdf from [kg/cell] to [kg/m2]
    cell_area = ds.attrs['pixelsize_m'] ** 2 # use last open dataset (same cell area for all)
    print('cell_area: ', cell_area)
    for var in ds.data_vars:
        if method == 'SSC' and 'pdf' in var: # pdf fields
            # normalize by cell area to get [kg/m2] for TSS, or [m3/m2]=[m] for deposition
            if 'active' in var:
                # scan the vertical depth band considered
                vert_band = var[1+var.find('['):var.find(']')]
                vert_band = vert_band.replace('_',',')
                vert_band = vert_band.replace('seafloor','0')
                vert_band = vert_band.replace('mid','0')
                if vert_band != 'all' : # cover all cases except 'all'
                    # define depth_band_thickness
                    vert_band = eval(vert_band) 
                    depth_band_thickness  = vert_band[0] - vert_band[1] # it is expected that depth is given as [shallower,deeper]                    
                else:
                    if 'water_depth' in ds:
                        # tile water depth matrix in time-dimension to fit with xr_dataset[var][:,:,:] dimensions
                        depth_band_thickness = np.tile(ds['water_depth'],[ds['time'].shape[0],1,1])
                    else :
                        print('Water depth information for processing level %s' % (vert_band))
                        print('Add water depth info using method : add_waterdepth_to_dataset()  ')
                        print('e.g. ds_combined = add_waterdepth_to_dataset (xr_dataset = ds_combined,config_obj = config_obj) ')
                        return

                # normalize by depth band thickness to get to kg/m3
                
                if vert_band == (0,-3): #surface
                    print(var,cell_area)
                    print(vert_band, depth_band_thickness)
                    #print('pdf_max: ', max(max(max(pdf_active_s))))
                    #print('pdf_min: ', min(min(min(pdf_active_s))))
                    pdf_active_s_m2 = pdf_active_s / cell_area # kg/m2
                    pdf_active_s_m3 = pdf_active_s_m2 / depth_band_thickness 
                elif vert_band == (0+1,0-1): # mid
                    print(var,cell_area)
                    print(vert_band, depth_band_thickness)
                    pdf_active_m_m2 = pdf_active_m / cell_area 
                    pdf_active_m_m3 = pdf_active_m_m2 / depth_band_thickness 
                elif vert_band == (0+2,0): # bottom
                    print(var,cell_area)
                    print(vert_band,depth_band_thickness)
                    pdf_active_b_m2 = pdf_active_b / cell_area 
                    pdf_active_b_m3 = pdf_active_b_m2 / depth_band_thickness    

        elif method == 'DEP' and 'settled' in var or 'stranded' in var :
            print(var,cell_area)
            pdf_settled = pdf_settled / cell_area #m
                # no normalization by vertical depth band required          
   
    if method == 'SSC':
        print(np.shape(pdf_active_s_m3))
        return lon,lat,h,pdf_active_s_m3,pdf_active_m_m3,pdf_active_b_m3
    if method == 'DEP':
        return lon,lat,h,pdf_settled    

def ssc_ts (fdir,mass_sand,mass_silt,mass_clay,nb_parts,hoc,dredge,coords,method,stats,switch):

    def make_ts_plot(time,ssc_s,ssc_m,ssc_b,dredge):
        print('Plotting TS!')           
        # make plt
        d0 = dt.datetime(1970,1,1,0,0,0)
        date = []
        for t in time: 
            date.append(d0+dt.timedelta(seconds = t))

        t = date[0:len(ssc_s[:,0])]
        # organise figure
        #fig = plt.figure(figsize=(13,8))  # a new figure window
        fig, (ax1,ax2,ax3) = plt.subplots(3, figsize=(13,8), sharex=True,gridspec_kw={'hspace': 0.1})

        print(np.shape(t),np.shape(ssc_s))
        print(np.max(ssc_s))
        c1=ax1.plot(t,ssc_s[:])
        ax1.legend((c1), ('site1', 'site2','site3','site4','site5','site6','site7','site8'),
        loc='upper right',ncol=4,fancybox=True)
        ax1.text(0.05, 0.95, 'surface',
        verticalalignment='top', horizontalalignment='left',
        transform=ax1.transAxes,
        color='black', fontsize=11)
        ax1.set_ylabel('SSC [mg/L]',fontsize = 14)
        ax1.grid(True)
        
        #
        c2=ax2.plot(t,ssc_m[:])
        ax2.legend((c2), ('site1', 'site2','site3','site4','site5','site6','site7','site8'),
        loc='upper right',ncol=4,fancybox=True)
        ax2.text(0.05, 0.95, 'mid-depth',
        verticalalignment='top', horizontalalignment='left',
        transform=ax2.transAxes,
        color='black', fontsize=11)
        ax2.set_ylabel('SSC [mg/L]',fontsize = 14)
        ax2.grid(True)
        #
        c3=ax3.plot(t,ssc_b[:])
        ax3.legend((c3), ('site1', 'site2','site3','site4','site5','site6','site7','site8'),
        loc='upper right',ncol=4,fancybox=True)
        ax3.text(0.05, 0.95, 'bottom',
        verticalalignment='top', horizontalalignment='left',
        transform=ax3.transAxes,
        color='black', fontsize=11)
        ax3.set_ylabel('SSC [mg/L]',fontsize = 14)
        ax3.grid(True)

        # display
        print('Saving figure ')
        outname = 'SSC_%s.png' %(dredge)
        plt.savefig(outname,bbox_inches = 'tight',
        pad_inches = 0.1)

    count = 0
    for key in ['sand','silt','clay']:
        fname_s = os.path.join(fdir,'opendrift_nelson_%s_%s_*.nc' % (dredge,key))
        flist = sorted(glob.glob(fname_s))
        for file in flist:
            dx = int(file[len(file)-5:len(file)-3])
            cell_volume = dx * dx * hoc   
            cell_area = dx * dx        

        if method == 'SSC':
            [lon,lat,h,time,pdf_s,pdf_m,pdf_b,nb_part_act]=get_vars(flist,method,stats,switch)
            [pdf_s,pdf_m,pdf_b]=get_points(lon,lat,time,pdf_s,pdf_m,pdf_b,coords)
            print(np.max(pdf_s))

        if count == 0:
            ssc_s_all = np.zeros([len(time),np.shape(coords)[1]])
            ssc_m_all = np.zeros([len(time),np.shape(coords)[1]])
            ssc_b_all = np.zeros([len(time),np.shape(coords)[1]])
            count+=1

        #SC calc            
        # get weight per particle
        if key == 'sand':
            mass_per_part = mass_sand/nb_parts
            print('Mass per particle (kg): ',mass_per_part)
        elif key == 'silt':
            mass_per_part = mass_silt/nb_parts
        elif key == 'clay':
            mass_per_part = mass_clay/nb_parts
                    
        # SSC
        SSC_b = do_calc(cell_volume,mass_per_part,pdf_b)
        SSC_m = do_calc(cell_volume,mass_per_part,pdf_m)
        SSC_s = do_calc(cell_volume,mass_per_part,pdf_s)

        print(key,'max SSC b / s: ',np.amax(SSC_b),np.amax(SSC_s))
   
        print(np.shape(pdf_b))
        print(np.shape(pdf_s))
        print(np.shape(SSC_b))
        
        ssc_s_all = ssc_s_all+SSC_s[0:len(ssc_s_all[:,0]),:]
        ssc_m_all = ssc_m_all+SSC_m[0:len(ssc_m_all[:,0]),:]
        ssc_b_all = ssc_b_all+SSC_b[0:len(ssc_b_all[:,0]),:]

    make_ts_plot(time,ssc_s_all,ssc_m_all,ssc_b_all,dredge)

def make_map(flist = None,
             statslist = None,
             method = None,
             topo = None,
             dredge_area = None,
             dredge=None,
             area = None,
             bulk_vol = None):
    tickMinutes = .5
    def plotCommon(ax):
        ax.xaxis.set_major_locator(ticker.MultipleLocator(tickMinutes/60.0))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(tickMinutes/60.0))
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(utils.formatDegreesE))
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(utils.formatDegreesS))

        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(11)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(11)

    def var2nan(values,var):
        """Replace every 0 with 'nan' and return a copy."""
        values[ values<=var ]=np.nan
        return values


    for stat in statslist:
        print(stat)
        outname = '%s_%s_%s_%s.png'%(method,dredge,stat,area)
        count = 0
        # organise figure
        if method == 'DEP':
            fig, axs = plt.subplots(nrows=len(flist),figsize=(8.27,11.69), sharex=True,gridspec_kw={'hspace': 0.2})
            tickMinutes = .7
            vmin = 0.0001 # .1mm
            vmax = 2.0
            dredge=160

        elif method == 'SSC':     
            #fig = plt.figure(figsize=(8.27,11.69))
            fig, axs = plt.subplots(len(flist),3,figsize=(8.27,11.69), sharey=True,sharex=True,gridspec_kw={'hspace': 0.05})
            tickMinutes = .8

            if dredge == 'TSHD':
                if stat == 'mean':
                    vmin = 1.0
                    vmax = 100.
                elif stat == 'perc':
                    vmin = 1.0
                    vmax = 1000.
                elif stat == 'P90':
                    vmin = 0.0
                    vmax = 1000.
            if dredge == 'BHD':
                if stat == 'mean':
                    vmin = 1.0
                    vmax = 100.
                elif stat == 'perc':
                    vmin = 1.0
                    vmax = 1000.        

        #import pdb;pdb.set_trace()    
        for i_file,file in enumerate(flist):
            print(count,file)
            ds = xarray.open_dataset(file)
            
            if method == 'SSC':
                if 'exceed' in stat:
                    count = 0
                ax1 = axs[count][0]
                ax2 = axs[count][1]
                ax3 = axs[count][2]
                #except: import pdb;pdb.set_trace()
            elif method == 'DEP':
                ax1 = axs[count]

            if topo is not None:
                datafile = gdal.Open(topo)
                bnd1 = datafile.GetRasterBand(1).ReadAsArray()
                bnd2 = datafile.GetRasterBand(2).ReadAsArray()
                bnd3 = datafile.GetRasterBand(3).ReadAsArray()
                nx = datafile.RasterXSize # Raster xsize
                ny = datafile.RasterYSize # Raster ysize
                #
                img = np.dstack((bnd1, bnd2, bnd3))
                gt = datafile.GetGeoTransform()
                proj = datafile.GetProjection()

                print("Geotransform",gt)
                print("proj=", proj)
                xres = gt[1]
                yres = gt[5]

                # get the edge coordinates and add half the resolution 
                # to go to center coordinates
                xmin = gt[0] + xres * 0.5
                xmax = gt[0] + (xres * nx) - xres * 0.5
                ymin = gt[3] + (yres * ny) + yres * 0.5
                ymax = gt[3] - yres * 0.5
                print("xmin=", xmin,"xmax=", xmax,"ymin=",ymin, "ymax=", ymax)
                # create a grid of lat/lon coordinates in the original projection
                (lon_source,lat_source) = np.mgrid[xmin:xmax+xres:xres, ymax+yres:ymin:yres]
                print(xmin,xmax+xres,xres, ymax+yres,ymin,yres)
                if method == 'SSC':
                    m1 = Basemap(projection='cyl',llcrnrlat=ymin,urcrnrlat=ymax,
                        llcrnrlon=xmin,urcrnrlon=xmax , resolution='f', ax=ax1)
                    m2 = Basemap(projection='cyl',llcrnrlat=ymin,urcrnrlat=ymax,
                        llcrnrlon=xmin,urcrnrlon=xmax , resolution='f', ax=ax2)  
                    m3 = Basemap(projection='cyl',llcrnrlat=ymin,urcrnrlat=ymax,
                        llcrnrlon=xmin,urcrnrlon=xmax , resolution='f', ax=ax3)                
        
                    # project in the original Basemap

                    x,y = m1(lon_source, lat_source)
                    print("shape lon and lat_source: ", lon_source.shape, lat_source.shape,bnd1.T.shape)
                    m1.imshow(img, origin='upper', ax=ax1)
                    m2.imshow(img, origin='upper', ax=ax2)
                    m3.imshow(img, origin='upper', ax=ax3)

                elif method == 'DEP':
                    #basemap
                    map = Basemap(projection='cyl',llcrnrlat=ymin,urcrnrlat=ymax,
                    llcrnrlon=xmin,urcrnrlon=xmax , resolution='f', ax=ax1)
                    # project in the original Basemap
                    x,y = map(lon_source, lat_source)
                    print("shape lon and lat_source: ", lon_source.shape, lat_source.shape,bnd1.T.shape)
                    map.imshow(img, origin='upper', ax=ax1)

            if method == 'SSC' and 'exceed' in stat:
                thresh = [10,50,100]
                fig, axs = plt.subplots(len(thresh),3,figsize=(8.27,8.27), sharey=True,sharex=True,gridspec_kw={'hspace': 0.01})
                count = 0
                vmin = 1.0
                vmax = 25.

                # start plotting
                for it,thr in enumerate(thresh):
                    print(it,count)

                    # sea x and y labels
                    ax1 = axs[count][0]
                    ax2 = axs[count][1]
                    ax3 = axs[count][2]
                    ax1.set_ylabel('Latitude $^{\circ}$S',fontsize = 14)

                    if count == 0:                   
                        # set title
                        ax1.set_title('Surface',fontsize=14)
                        ax2.set_title('Middle',fontsize=14)
                        ax3.set_title('Bottom',fontsize=14)

                    # set x and y limits
                    #import pdb;pdb.set_trace()
                    ax1.set_xlim(np.min(ds['lon']),np.max(ds['lon']))
                    ax1.set_ylim(np.min(ds['lat']), np.max(ds['lat']))
                    ax2.set_xlim(np.min(ds['lon']),np.max(ds['lon']))
                    ax2.set_ylim(np.min(ds['lat']), np.max(ds['lat']))
                    ax3.set_xlim(np.min(ds['lon']),np.max(ds['lon']))
                    ax3.set_ylim(np.min(ds['lat']), np.max(ds['lat']))
                    #
                    ax1.set_aspect('equal')
                    ax2.set_aspect('equal')
                    ax3.set_aspect('equal')

                    if count == len(thresh)-1:
                        ax1.set_xlabel('Longitude $^{\circ}$E',fontsize = 14)
                        ax2.set_xlabel('Longitude $^{\circ}$E',fontsize = 14)  
                        ax3.set_xlabel('Longitude $^{\circ}$E',fontsize = 14)

                    stat = 'exceed_'+str(thr)
                    var = ds['ssc_'+stat][0,:,:]
                    
                    var1 = var2nan(np.array(var),0.0)
                    c = ax1.pcolormesh(ds['lon'],ds['lat'],var1,vmin=vmin,vmax=vmax,cmap=cmo.matter)
                    ax1.plot(ds.attrs['center_point'][0], ds.attrs['center_point'][1],marker='.',markersize=1.0,color='w')
                    c.set_clim(vmin, vmax)
                    # circle / radius around point
                    draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),250/111000,color='k',fill=False,linestyle='--')
                    ax1.add_artist(draw_circle)
                    draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),500/111000,color='k',fill=False,linestyle='--')
                    ax1.add_artist(draw_circle)
                    draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),1000/111000,color='k',fill=False,linestyle='--')
                    ax1.add_artist(draw_circle)
                    #m1.readshapefile(dredge_area, 'dredgeA')
                    plotCommon(ax1)
                    # plot2
                    var = ds['ssc_'+stat][1,:,:]
                    var2 = var2nan(np.array(var),0.0)
                    c = ax2.pcolormesh(ds['lon'],ds['lat'],var2,vmin=vmin,vmax=vmax,cmap=cmo.matter)
                    ax2.plot(ds.attrs['center_point'][0], ds.attrs['center_point'][1],marker='.',markersize=1.0,color='w')
                    c.set_clim(vmin, vmax)
                    # circle / radius around point
                    draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),250/111000,color='k',fill=False,linestyle='--')
                    ax2.add_artist(draw_circle)
                    draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),500/111000,color='k',fill=False,linestyle='--')
                    ax2.add_artist(draw_circle)
                    draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),1000/111000,color='k',fill=False,linestyle='--')
                    ax2.add_artist(draw_circle)
                    #m2.readshapefile(dredge_area, 'dredgeA')
                    plotCommon(ax2)
                    # plot 3
                    var = ds['ssc_'+stat][2,:,:]
                    var3 = var2nan(np.array(var),0.0)
                    c = ax3.pcolormesh(ds['lon'],ds['lat'],var3,vmin=vmin,vmax=vmax,cmap=cmo.matter)
                    ax3.plot(ds.attrs['center_point'][0], ds.attrs['center_point'][1],marker='.',markersize=1.0,color='w')
                    c.set_clim(vmin, vmax)
                    # circle / radius around point
                    draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),250/111000,color='k',fill=False,linestyle='--')
                    ax3.add_artist(draw_circle)
                    draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),500/111000,color='k',fill=False,linestyle='--')
                    ax3.add_artist(draw_circle)
                    draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),1000/111000,color='k',fill=False,linestyle='--')
                    ax3.add_artist(draw_circle)
                    #m3.readshapefile(dredge_area, 'dredgeA')
                    #colorbar # hard coded..
                    if len(thresh) == 3:
                        if count == 0:
                            cax = fig.add_axes([0.92, 0.65, 0.02, 0.2])
                        elif count == 1:
                            cax = fig.add_axes([0.92, 0.395, 0.02, 0.2])
                        elif count ==2:
                            cax = fig.add_axes([0.92, 0.135, 0.02, 0.2])
    
                    cbar = fig.colorbar(c, cax=cax)
                    cbar.ax.set_ylabel('% time > '+str(thr)+'[mg/L]',size=11) 
                    #map.readshapefile(dredge_area, 'dredgeA')
                    plotCommon(ax3)
                    count+=1
                
                # display
                print('Saving figure ')
                print(i_file)
                outname = '%s_%s_%s_%s.png'%(method,dredge,i_file,area)
                plt.savefig(outname,bbox_inches = 'tight',
                pad_inches = 0.1)    

            if method == 'SSC' and stat == 'mean':
                # sea x and y labels
                ax1.set_ylabel('Latitude $^{\circ}$S',fontsize = 14)

                if count == 0:                   
                    # set title
                    ax1.set_title('Surface',fontsize=14)
                    ax2.set_title('Mid-water',fontsize=14)
                    ax3.set_title('Bottom',fontsize=14)

                if count == len(flist)-1:
                    ax1.set_xlabel('Longitude $^{\circ}$E',fontsize = 14)
                    ax2.set_xlabel('Longitude $^{\circ}$E',fontsize = 14)  
                    ax3.set_xlabel('Longitude $^{\circ}$E',fontsize = 14)

                # set x and y limits
                #import pdb;pdb.set_trace()
                ax1.set_xlim(np.min(ds['lon']),np.max(ds['lon']))
                ax1.set_ylim(np.min(ds['lat']), np.max(ds['lat']))
                ax2.set_xlim(np.min(ds['lon']),np.max(ds['lon']))
                ax2.set_ylim(np.min(ds['lat']), np.max(ds['lat']))
                ax3.set_xlim(np.min(ds['lon']),np.max(ds['lon']))
                ax3.set_ylim(np.min(ds['lat']), np.max(ds['lat']))
                #
                ax1.set_aspect('equal')
                ax2.set_aspect('equal')
                ax3.set_aspect('equal')               
                # start plotting
                # plot1 
                #import pdb;pdb.set_trace()
                #ds['ssc_'+stat][0,:,:]<=vmin = nan #var1 =var2nan(var1,vmin) 
                var = ds['ssc_'+stat][0,:,:]
                var1 = var2nan(np.array(var),vmin)
                c = ax1.pcolormesh(ds['lon'],ds['lat'],var1,vmin=vmin,vmax=vmax,cmap=cmo.matter)
                ax1.plot(ds.attrs['center_point'][0], ds.attrs['center_point'][1],marker='.',markersize=1.0,color='w')
                c.set_clim(vmin, vmax)
                # circle / radius around point
                draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),250/111000,color='k',fill=False,linestyle='--')
                ax1.add_artist(draw_circle)
                draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),500/111000,color='k',fill=False,linestyle='--')
                ax1.add_artist(draw_circle)
                draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),1000/111000,color='k',fill=False,linestyle='--')
                ax1.add_artist(draw_circle)
                #m1.readshapefile(dredge_area, 'dredgeA')
                # colorbar
                #cbar = fig.colorbar(c, ax=ax1)
                #cbar.ax.set_ylabel('SSC [mg/L]',size=11) 
                plotCommon(ax1)
                # plot2
                #var2 = var2nan(var2,vmin)
                var = ds['ssc_'+stat][1,:,:]
                var2 = var2nan(np.array(var),vmin)
                c = ax2.pcolormesh(ds['lon'],ds['lat'],var2,vmin=vmin,vmax=vmax,cmap=cmo.matter)
                ax2.plot(ds.attrs['center_point'][0], ds.attrs['center_point'][1],marker='.',markersize=1.0,color='w')
                c.set_clim(vmin, vmax)
                # circle / radius around point
                draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),250/111000,color='k',fill=False,linestyle='--')
                ax2.add_artist(draw_circle)
                draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),500/111000,color='k',fill=False,linestyle='--')
                ax2.add_artist(draw_circle)
                draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),1000/111000,color='k',fill=False,linestyle='--')
                ax2.add_artist(draw_circle)
                #m2.readshapefile(dredge_area, 'dredgeA')
                # colorbar
                #cbar = fig.colorbar(c, ax=ax2)
                #cbar.ax.set_ylabel('SSC [mg/L]',size=11) 
                plotCommon(ax2)
                # plot 3
                #var3 = var2nan(var3,vmin)
                var = ds['ssc_'+stat][2,:,:]
                var3 = var2nan(np.array(var),vmin)
                c = ax3.pcolormesh(ds['lon'],ds['lat'],var3,vmin=vmin,vmax=vmax,cmap=cmo.matter)
                ax3.plot(ds.attrs['center_point'][0], ds.attrs['center_point'][1],marker='.',markersize=1.0,color='w')
                c.set_clim(vmin, vmax)
                # circle / radius around point
                draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),250/111000,color='k',fill=False,linestyle='--')
                ax3.add_artist(draw_circle)
                draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),500/111000,color='k',fill=False,linestyle='--')
                ax3.add_artist(draw_circle)
                draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),1000/111000,color='k',fill=False,linestyle='--')
                ax3.add_artist(draw_circle)
                #m3.readshapefile(dredge_area, 'dredgeA')
                #colorbar # hard coded..
                if len(flist) == 5:
                    if count == 0:
                        cax = fig.add_axes([0.92, 0.75, 0.02, 0.11])
                    elif count == 1:
                        cax = fig.add_axes([0.92, 0.6, 0.02, 0.11])
                    elif count ==2:
                        cax = fig.add_axes([0.92, 0.44, 0.02, 0.11])
                    elif count == 3:
                        cax = fig.add_axes([0.92, 0.29, 0.02, 0.11])
                    elif count == 4:
                        cax = fig.add_axes([0.92, 0.13, 0.02, 0.11])

                cbar = fig.colorbar(c, cax=cax)
                cbar.ax.set_ylabel('SSC [mg/L]',size=11) 
                #map.readshapefile(dredge_area, 'dredgeA')
                plotCommon(ax3)         

            if method == 'SSC' and stat == 'perc':
                #import pdb; pdb.set_trace()
                 # sea x and y labels
                ax1.set_ylabel('Latitude $^{\circ}$S',fontsize = 14)

                if count == 0:                   
                    # set title
                    ax1.set_title('Surface',fontsize=14)
                    ax2.set_title('Mid-water',fontsize=14)
                    ax3.set_title('Bottom',fontsize=14)

                if count == len(flist)-1:
                    ax1.set_xlabel('Longitude $^{\circ}$E',fontsize = 14)
                    ax2.set_xlabel('Longitude $^{\circ}$E',fontsize = 14)  
                    ax3.set_xlabel('Longitude $^{\circ}$E',fontsize = 14)

                # set x and y limits
                #import pdb;pdb.set_trace()
                ax1.set_xlim(np.min(ds['lon']),np.max(ds['lon']))
                ax1.set_ylim(np.min(ds['lat']), np.max(ds['lat']))
                ax2.set_xlim(np.min(ds['lon']),np.max(ds['lon']))
                ax2.set_ylim(np.min(ds['lat']), np.max(ds['lat']))
                ax3.set_xlim(np.min(ds['lon']),np.max(ds['lon']))
                ax3.set_ylim(np.min(ds['lat']), np.max(ds['lat']))
                #
                ax1.set_aspect('equal')
                ax2.set_aspect('equal')
                ax3.set_aspect('equal')               
                # start plotting
                # plot1 
                #P90
                var = ds['ssc_P90'][0,:,:]
                #var1 = var2nan(np.array(var),vmin)
                var1 = np.where(var > 0.0,1.0,var)
                c1 = ax1.contour(ds['lon'],ds['lat'],var1,1,label='P90',colors='grey',linewidths=0.5)
                # P95
                var = ds['ssc_P95'][0,:,:]
                #var2 = var2nan(np.array(var),vmin)
                var2 = np.where(var > 0.0,1.0,var)
                c2 = ax1.contour(ds['lon'],ds['lat'],var2,1,label='P95',colors='black',linewidths=0.5)
                # P99
                var = ds['ssc_P99'][0,:,:]
                #var3 = var2nan(np.array(var),vmin)
                var3 = np.where(var > 0.0,1.0,var)
                c3 = ax1.contour(ds['lon'],ds['lat'],var3,1,label='P99',colors='red',linewidths=0.5)
                #
                ax1.plot(ds.attrs['center_point'][0], ds.attrs['center_point'][1],marker='.',markersize=1.0,color='w')
                # circle / radius around point
                draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),250/111000,color='k',fill=False,linestyle='--')
                ax1.add_artist(draw_circle)
                draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),500/111000,color='k',fill=False,linestyle='--')
                ax1.add_artist(draw_circle)
                draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),1000/111000,color='k',fill=False,linestyle='--')
                ax1.add_artist(draw_circle)
                #m1.readshapefile(dredge_area, 'dredgeA')
                plotCommon(ax1)
                # plot2
                var = ds['ssc_P90'][1,:,:]
                var1 = np.where(var > 0.0,1.0,var)
                #var1 = var2nan(np.array(var),vmin)
                c1 = ax2.contour(ds['lon'],ds['lat'],var1,1,label='P90',colors='grey',linewidths=0.5)
                # P95
                var = ds['ssc_P95'][1,:,:]
                var2 = np.where(var > 0.0,1.0,var)
                #var2 = var2nan(np.array(var),vmin)
                c2 = ax2.contour(ds['lon'],ds['lat'],var2,1,label='P95',colors='black',linewidths=0.5)
                # P95
                var = ds['ssc_P99'][1,:,:]
                var3 = np.where(var > 0.0,1.0,var)
                #var3 = var2nan(np.array(var),vmin)
                c3 = ax2.contour(ds['lon'],ds['lat'],var3,1,label='P99',colors='red',linewidths=0.5)
                #
                ax2.plot(ds.attrs['center_point'][0], ds.attrs['center_point'][1],marker='.',markersize=1.0,color='w')
                # circle / radius around point
                draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),250/111000,color='k',fill=False,linestyle='--')
                ax2.add_artist(draw_circle)
                draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),500/111000,color='k',fill=False,linestyle='--')
                ax2.add_artist(draw_circle)
                draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),1000/111000,color='k',fill=False,linestyle='--')
                ax2.add_artist(draw_circle)
                #m2.readshapefile(dredge_area, 'dredgeA')
                plotCommon(ax2)
                # plot 3
                #var3 = var2nan(var3,vmin)
                var = ds['ssc_P90'][2,:,:]
                var1 = np.where(var > 0.0,1.0,var)
                #var1 = var2nan(np.array(var),vmin)
                c1 = ax3.contour(ds['lon'],ds['lat'],var1,1,label='P90',colors='grey',linewidths=0.5)
                # P95
                var = ds['ssc_P95'][2,:,:]
                var2 = np.where(var > 0.0,1.0,var)
                #var2 = var2nan(np.array(var),vmin)
                c2 = ax3.contour(ds['lon'],ds['lat'],var2,1,label='P95',colors='black',linewidths=0.5)
                # P95
                var = ds['ssc_P99'][2,:,:]
                var3 = np.where(var > 0.0,1.0,var)
                #var3 = var2nan(np.array(var),vmin)
                c3 = ax3.contour(ds['lon'],ds['lat'],var3,1,label='P99',colors='red',linewidths=0.5)
                ax3.plot(ds.attrs['center_point'][0], ds.attrs['center_point'][1],marker='.',markersize=1.0,color='w')
                # circle / radius around point
                draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),250/111000,color='k',fill=False,linestyle='--')
                ax3.add_artist(draw_circle)
                draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),500/111000,color='k',fill=False,linestyle='--')
                ax3.add_artist(draw_circle)
                draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),1000/111000,color='k',fill=False,linestyle='--')
                ax3.add_artist(draw_circle)
                #m3.readshapefile(dredge_area, 'dredgeA')
                plotCommon(ax3)   

            elif method == 'DEP':
                
                ax1.set_ylabel('Latitude $^{\circ}$S',fontsize = 14)
                if count == len(flist)-1:
                    ax1.set_xlabel('Longitude $^{\circ}$E',fontsize = 14)
                #
                ax1.set_xlim(np.min(ds['lon']),np.max(ds['lon']))
                ax1.set_ylim(np.min(ds['lat']), np.max(ds['lat']))
                ax1.set_aspect('equal')
                #
                print(ds.attrs['center_point'])
                rad = 160/111000
                X,Y = np.meshgrid(ds['lon'],ds['lat'])
                coords=np.array((X.flatten(),Y.flatten())).T
                ci=mpl.patches.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),rad)
                cipath = mpl.path.Path(ci.get_verts())
                validcoords=cipath.contains_points(coords)
                n_in = validcoords.sum()
                #import pdb; pdb.set_trace()
                if stat == 'perc':
                    #import pdb; pdb.set_trace()
                    mask = validcoords.reshape(np.squeeze(ds['ssc_tot_settled_P90']).shape)
                    mask=mask*(bulk_vol/(np.pi*(160.**2)))#(n_in*(ds.attrs['pixelsize_m']**2))) # vol / circle area in m (pi r2)
                    print(mask.max(),mask.min())
                    print(vmin)
                    #var = ds['ssc_tot_settled_P90'][0,:,:]
                    var =1.5*var2nan(np.array(ds['ssc_tot_settled_P90'][0,:,:]+mask),0.0) # Bulking coeff
                    var1 = np.where(var > 0.0001,1.0,var)
                    #c = ax1.pcolormesh(ds['lon'],ds['lat'],var1,vmin = 0.0000005,cmap=cmo.matter)
                    c1 = ax1.contour(ds['lon'],ds['lat'],var1,1,label='P90',colors='grey',linewidths=0.5,zorder=0.5)
                    #
                    mask = validcoords.reshape(np.squeeze(ds['ssc_tot_settled_P95']).shape)
                    mask=mask*(bulk_vol/(np.pi*(160.**2)))#(n_in*(ds.attrs['pixelsize_m']**2))) # vol / circle area in m (pi r2)
                    print(mask.max(),mask.min())
                    print(vmin)
                    #var = ds['ssc_tot_settled_P95'][0,:,:]
                    var =1.5*var2nan(np.array(ds['ssc_tot_settled_P95'][0,:,:]+mask),0.0) # Bulking coeff
                    var2 = np.where(var > 0.0001,1.0,var)                    
                    c2 = ax1.contour(ds['lon'],ds['lat'],var2,1,label='P95',colors='black',linewidths=0.5,zorder=0.1)
                    #
                    mask = validcoords.reshape(np.squeeze(ds['ssc_tot_settled_P99']).shape)
                    mask=mask*(bulk_vol/(np.pi*(160.**2)))#(n_in*(ds.attrs['pixelsize_m']**2))) # vol / circle area in m (pi r2)
                    print(mask.max(),mask.min())
                    print(vmin)
                    #var = ds['ssc_tot_settled_P99'][0,:,:]
                    var =1.5*var2nan(np.array(ds['ssc_tot_settled_P99'][0,:,:]+mask),0.0) # Bulking coeff
                    var3 = np.where(var > 0.0001,1.0,var)                    
                    c3 = ax1.contour(ds['lon'],ds['lat'],var3,1,label='P99',colors='red',linewidths=0.5)
                    #
                    #ax1.legend([c1.legend_elements()[0],c2.legend_elements()[0],c3.legend_elements()[0]],['P90','P95','P99'],loc='best')
                # P95
                if stat == 'mean':    
                    # inflate all locations within a radious of the centre point
                    print(ds.attrs['center_point'])
                    rad = 160/111000
                    X,Y = np.meshgrid(ds['lon'],ds['lat'])
                    coords=np.array((X.flatten(),Y.flatten())).T
                    ci=mpl.patches.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),rad)
                    cipath = mpl.path.Path(ci.get_verts())
                    validcoords=cipath.contains_points(coords)
                    n_in = validcoords.sum()
                    mask = validcoords.reshape(np.squeeze(ds['ssc_tot_settled_'+stat]).shape)
                    mask=mask*(bulk_vol/(np.pi*(160.**2)))#(n_in*(ds.attrs['pixelsize_m']**2))) # vol / circle area in m (pi r2)
                    print(mask.max(),mask.min())
                    print(vmin)
                    var =1.5*var2nan(np.array(ds['ssc_tot_settled_'+stat][0,:,:]+mask),vmin) # Bulking coeff
                    #var =1.5*var2nan(var,vmin) # Bulking coeff
                    c = ax1.pcolormesh(ds['lon'],ds['lat'],var,vmin=0.0,vmax=vmax,cmap=cmo.matter)
                    c.set_clim(vmin,vmax)
                    if len(flist) == 5:
                        if count == 0:
                            cax = fig.add_axes([0.7, 0.74, 0.02, 0.144])
                        elif count == 1:
                            cax = fig.add_axes([0.7, 0.58, 0.02, 0.144])
                        elif count ==2:
                            cax = fig.add_axes([0.7, 0.42, 0.02, 0.144])
                        elif count == 3:
                            cax = fig.add_axes([0.7, 0.26, 0.02, 0.144])
                        elif count == 4:
                            cax = fig.add_axes([0.7, 0.10, 0.02, 0.144])

                    cbar = fig.colorbar(c, cax=cax)
                    cbar.ax.set_ylabel('Deposition Thickness [m]\n(bulking ratio %s) ' %(str(1.5)),size=10) 
                    #map.readshapefile(dredge_area, 'dredgeA')

                    # plot 10 cm and 1 mm contours
                    #CS = ax.contour(ds['lon'],ds['lat'],var,levels=[0.005,0.1],colors=['grey','black'])
                    # show area

                draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),250/111000,color='k',fill=False,linestyle='--')
                ax1.add_artist(draw_circle)
                draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),500/111000,color='k',fill=False,linestyle='--')
                ax1.add_artist(draw_circle)
                draw_circle=plt.Circle((ds.attrs['center_point'][0],ds.attrs['center_point'][1]),1000/111000,color='k',fill=False,linestyle='--')
                ax1.add_artist(draw_circle)
                ax1.plot(ds.attrs['center_point'][0], ds.attrs['center_point'][1], marker='.', markersize=1.0,color='w')

                plotCommon(ax1) 
            
            count +=1

        # display
        print('Saving figure ')
        plt.savefig(outname,bbox_inches = 'tight',
        pad_inches = 0.1)

def initialize_combined_dataset(sed_class=None,
                                dredge=None,
                                area=1,
                                sitename='nelson',
                                pixelsize = 20,
                                fdir=None):
    # The function initializes the xarray dataset to be used to save "combined" pdfs including contributions 
    # from all different sediment classes & # release dates
    # 
    # This is done by looping through all different processed discrete files (one per release, per class)
    # and taking the one with the longest time vector, all values are then set to zero
    # We force the dataset to include both pdf_settled, and pdf_stranded

    # Find the longest simulation, and use that as template file
    print('Initialize combined dredging/disposal_dataset')
    ntime = 0
    for i_class,sed_class in enumerate(sed_class):
            fname = np.sort(glob.glob(fdir+'/processed_pdfs/opendrift_%s_%s_%s_%s_*_processed_dx%i.nc' %(sitename,dredge,sed_class,area,pixelsize)))
            for i_file, file in enumerate(np.sort(fname)): 
                print(file)
                ds = xarray.open_dataset(file)
                if ds['time'].shape[0] > ntime:
                    fname_longest = file
                    ntime = ds['time'].shape[0]

    ds_tmp = xarray.open_dataset(fname_longest)
    # hard copy of the dataset as template for combined file
    ds_combined = ds_tmp.copy(deep = True)
    # clean up
    del(ds_tmp)
    del(ds)
    gc.collect()
    # set all to variables to 0.0
    for var in ds_combined.data_vars :
        
        if 'pdf' in var: # pdf fields
            ds_combined[var][:,:,:] = 0.0
            if 'active' in var:
                # update attributes
                ds_combined[var].attrs = {'units': '[kg/m3]',
                                          'standard_name' : var + '_TSS_combined'}
            elif 'settled' in var:
                ds_combined[var].attrs = {'units': '[m]',
                                          'standard_name' : var + '_thickness_combined_cumulative'}
            elif 'stranded' in var:
                ds_combined[var].attrs = {'units': '[m]',
                                          'standard_name' : var + '_thickness_combined_cumulative'}
            else:
                import pdb;pdb.set_trace()

        elif 'nb_part' : # number of particles
            ds_combined[var][:] = 0.0
            if 'active' in var:
                # update attributes
                ds_combined[var].attrs = {'units': '[kg]',
                                          'standard_name' : var + '_mass_combined'}
            else: # for settled and stranded particles, we track total volume rather than mass
                ds_combined[var].attrs = {'units': 'm3',
                                          'standard_name' : var + '_volume_combined'}

    # Make sure the dataset includes both settled and stranded pdfs, 
    # (all other pdf relating to 'active' particles will always be in all files, so already included)
    if 'pdf_stranded' not in ds_combined.data_vars :
        ds_combined['pdf_stranded'] = (['time','y','x'],  np.zeros( (ds_combined['time'].shape[0],ds_combined['lat'].shape[0],ds_combined['lon'].shape[0]) ) )
        ds_combined['pdf_stranded'].attrs = {'units': '[m]',
                                             'standard_name' : 'pdf_stranded' + '_thickness_combined_cumulative'}

        ds_combined['nb_part_stranded'] = (['time'],  np.zeros( ds_combined['time'].shape[0]) )
        ds_combined['nb_part_stranded'].attrs = {'units': '[m3]',
                                                'standard_name' : var + '_volume_combined'}
    if 'pdf_settled' not in ds_combined.data_vars :
        ds_combined['pdf_settled'] = (['time','y','x'],  np.zeros( (ds_combined['time'].shape[0],ds_combined['lat'].shape[0],ds_combined['lon'].shape[0]) ) )
        ds_combined['pdf_settled'].attrs = {'units': '[m]',
                                            'standard_name' : 'pdf_settled' + '_thickness_combined_cumulative'}

        ds_combined['nb_part_settled'] = (['time'],  np.zeros( ds_combined['time'].shape[0]) )
        ds_combined['nb_part_settled'].attrs = {'units': '[m3]',
                                                'standard_name' : var + '_volume_combined'}

    return ds_combined

def populate_dataarray_from_pdfs(xr_dataset= None, 
                                sed_type=None,
                                area=1,
                                site_name = 'nelson',
                                pixelsize = 20,
                                fdir = None,
                                dredge = None,
                                dredge_area = 1,
                                mass_sand = None,
                                mass_silt = None,
                                mass_clay = None,
                                mass_gravel = None):
    
    print('Populate combined dredging / drillcutting_dataset')

    ################################################################################################
    # Loop through releases and sediment classes, and add contributions to total TSS and deposition
    ################################################################################################      
    
    for i_sed,sed_class in enumerate(sed_type):
            print(i_sed,sed_class)
            fname = glob.glob(fdir+'/processed_pdfs/opendrift_%s_%s_%s_%s_*_processed_dx%i.nc' %(site_name,dredge,sed_class,area,pixelsize))
            for i_file, file in enumerate(np.sort(fname)): 
                print(file)
                ds = xarray.open_dataset(file)
                ntime = ds['time'].shape[0]
                # work out particle load (i.e. mass/particle) for that discrete release/sed class
                nb_parts = ds.attrs['nb_part_total'] # total number of particle release through simulation
                mass_per_part = eval('mass_'+sed_class)[dredge_area-1] / nb_parts
                if sed_class == 'sand':
                    sed_vol = eval('mass_'+sed_class)[dredge_area-1] / 1600.0
                elif sed_class == 'silt' or 'clay':
                    sed_vol = eval('mass_'+sed_class)[dredge_area-1] / 500.0

                vol_per_part = sed_vol / nb_parts
                print('Mass: ', mass_per_part)
                print('Volume: ', vol_per_part)

                for var in ds.data_vars :
                    print('Adding contribution of %s - %s to variable : %s' % (sed_class,area,var) )
                    if 'pdf' in var: # pdf fields
                        if 'active' in var:
                            # add contribution for the correct time period 
                            # [kg/cell] 
                            xr_dataset[var][0:ntime,:,:] = xr_dataset[var][0:ntime,:,:] + (ds[var].copy(deep=True)* mass_per_part)
                        else:
                            # For cumulative variables such as settled/stranded/died etc...we need to expand the datacube ds[var] 
                            # to xr_dataset[var] shape. The new array is backfilled using the last good values before deactivation 
                            if area == 1: # assume it takes 
                                days = 11.27
                            elif area == 2:
                                days = 5.6
                            elif area == 3:
                                days = 13.405

                            dt = (ds['time'][1]-ds['time'][0]).dt.seconds
                            last_i = np.round(days*dt)
                            # gets the last timestep for the appropriate lenght of time it will take to dredge the area specified
                            last_step = ds[var][last_i,:,:].copy(deep=True)
                            # tile the last matrix in time dimension to backfill array with last value before activation
                            last_step_tiled = np.tile(last_step,[xr_dataset[var].shape[0]-ds[var].shape[0],1,1]) 
                            backfilled_datacube = np.append(ds[var].copy(deep=True),last_step_tiled, axis = 0) # axis 0 is time dimension
                            if 'settled' in var or 'stranded' in var:
                                # [m3/cell] 
                                xr_dataset[var] = xr_dataset[var] + (backfilled_datacube * volume_per_part)
                            else: # any other case # use mass by default for now
                                xr_dataset[var] = xr_dataset[var] + (backfilled_datacube * mass_per_part)

                    elif 'nb_part' in var : # 
                        if 'active' in var:
                            # convert nb of particle to total mass
                            xr_dataset[var][0:ntime] = xr_dataset[var][0:ntime] + (ds[var].copy(deep=True) * mass_per_part)
                        else:
                            # For cumulative variables such as settled/stranded/died etc...we need to expand the vector or particle
                            # number ds[var]. The new array is backfilled using the last good values before deactivation 
                            if area == 1: # assume it takes 
                                days = 11.27
                            elif area == 2:
                                days = 5.6
                            elif area == 3:
                                days = 13.405

                            dt = (ds['time'][1]-ds['time'][0]).dt.seconds
                            last_i = np.round(days*dt)
                            last_step = ds[var][last_i].copy(deep=True)
                            last_step_tiled = np.tile(last_step,[xr_dataset[var].shape[0]-ds[var].shape[0]]) 
                            backfilled_vector = np.append(ds[var].copy(deep=True),last_step_tiled, axis = 0) # axis 0 is time dimension
                            xr_dataset[var] = xr_dataset[var] + (backfilled_vector * mass_per_part)
                # rename so that all different locations modelled can be contained within one file
                xr_dataset[var+'_'+'location'+'_'+str(area)]=xr_dataset[var]
                #xr_dataset = xr_dataset.drop([var])
                print(var, ' to ')
                var = var+'_'+'location'+'_'+str(area)
                print(var)

                del(ds)
                gc.collect()                  

def distance(s_lat, s_lng, e_lat, e_lng):
    #import pdb;pdb.set_trace()
    # https://gist.github.com/rochacbruno/2883505
    # approximate radius of earth in km
    R = 6373.0
    
    s_lat = s_lat*np.pi/180.0                      
    s_lng = np.deg2rad(s_lng)     
    e_lat = np.deg2rad(e_lat)                       
    e_lng = np.deg2rad(e_lng)  
    
    d = np.sin((e_lat - s_lat)/2)**2 + np.cos(s_lat)*np.cos(e_lat) * np.sin((e_lng - s_lng)/2)**2
    
    return 2 * R * np.arcsin(np.sqrt(d)) 

def stats_dataset(fdir=None,
                  stats=None,
                  method=None,
                  file_list = None,
                  output_filename = None,
                  coords=None):

    # thresholds and options - hardcoded for now
    deposition_thick_threshold_mm = [0.0025,0.05,1] # in mm
    spreading_radius = [160] #radius [m] to apply thickness of sediment

    if stats is None or method == 'DEP':
        ssc_stats_to_process = ['mean','max','P90','P95','P99']
    else:
        ssc_stats_to_process = stats
    
    ssc_threshold_mg_per_liters = [10,50,100] # in mg/L    

    # create a xarray dataset where stats will be saved
    ds = xarray.open_dataset(file_list)
    ds_stats = xarray.Dataset(coords={'lon': (['x'], ds['lon'].data),
                    'lat': (['y'], ds['lat'].data),
                    'lon_corner': (['x_corner'], ds['lon_corner'].data), # saved for plotting purposes only
                    'lat_corner': (['y_corner'], ds['lat_corner'].data), # saved for plotting purposes only
                    'time': range(0,len(ds.data_vars))})  # time is actually event number here
    ds_stats.attrs = ds.attrs # save attributes
    ds_stats.attrs['center_point'] = coords
    ds_stats.attrs['file_list'] = fdir
    # make up a grid for distance computations
    xx,yy = np.meshgrid(ds_stats['lon'].data,ds_stats['lat'].data)
    xx = xx.ravel()
    yy = yy.ravel()
    for stat_param in ssc_stats_to_process:
        for var in ds.data_vars:
            if 'ssc' in var:
                if stat_param != 'exceed':
                    ds_stats['ssc_' + stat_param] = (['time','y','x'],  np.zeros( (len(ds.data_vars),ds['lat'].shape[0],ds['lon'].shape[0]) ) )
                    ds_stats['ssc_' + stat_param].attrs = {'units' : 'mg/L'}
                else: # make variable for each exceedence threshold
                    for thresh in ssc_threshold_mg_per_liters:
                        ds_stats['ssc_' + stat_param + '_'+str(thresh)] = (['time','y','x'],  np.zeros( (len(ds.data_vars),ds['lat'].shape[0],ds['lon'].shape[0]) ) )
                        ds_stats['ssc_' + stat_param + '_'+str(thresh)].attrs = {'units' : '%'}
            if 'settled' in var:
                ds_stats[var+'_'+stat_param] = (['time','y','x'],  np.zeros( (len(ds.data_vars),ds['lat'].shape[0],ds['lon'].shape[0]) ) )
                ds_stats[var+'_'+stat_param].attrs = ds[var].attrs

    # loop through all files, and then variables to process stats
    for i,var in enumerate(ds.data_vars):
        for stat_param in ssc_stats_to_process:            
            if 'ssc' in var: #convert to mg/L
                print(var + '_' + stat_param)
                if stat_param == 'mean':
                    ds_stats['ssc' + '_' + stat_param][i,:,:] = 1e3*ds[var].mean(dim='time')
                elif stat_param == 'median':
                    ds_stats['ssc' + '_' + stat_param][i,:,:] = 1e3*ds[var].median(dim='time')
                elif stat_param == 'P90':
                    ds_stats['ssc' + '_' + stat_param][i,:,:] = 1e3*ds[var].quantile(q=0.90,dim = 'time',interpolation = 'higher')
                elif stat_param == 'P95':
                    ds_stats['ssc' + '_' + stat_param][i,:,:] = 1e3*ds[var].quantile(q=0.95,dim = 'time',interpolation = 'higher')
                elif stat_param == 'P99':
                    ds_stats['ssc' + '_' + stat_param][i,:,:] = 1e3*ds[var].quantile(q=0.99,dim = 'time',interpolation = 'higher')                                    
                elif stat_param == 'max':
                    ds_stats['ssc' + '_' + stat_param][i,:,:] = 1e3*ds[var].max(dim='time')
                elif stat_param == 'exceed':
                    for thresh in ssc_threshold_mg_per_liters :
                        thresh_str = str(thresh)
                        dt_hours=float(ds['time'][1].data-ds['time'][0])/1e9/3600
                        len_time = ds[var].shape[0]
                        ds_stats['ssc' + '_' + stat_param + '_' + thresh_str][i,:,:] = 1e2*((1e3*ds[var]>thresh).sum(dim='time')/len_time)#*dt_hours # sum of hours > thresh = % exceed

            if 'settled' in var: # stranded or settled, in that case want the last time step
                #import pdb; pdb.set_trace()
                if stat_param == 'mean':
                    ds_stats[var + '_' + stat_param][i,:,:] = ds[var].mean(dim='time')
                elif stat_param == 'max':
                    ds_stats[var + '_'+stat_param][i,:,:] = ds[var].max(dim='time')
                elif stat_param == 'P90':
                    ds_stats[var + '_' + stat_param][i,:,:] = ds[var].quantile(q=0.90,dim='time',interpolation = 'higher')     
                elif stat_param == 'P95':
                    ds_stats[var + '_' + stat_param][i,:,:] = ds[var].quantile(q=0.95,dim='time',interpolation = 'higher')
                elif stat_param == 'P99':
                    ds_stats[var + '_' + stat_param][i,:,:] = ds[var].quantile(q=0.99,dim='time',interpolation = 'higher')    
                                       
    #import pdb; pdb.set_trace()

    if output_filename is not None:
        ds_stats.to_netcdf(path = output_filename)

class OpenDriftPostProcess(object):
    """Wraps all post-processing operations on Opendrift netcdf output files.

    Arguments:
    opendrift_output_file: Opendrift output netcdf file
    """
    def __init__(self,
                 opendrift_output_file = None,
                 **kwargs):
        
        self.opendrift_output_file = opendrift_output_file
        # define some details for export
        path, fname = os.path.split(self.opendrift_output_file)
        self.opendrift_output_file_path = path
        self.opendrift_output_file_fname = fname
        self.processed_path = os.path.join(self.opendrift_output_file_path,'processed_pdfs')
    
    def load_opendrift_object(self):    
        # load the output file as a native opendrift object.
        # This allows accessing all built-in functions already available
        self.opendrift_object = opendrift.open(self.opendrift_output_file)
        
    def load_xarray_dataset(self):
        # load as xarray dataset for future operations
        self.opendrift_xr = xarray.open_dataset(self.opendrift_output_file)

        super(OpenDriftPostProcess, self).__init__()

    def plot_simple(self,**kwargs):
        import matplotlib.pyplot as plt
        lon = self.opendrift_object.history['lon']
        lat = self.opendrift_object.history['lat']
        lon360 = self.to_longitude_0_360(lon)
        plt.plot(lon360, lat,'k.')
        # adding a frame to test extents
        if 'frame' in kwargs.keys():
            plt.plot([ kwargs['frame'][0],kwargs['frame'][0],kwargs['frame'][1],kwargs['frame'][1],kwargs['frame'][0] ],
            [ kwargs['frame'][2],kwargs['frame'][3],kwargs['frame'][3],kwargs['frame'][2],kwargs['frame'][2] ],'r')

        plt.ion()
        plt.show()

    def to_longitude_0_360(self,longitude):
        '''
        converts  longitude to different convention
        from [-180<longitude<180] to [0<longitude<360]
        '''
        # lon360 = lon180.deepcopy() # convert lon to [0-360] convention
        longitude[np.where(longitude<0)] = 360+ longitude[np.where(longitude<0)]
        return longitude

    def compute_density_array(self, pixelsize_m = 1000.0, 
                                    frame = None, # None or [lon1,lon2,lat1,lat2] 
                                    weight_name=None,  # use a variable to weight the pdf, expects variable name e.g. 'oil_mass''
                                    center_point = None, # None or [lon,lat] used to adjust grid so that middle point is right on the user-input location
                                    export_to_netcdf = True,
                                    pdf_options = None, 
                                     **kwargs):
        ''' adapted from Opendrift native "get_density_array" in basemodel.py
         
        Same as compute_density_array but using xarray only for particle positions and subsetting

        Parameters :

            pixelsize_m : 1000.0 grid pixel size in meters 
            frame : None,  None or [lon1,lon2,lat1,lat2] 
            weight:None,  None or use a variable to weight the pdf e.g. ['oil_mass'], 
            center_point : None, # None or [lon,lat] used to adjust grid so that middle point is right on the user-input location
            export_to_netcdf : True, # save to netcdf fo;e
            pdf_options : None, or dictionary with following items:
                        pdf_options = {'pdf_method' : 'numpy.histogram2d', or 'kde2D_sklearn', or 'histo2d_fast',
                                       'vertical_levels' : [['all'],[-5,-10],[0,-50],['mid+5','mid-5'],['seafloor+10','seafloor'],
                                       'time_frame': [datetime_start, datetime_end], # to look at one particular period      - to implement
                                       'output_timestep_hours': 6.0, # if needed to be different from native opendrift file  - to implement
                                       }

                                       # 'vertical_levels' = vertical levels to use for computing pdf, use either numeric range or keyword-based ['mid+5','mid-5']
                                       #    vertical level are expected to be input as [upper_limit,lower_limit] 
                                       #    ['all'] will consider all active particle in water column
             **kwargs : additional keyword arguments

                        normalize_pdf : True/False Normalize pdf by cell surface area and, only for active particle, vertical depth band. - Optional

        Returns:

            ds : xarray dataset with processed pdfs,
            write netcdf if export_to_netcdf is True
 
        '''
        # save function arguments for future use in submethods
        self.pixelsize_m = pixelsize_m
        self.frame = frame
        self.weight_name = weight_name
        self.center_point = center_point
        self.export_to_netcdf = export_to_netcdf
        self.pdf_options = pdf_options
        
        #check if file already exists
        self.fname_processed = self.opendrift_output_file_fname[:-3] + '_processed_dx%s' % (int(pixelsize_m)) + '.nc'
        if os.path.exists(os.path.join(self.processed_path,self.fname_processed)):
            print(os.path.join(self.processed_path,self.fname_processed) + ' already processed')# - Loading as xarray dataset')
            # self.opendrift_xr_processed = xarray.open_dataset(os.path.join(self.processed_path,self.fname_processed))
            #  >> this seems to be too much memory use when loading many files
            return

        # only load these objects if we need to process
        # self.load_opendrift_object()
        self.load_xarray_dataset()
        
        # using xarray
        lon = self.opendrift_xr['lon']
        lat = self.opendrift_xr['lat']
        times = self.opendrift_xr['time']
        z = self.opendrift_xr['z']

        if (lon<0).any():
            # convert longitude to [0-360] convention
            self.opendrift_xr['lon'].data = self.to_longitude_0_360(self.opendrift_xr['lon'].data) 
            
        lon_array,lat_array = self.define_processing_grid(pixelsize_m = pixelsize_m,frame = frame, center_point = center_point)

        bins=(lon_array, lat_array)

        if pdf_options is None:
            # build dict with default item values
            pdf_options = {'pdf_method': 'numpy.histogram2d', 
                           'vertical_levels': [ ['all'] ]}       
        # save water depth at particle position if needed
        if 'mid' in str(pdf_options['vertical_levels']) or 'seafloor' in str(pdf_options['vertical_levels']):
            need_water_depth = True
            sea_floor_depth_below_sea_level = self.opendrift_xr['sea_floor_depth_below_sea_level']
        else:
            pass

        # save variable used as weight if applicable
        if weight_name is not None:
            weight_array = self.opendrift_xr[weight_name]
            # weight_array = self.opendrift_xr[weight].data

        status = self.opendrift_xr['status']
        # subset the particle positions based on their status [active,stranded,settled,died etc...]
        self.status_categories = self.opendrift_xr['status'].attrs['flag_meanings'].split() # self.opendrift_object.status_categories 
        self.status_id = self.opendrift_xr['status'].attrs['flag_values'] 
        #status_id =[]; for i,val in enumerate(self.opendrift_object.status_categories): status_id.append(self.opendrift_object.status_categories.index(val))

        self.create_xr_dataset_processed(x =lon_array, y= lat_array,time_vector = times, pdf_options = pdf_options)
        
        ########################################################################################
        # start looping through times  and populate the xr dataset
        ########################################################################################
        for itime in range(len(self.opendrift_xr['time'])):

            print(times[itime].data)

            lon_itime =  self.opendrift_xr['lon'][:,itime].copy(deep = True)
            lat_itime =  self.opendrift_xr['lat'][:,itime].copy(deep = True)
            z_itime = self.opendrift_xr['z'][:,itime].copy(deep = True)
            status_itime = self.opendrift_xr['status'][:,itime].copy(deep = True)
            water_depth_itime = self.opendrift_xr['sea_floor_depth_below_sea_level'][:,itime].copy(deep = True)
            # sort out weight if applicable
            if weight_name is not None:
                weight_itime = self.opendrift_xr[weight_name][:,itime].copy(deep = True)
            else:
                weight_itime = None
            # loop through different status 
            for istat,stat in enumerate(self.status_categories):
                # extract positions for that time itime, and that status istat
                # a)mask data for status different than istat
                lon_tmp = np.ma.masked_where(status_itime != self.status_id[istat],lon_itime)
                lat_tmp = np.ma.masked_where(status_itime != self.status_id[istat],lat_itime)
                z_tmp = np.ma.masked_where(status_itime != self.status_id[istat],z_itime)
                water_depth_tmp = np.ma.masked_where(status_itime != self.status_id[istat],water_depth_itime)
                # b)keep only valid lon/lat for that status
                lon_tmp = lon_tmp[~lon_tmp.mask]
                lat_tmp = lat_tmp[~lat_tmp.mask]
                z_tmp = z_tmp[~z_tmp.mask]
                water_depth_tmp = water_depth_tmp[~water_depth_tmp.mask]

                # add weight if required
                if weight_itime is not None:
                    weight_tmp = np.ma.masked_where(status_itime != self.status_id[istat],weight_itime)
                    weight_tmp = weight_tmp[~weight_tmp.mask]
                else:
                    weight_tmp = None   

                if len(lon_tmp.shape) == 2:
                    # this happens when the mask is not an array and rather True or False, in that case the returned
                    # arrays lon_tmp,lat_tmp have shapes (1,nb_part_total) which doesnt work further down the track
                    # flatten array
                    lon_tmp = lon_tmp.flatten()
                    lat_tmp = lat_tmp.flatten()
                    z_tmp = z_tmp.flatten()
                    water_depth_tmp = water_depth_tmp.flatten()
                    if weight_itime is not None:  weight_tmp = weight_tmp.flatten()

                # save the total number of particles with that status, at that time
                nb_part_itime_istat = len(lon_tmp) # this is the TOTAL number of particle with that status, at that time  ***regardless of frame, vertical level etc..***
                self.opendrift_xr_processed['nb_part_' + stat][itime] = nb_part_itime_istat
                # print(stat)
                # print(nb_part_itime_istat)

                # proba density function computation 
                if stat == 'active':
                    # subset particle cloud if applicable based on pdf_options, only for 'active' particle
                    for lev in pdf_options['vertical_levels']:
                        suffix = str(lev).replace(" ", "").replace("'", "").replace(",", "_")
                        lon_sub,lat_sub,z_sub,weight_sub = self.subset_active_particles_vertical(x = lon_tmp.data,
                                                                                                 y = lat_tmp.data,
                                                                                                 z = z_tmp.data,
                                                                                                 water_depth = water_depth_tmp.data,
                                                                                                 weight = weight_tmp,
                                                                                                 vert_level = lev) 
                        # compute proba density function for that subset
                        pdf, nb_part_sub = self.compute_pdf_core(lon = lon_sub,lat = lat_sub, weights = weight_sub, bins = bins,pdf_method = pdf_options['pdf_method'])

                        self.opendrift_xr_processed['pdf_' + stat + suffix][itime,:,:] = pdf.copy() # write to xr dataset at time itime
                        del lon_sub,lat_sub,z_sub,weight_sub,nb_part_sub,pdf
                else :
                    # compute proba density function
                    pdf, nb_part = self.compute_pdf_core(lon = lon_tmp,lat = lat_tmp, weights = weight_tmp, bins = bins,pdf_method = pdf_options['pdf_method'])
                    self.opendrift_xr_processed['pdf_' + stat ][itime,:,:] = pdf.copy() # write to xr dataset at time itime
                    del pdf, nb_part,lon_tmp,lat_tmp,weight_tmp

            del  lon_itime ,lat_itime ,z_itime ,status_itime ,water_depth_itime ,nb_part_itime_istat     
        ########################################################################################
        # end of time loop
        ########################################################################################
        
        # for all pdf and particle numbers other than 'active', we choose to save the cumulative pdf
        # and number of particles rather than timeseries of each settlement/stranding etc...
        for istat,stat in enumerate(self.status_categories):
            if stat != 'active':
                self.opendrift_xr_processed['pdf_' + stat ] =  np.cumsum(self.opendrift_xr_processed['pdf_' + stat ],axis = 0)
                self.opendrift_xr_processed['nb_part_' + stat ] =  np.cumsum(self.opendrift_xr_processed['nb_part_' + stat ])
        
        # add water depth information to the dataset 
        # the depth information is retrieved from the particles' [sea_floor_depth_below_sea_level'] info
        # (only at grid cells where pdf ~=0 though).
        self.add_waterdepth_info(bins = bins)

        if 'normalize_pdf' in kwargs:
            # Normalize pdf by cell surface area and, only for active particle, vertical depth band. - Optional
            if kwargs['normalize_pdf'] :
                self.normalize_pdf_by_surface_depth()
                
        # Export to netcdf 
        if export_to_netcdf :
            if not os.path.exists(self.processed_path):
                os.mkdir(self.processed_path)
            self.opendrift_xr_processed.to_netcdf(path = os.path.join(self.processed_path,self.fname_processed))
            print('pdfs saved to %s' % (os.path.join(self.processed_path,self.fname_processed)))

    def define_processing_grid(self,pixelsize_m ,frame , center_point ):
        # work out processing grid based on input pixel size, frame, center_point
        lon = self.opendrift_xr['lon']
        lat = self.opendrift_xr['lat']
        if frame is None :
            # define pixel size in degrees based on particle cloud extent
            deltalat = pixelsize_m/111000.0  # m to degrees
            deltalon = deltalat/np.cos(np.radians((np.nanmin(lat) +
                                               np.nanmax(lat))/2)) 
            # no frame input, use full extent of particle cloud
            lat_array = np.arange(np.nanmin(lat)-deltalat,
                                  np.nanmax(lat)+deltalat, deltalat)
            lon_array = np.arange(np.nanmin(lon)-deltalon,
                                  np.nanmax(lon)+deltalon, deltalon)
            
            # save frame including the entire particle cloud
            frame = [np.nanmin(lon)-deltalon,np.nanmax(lon)+deltalon,np.nanmin(lat)-deltalat,np.nanmax(lat)+deltalat]
        else : # use the user-input frame [lon1,lon2,lat1,lat2]
            # note that the item of lat_array may not be exactly equal to lon2,lat2 due to 
            # the way np.arrange works

            # define pixel size in degrees based on user-input frame      
            deltalat = pixelsize_m/111000.0  # m to degrees
            deltalon = deltalat/np.cos(np.radians((frame[2] + frame[3])/2))

            lat_array = np.arange(frame[2],frame[3], deltalat)
            lon_array = np.arange(frame[0],frame[1], deltalon)

            # could use np.linspace rather than arrange ?
            # lat_array = np.linspace(frame[2],frame[3], 1+int((frame[3]-frame[2])/deltalat))
            # lon_array = np.linspace(frame[0],frame[1], 1+int((frame[1]-frame[0])/deltalon))

            if center_point is not None :
               # adjust grid so that middle point is right on the user-input location
               # this assumes that the frame is already more or less centered on the center_point
               # this is done by inputting  frame = [center_point[0]-1,center_point[0]+1,center_point[1]-1,center_point[1]+1]
               id_lon = int(np.floor(lon_array.shape[0]/2))
               id_lat = int(np.floor(lat_array.shape[0]/2))
               lon_middle_point = lon_array[id_lon]
               lat_middle_point = lat_array[id_lat]
               dx_lon = center_point[0] - lon_middle_point
               dx_lat = center_point[1] - lat_middle_point
               lon_array = lon_array + dx_lon
               lat_array = lat_array + dx_lat
               # lon_array[id_lon] == center_point[0]
               # lat_array[id_lat] == center_point[1]

        return lon_array,lat_array

    def create_xr_dataset_processed(self,x,y,time_vector,pdf_options):

        ########################################################################################        
        # convert to xarray dataset for easy netcdf export
        # http://xarray.pydata.org/en/stable/data-structures.html#creating-a-dataset
        ########################################################################################        
        # 
        # Grid Coordinates : note that the pdf will be of size [n-1,m-1] where n and m are length of lon_array,lat_array
        # i.e. the pdf is given at the center of the cell, not corners
        # 
        # these are the corners of each grid cell 
        self.lon_corners = x
        self.lat_corners = y
        # define cell centers (array will be of length [n-1,m-1] where n and m are length of lon_array,lat_array)
        # The PDF are defined at the center of each cell 
        self.dx = np.diff(x[0:2])
        self.dy = np.diff(y[0:2])
        self.lon_center = x[:-1] + self.dx/2
        self.lat_center = y[:-1] + self.dy/2        
        # coordinates saved to file
        lon_grid_center,lat_grid_center = np.meshgrid(self.lon_center,self.lat_center)
        # ds = xarray.Dataset(coords={'lon': (['x', 'y'], lon_grid_center),
        #                 'lat': (['x', 'y'], lat_grid_center),
        #                 'time': times}) 
        ds = xarray.Dataset(coords={'lon': (['x'], self.lon_center),
                        'lat': (['y'], self.lat_center),
                        'lon_corner': (['x_corner'], self.lon_corners), # saved for plotting purposes only
                        'lat_corner': (['y_corner'], self.lat_corners), # saved for plotting purposes only
                        'time': time_vector}) 

        # add some information about the opendrift simulation, and options used
        if self.weight_name is None: self.weight_name='no weighting'
        if self.center_point is None: self.center_point='None'

        ds.attrs = {'opendrift_output_file' : self.opendrift_output_file,
                    'pixelsize_m': self.pixelsize_m,
                    'pdf_method' : pdf_options['pdf_method'],
                    'vertical_levels' : str(pdf_options['vertical_levels']),
                    'frame':  self.frame,
                    'weight': self.weight_name,
                    'center_point': self.center_point,
                    'nb_part_total' : self.opendrift_xr.dims['trajectory'] 
                    # note we have access to all native opendrift config here self.opendrift_xr.attrs
                    } 
        # now add variables
        mat_pdf_tmp  = np.zeros((len(time_vector), len(self.lat_corners) - 1,
                                 len(self.lon_corners) - 1)) #.astype(int) 
        mat_nbpart_tmp  = np.zeros(len(time_vector)) #.astype(int) 
        # matrices to be used to initialize variables

        for istat,stat in enumerate(self.status_categories):
            if stat == 'active':
                for lev in pdf_options['vertical_levels']:
                    suffix = str(lev).replace(" ", "").replace("'", "").replace(",", "_")
                    # save pdf fields for that status and vertical level
                    ds['pdf_' + stat + suffix] = (['time','y','x'],  mat_pdf_tmp.copy())
                    # save total number of particle for that status
                    ds['pdf_'+ stat + suffix].attrs = {'units': '[nb_part_' + stat + '_per_cell]',
                                                       'standard_name' : 'pdf_' + stat }
            else : # 'settled', 'stranded', 'died' etc...
                # save pdf fields
                ds['pdf_' + stat] = (['time','y','x'],  mat_pdf_tmp.copy())
                ds['pdf_'+ stat].attrs = {'units': '[nb_part_' + stat + '_per_cell]',
                                          'standard_name' : 'pdf_' + stat }
            # save total number of particles for that status
            ds['nb_part_' + stat] = (['time'],mat_nbpart_tmp.copy())
            ds['nb_part_' + stat].attrs = {'units': '[nb_particles_' + stat + '_total]', 
                                           'standard_name' : '[nb_particles_' + stat + '_total]'}
            # # save fraction of particles in frame
            # ds['nb_part_' + stta] = (['time'],pdf_dict['nb_part_' + stat])
            # ds['nb_part_' + stat].attrs = {'units': '[nb_particles_' + stat + '_total]', 
            #                                'standard_name' : '[nb_particles_' + stat + '_total]'}

        self.opendrift_xr_processed = ds
    
    def compute_pdf_core(self,lon,lat, bins,weights = None,pdf_method = 'numpy.histogram2d'): 

        if pdf_method == 'numpy.histogram2d':
            # numpy.histogram2d returns the number of particle (or weight) in each bin
            pdf , dummy, dummy = np.histogram2d(lon, lat ,weights = weights, bins = bins) # return a matrix of size [nlon,nlat]
            pdf = pdf.T # save as [nlat,nlon] - more common/intuitive
        elif pdf_method == 'kde2D_sklearn':
            pass # see method below - make sure that pdf(x,y) are provided for the same bins as np.histogram2d i.e. check cell corners vs center 
        elif pdf_method == 'histo2d_fast':
            pass
            # see https://github.com/astrofrog/fast-histogram  
        nb_part = len(lon)

        return pdf,nb_part

    def subset_active_particles_vertical(self,x,y,z,water_depth,weight,vert_level = None):

        ''' the function subset the active particle clouds based on the user-input vert_level list
        
        ***Note*** we could consider passing the id_time and read directly from the base object 

        args:

            - x,y,z  : are masked arrays for a given status, at given time ['active','settled',etc...]
                * particle depth 'z' are negative down

            - vert_level : can be ['all'],[-5,-10],['mid+10','mid-10'],['seafloor+2','seafloor']...
                       vertical level are expected to be input as [upper_limit,lower_limit] 
        
        returns:

        subset x,y,z,weight

        '''
        if hasattr(weight,'data'): weight = weight.data
        # if hasattr(z,'data'): z = z.data

        if vert_level[0] == 'all':
            # no subsetting needed - use all active particles
            return x,y,z,weight

        if len(vert_level)>1:
            if 'mid' in str(vert_level): # use band around mid-depth
                # work out the actual vertical band consider
                # it is expected that vert_level = ['mid+depth1','mid-depth2'] e.g. ['mid+10','mid-10'],['mid+14','mid-10']
                mid_depth = -np.abs(water_depth/2)# make sure depth is negative down, as particle z
                d1 = float(vert_level[0][3:]) # upper limit
                d2 = float(vert_level[1][3:]) # lower limit
                upper_limit = mid_depth+d1
                lower_limit = mid_depth+d2 # d2 is expected to be negative
                id_part = np.where(np.logical_and(z<=upper_limit, z>=lower_limit))                
            elif 'seafloor' in str(vert_level):
                # it is expected that vert_level = ['seafloor+10','seafloor+1'], or ['seafloor+10','seafloor']
                water_depth = -np.abs(water_depth)
                d1 = float(vert_level[0][8:]) # upper limit
                if vert_level[1][8:] == '':
                    d2 = 0.0
                else :
                    d2 = float(vert_level[1][8:]) # lower limit
                upper_limit = water_depth+d1
                lower_limit = water_depth+d2 # d2 is expected to be negative
                id_part = np.where(np.logical_and(z<=upper_limit, z>=lower_limit))                 
            else: # subset using numeric depth range
                upper_limit = vert_level[0]
                lower_limit = vert_level[1]
                id_part = np.where(np.logical_and(z<=upper_limit, z>=lower_limit))

            # sort out weight if needed
            
            if weight is not None: 
                weight_sub = weight[id_part]
            else:
                weight_sub = weight

            return x[id_part],y[id_part],z[id_part],weight_sub

    def normalize_pdf_by_surface_depth(self):
        # Normalize pdf by cell surface area and, only for active particle, vertical depth band. - Optional

        # convert pdf from [kg/cell] to [kg/m2]
        cell_area = self.opendrift_xr_processed.attrs['pixelsize_m'] ** 2 # use last open dataset (same cell area for all)

        for var in self.opendrift_xr_processed.data_vars:
            if 'pdf' in var: # pdf fields
                # normalize by cell area to get [particle/m2]
                self.opendrift_xr_processed[var][:,:,:] = self.opendrift_xr_processed[var][:,:,:] / cell_area
                # update units - nb part per m2
                self.opendrift_xr_processed[var].attrs['units'] = self.opendrift_xr_processed[var].units.replace('_per_cell','_per_m2')
                
                # For active/suspended particles, we need to normalize by depth to get concentrations in [nb_part/m3]
                if 'active' in var: 
                    # scan the vertical depth band considered
                    vert_band = var[1+var.find('['):var.find(']')]
                    vert_band = vert_band.replace('_',',')
                    vert_band = vert_band.replace('seafloor','0')
                    vert_band = vert_band.replace('mid','0')
                    if vert_band != 'all' : # cover all cases except 'all'
                        # define depth_band_thickness
                        vert_band = eval(vert_band) 
                        depth_band_thickness  = vert_band[0] - vert_band[1] # it is expected that depth is given as [shallower,deeper]
                        self.opendrift_xr_processed[var].attrs['depth_band_thickness'] = depth_band_thickness # add to dataset for sanity checks
                    else:
                        if 'water_depth' in self.opendrift_xr_processed:
                            # tile water depth matrix in time-dimension to fit with xr_dataset[var][:,:,:] dimensions
                            depth_band_thickness = np.tile(self.opendrift_xr_processed['water_depth'],[self.opendrift_xr_processed['time'].shape[0],1,1])
                            self.opendrift_xr_processed[var].attrs['depth_band_thickness'] = 'using variable: water_depth_center'
                        else :
                            print('Water depth information for processing level %s' % (vert_band))
                            print('Add water depth info using method : add_waterdepth_to_dataset()  ')
                            print('e.g. ds_combined = add_waterdepth_to_dataset (xr_dataset = ds_combined,config_obj = config_obj) ')
                            return
                    # normalize by depth band thickness to get to kg/m3
                    self.opendrift_xr_processed[var][:,:,:] = self.opendrift_xr_processed[var][:,:,:] / depth_band_thickness 
                    # update units - nb part per m3
                    self.opendrift_xr_processed[var].attrs['units'] = self.opendrift_xr_processed[var].units.replace('_per_m2','_per_m3')
        print('Normalized pdfs by cell surface (and depth for active particles)')
                            
    def add_waterdepth_info(self,bins):
        # add water depth information to the processed dataset
        # the depth information is retrieved from the particles' [sea_floor_depth_below_sea_level'] info
        # (only at grid cells where pdf ~=0 though). 
        
        print('Estimating water depth at grid cells and adding to dataset')
        pdf , dummy, dummy = np.histogram2d(self.opendrift_xr['lon'][:,:].data.flatten(), self.opendrift_xr['lat'][:,:].data.flatten(), weights = None, bins = bins)
        pdf[pdf ==0] = np.nan
        pdf_depth , dummy, dummy = np.histogram2d(self.opendrift_xr['lon'][:,:].data.flatten(), self.opendrift_xr['lat'][:,:].data.flatten(), weights = self.opendrift_xr['sea_floor_depth_below_sea_level'][:,:].data.flatten(), bins = bins)
        water_depth = pdf_depth.T / pdf.T

        # now add water depth to xr_dataset
        # self.opendrift_xr_processed['water_depth_corner'] = (['y_corner','x_corner'],water_depth_corner)
        # self.opendrift_xr_processed['water_depth_corner'].attrs = {'units': 'meters', 
        #                                                  'standard_name' : 'sea_floor_depth_below_sea_level'}
        self.opendrift_xr_processed['water_depth'] = (['y','x'],water_depth)
        self.opendrift_xr_processed['water_depth'].attrs = {'units': 'meters', 
                                                            'standard_name' : 'sea_floor_depth_below_sea_level'}
        # water_depth is of size consistent with [lon,lat] i.e. consistent with all other saved pdfs
        # and of size (M-1,N-1) relative to [lon_corner,lat_corner]

        # plot using 
        # import matplotlib.pyplot as plt
        # plt.ion()
        # plt.pcolor(self.opendrift_xr_processed['lon_corner'].data,self.opendrift_xr_processed['lat_corner'].data,pdf_depth.T/pdf.T)

    def animate_density_array(self):
        ds = self.opendrift_xr_processed
        import matplotlib.pyplot as plt
        plt.ion()
        for it,time_str in enumerate(ds['time']) : 
            plt.pcolormesh(ds['lon_corner'],ds['lat_corner'],ds['pdf_active[all]'].data[it,:,:])
            plt.plot(ds.attrs['center_point'][0],ds.attrs['center_point'][1])
            plt.title(time_str.values)
            plt.pause(0.1)

        if not hasattr(self,'opendrift_xr_processed'):
            print('no processed data available - run get_density_array() first')
    
    def plot_density_array(self):
        ds = self.opendrift_xr_processed
        import matplotlib.pyplot as plt
        plt.ion()
        # mean
        fig, ax = plt.subplots()
        im = ax.pcolormesh(ds['lon_corner'],ds['lat_corner'], ds['pdf_active[all]'].mean(dim = 'time'))
        cbar = fig.colorbar(im, ax=ax)
        cbar.ax.set_ylabel(ds['pdf_active'].attrs['units'])
        ax.plot(ds.attrs['center_point'][0],ds.attrs['center_point'][1],'k.')
        plt.title('MEAN')
        # P90
        fig1, ax1 = plt.subplots()
        ax1.pcolormesh(self.lon_grid_corners,self.lat_grid_corners, ds['pdf_active'].quantile(q = 0.9, dim = 'time', interpolation = 'nearest').T)
        ax1.title('P90')

        if not hasattr(self,'opendrift_xr_processed'):
            print('no processed data available - run get_density_array() first')

#################################################################
# end of OpenDriftPostProcess class methods 
#################################################################

