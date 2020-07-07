import sys, os
import glob
import netCDF4
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from mpl_toolkits.basemap import Basemap
from osgeo import gdal

import datetime as dt

import utils

import cmocean 
import cmocean.cm as cmo

def do_calc(cell_volume,mass_per_part,pdf):
    # do SSC calcs
    mass_per_cell = mass_per_part * pdf
    # SSC
    SSC = (mass_per_cell / cell_volume) *1000 # mg/l
    return SSC

def do_stats_1(pdf_s,pdf_m,pdf_b,stats):

    if stats == 'mean':
        print('Calculating Mean', np.shape(pdf_b), 'to ')
        if len(pdf_s[0,0,:]) <=12:
            pdf_bM = np.mean(pdf_b,axis=2)
            pdf_mM = np.mean(pdf_m,axis=2)
            pdf_sM = np.mean(pdf_s,axis=2)            
        else:
            pdf_bM = np.mean(pdf_b,axis=0)
            pdf_mM = np.mean(pdf_m,axis=0)
            pdf_sM = np.mean(pdf_s,axis=0)

        print(np.shape(pdf_bM))

        return pdf_sM,pdf_mM,pdf_bM

    elif stats == 'P95':
        print('Calculating P95', np.shape(pdf_b), 'to ')
        if len(pdf_s[0,0,:]) <=12:
            #pdf_b95 = np.percentile(pdf_b,95,axis=2)
            #pdf_m95 = np.percentile(pdf_m,95,axis=2)
            #pdf_s95 = np.percentile(pdf_s,95,axis=2)  
            pdf_b95 = np.mean(pdf_b,axis=2)
            pdf_m95 = np.mean(pdf_m,axis=2)
            pdf_s95 = np.mean(pdf_s,axis=2)
        else:
            pdf_b95 = np.percentile(pdf_b,95,axis=0)
            pdf_m95 = np.percentile(pdf_m,95,axis=0)
            pdf_s95 = np.percentile(pdf_s,95,axis=0)

        print(np.shape(pdf_b95))
            
        return pdf_s95,pdf_m95,pdf_b95


    elif stats == 'P90':
        print('Calculating P90', np.shape(pdf_b), 'to ')
        if len(pdf_s[0,0,:]) <=12:
            #pdf_b90 = np.percentile(pdf_b,90,axis=2)
            #pdf_m90 = np.percentile(pdf_m,90,axis=2)
            #pdf_s90 = np.percentile(pdf_s,90,axis=2)   
            pdf_b90 = np.mean(pdf_b,axis=2)
            pdf_m90 = np.mean(pdf_m,axis=2)
            pdf_s90 = np.mean(pdf_s,axis=2)
        else:
            pdf_b90 = np.percentile(pdf_b,90,axis=0)
            pdf_m90 = np.percentile(pdf_m,90,axis=0)
            pdf_s90 = np.percentile(pdf_s,90,axis=0)

        print(np.shape(pdf_b90))

        return pdf_s90,pdf_m90,pdf_b90

def do_stats_2(var,stats):
    
    if stats == 'mean':
        print('Calculating Mean', np.shape(var), 'to ')
        varout = np.mean(var,axis=2)
    print(np.shape(varout))    
    return varout

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

def get_vars(flist,method,stats,switch):   

    if method == 'SSC':
        count = 0
        for file in flist:
            print(file)
            nc=netCDF4.Dataset(file)
            if count == 0:
                lon=nc.variables['lon'][:]
                lat=nc.variables['lat'][:]
                print(np.amin(lon),np.amax(lon),np.amin(lat),np.amax(lat))
                h=nc.variables['water_depth'][:]
                time = nc.variables['time'][:]  
                if (len(flist) > 6 and switch == 1): # memory issue work around, not ideal
                    if stats == 'mean':    
                        pdf_b = np.mean(nc.variables['pdf_active[seafloor+2_seafloor]'][:],axis=0)
                        pdf_m = np.mean(nc.variables['pdf_active[mid+1_mid-1]'][:],axis=0)
                        pdf_s = np.mean(nc.variables['pdf_active[0_-3]'][:],axis=0)
                    elif stats == 'P95':
                        pdf_b = np.percentile(nc.variables['pdf_active[seafloor+2_seafloor]'][:],95,axis=0)
                        pdf_m = np.percentile(nc.variables['pdf_active[mid+1_mid-1]'][:],95,axis=0)
                        pdf_s = np.percentile(nc.variables['pdf_active[0_-3]'][:],95,axis=0)
                    elif stats == 'P90':
                        pdf_b = np.percentile(nc.variables['pdf_active[seafloor+2_seafloor]'][:],90,axis=0)
                        pdf_m = np.percentile(nc.variables['pdf_active[mid+1_mid-1]'][:],90,axis=0)
                        pdf_s = np.percentile(nc.variables['pdf_active[0_-3]'][:],90,axis=0)            
                else:
                    pdf_b = nc.variables['pdf_active[seafloor+2_seafloor]'][:]
                    pdf_m = nc.variables['pdf_active[mid+1_mid-1]'][:]
                    pdf_s = nc.variables['pdf_active[0_-3]'][:]
                nb_part_act = nc.variables['nb_part_active'][:]
                nc.close()
                count +=1
            else:
                h=np.mean([h,nc.variables['water_depth'][:]])
                time = np.concatenate([time,nc.variables['time'][:]])
                if (len(flist) > 6 and switch == 1): 
                    if stats == 'mean':  
                        pdf_b = np.dstack([pdf_b,np.mean(nc.variables['pdf_active[seafloor+2_seafloor]'][:],axis=0)])
                        pdf_m = np.dstack([pdf_m,np.mean(nc.variables['pdf_active[mid+1_mid-1]'][:],axis=0)])
                        pdf_s = np.dstack([pdf_s,np.mean(nc.variables['pdf_active[0_-3]'][:],axis=0)])
                    elif stats == 'P95':
                        pdf_b = np.dstack([pdf_b,np.percentile(nc.variables['pdf_active[seafloor+2_seafloor]'][:],95,axis=0)])
                        pdf_m = np.dstack([pdf_m,np.percentile(nc.variables['pdf_active[mid+1_mid-1]'][:],95,axis=0)])
                        pdf_s = np.dstack([pdf_s,np.percentile(nc.variables['pdf_active[0_-3]'][:],95,axis=0)])           
                    elif stats == 'P90':
                        pdf_b = np.dstack([pdf_b,np.percentile(nc.variables['pdf_active[seafloor+2_seafloor]'][:],90,axis=0)])
                        pdf_m = np.dstack([pdf_m,np.percentile(nc.variables['pdf_active[mid+1_mid-1]'][:],90,axis=0)])
                        pdf_s = np.dstack([pdf_s,np.percentile(nc.variables['pdf_active[0_-3]'][:],90,axis=0)])
                else:
                    pdf_b = np.concatenate([pdf_b,nc.variables['pdf_active[seafloor+2_seafloor]'][:]],axis=0)
                    pdf_m = np.concatenate([pdf_m,nc.variables['pdf_active[mid+1_mid-1]'][:]],axis=0)
                    pdf_s = np.concatenate([pdf_s,nc.variables['pdf_active[0_-3]'][:]],axis=0)
                   
                nb_part_act = np.concatenate([nb_part_act,nc.variables['nb_part_active'][:]],axis=0)
                nc.close()
        return lon,lat,h,time,pdf_s,pdf_m,pdf_b,nb_part_act
            
    elif method == 'DEP':

        count = 0
        for file in flist:
            print(file)
            nc=netCDF4.Dataset(file)
            if count == 0:
                print(count)
                lon=nc.variables['lon'][:]
                lat=nc.variables['lat'][:]
                print(np.amin(lon),np.amax(lon),np.amin(lat),np.amax(lat))
                h=nc.variables['water_depth'][:]
                time = nc.variables['time'][:]
                print(len(time)) 
                pdf_settled = nc.variables['pdf_settled'][-1,:,:]
                nb_part_set = nc.variables['nb_part_settled'][:]
                nc.close()
                count+=1
            else:            
                h=np.mean([h,nc.variables['water_depth'][:]])
                time = np.concatenate([time,nc.variables['time'][:]])
                print(len(time)) 
                pdf_settled = np.dstack([pdf_settled,nc.variables['pdf_settled'][-1,:,:]])
                nb_part_set = np.concatenate([nb_part_set,nc.variables['nb_part_settled'][:]],axis=0)
                nc.close()
        return lon,lat,h,time,pdf_settled,nb_part_set  

def calc_ssc(fdir,sed,nb_parts,method,hoc,br,dredge,stats,key,switch):

    count = 0
    fname_s = os.path.join(fdir,'opendrift_nelson_%s_%s_*.nc' % (dredge,key))
    flist = sorted(glob.glob(fname_s))
    print(flist)
    if count == 0:
        for file in flist:
            dx = int(file[len(file)-5:len(file)-3])
            cell_volume = dx * dx * hoc   
            cell_area = dx * dx        
            print('Cell Volume: ', cell_volume)
            print('Cell Area: ', cell_area)
            count +=1
            
    if method == 'SSC':
        # extract variables
        [lon,lat,h,time,pdf_s,pdf_m,pdf_b,nb_part_act]=get_vars(flist,method,stats,switch)
        # convert to ssc
        mass_per_part = sed/nb_parts # get weight per particle
        print('Mass per particle (kg): ',mass_per_part)
        # SSC
        ssc_b = do_calc(cell_volume,mass_per_part,pdf_b)
        ssc_m = do_calc(cell_volume,mass_per_part,pdf_m)
        ssc_s = do_calc(cell_volume,mass_per_part,pdf_s)
        # do stats
        print(len(time))
        if len(np.shape(pdf_b))>2:
            [ssc_s,ssc_m,ssc_b]=do_stats_1(ssc_s,ssc_m,ssc_b,stats)

        return lon,lat,h,ssc_s,ssc_m,ssc_b

    elif method == 'DEP':
        # extract variables
        [lon,lat,h,time,pdf_settled,nb_part_set]=get_vars(flist,method,stats,switch)
        # convert to ssc
        # calculate volume per particle
        vol_per_part = (sed/nb_parts ) * br
        print('Vol per particle (kg): ',vol_per_part)
        # do calcs
        vol_per_cell = vol_per_part * pdf_settled # m3 per cell
        thickness = vol_per_cell / cell_area # m
        print(len(time))
        # do stats
        if len(np.shape(pdf_settled))>2:
                thickness=do_stats_2(thickness,stats)
        return lon,lat,h,thickness

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

def make_map(topo,points,var1,var2,var3,lon,lat,h,method,dredge_area,outname,stats,br):

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

    # organise figure
    if method == 'DEP':
        var = var1
        fig = plt.figure(figsize=(11.69,8.27))  # a new figure window (landscape A4)
        ax = fig.add_subplot(1, 1, 1)  # specify (nrows, ncols, axnum)
        tickMinutes = .5
    elif method == 'SSC':     
        #fig = plt.figure(figsize=(8.27,11.69))
        fig, (ax1,ax2,ax3) = plt.subplots(3,figsize=(8.27,11.69), sharex=True,gridspec_kw={'hspace': 0.05})
        tickMinutes = 1.

    #fig = plt.figure(figsize=(13,8))  # a new figure window

    # sort coords data
    xp, yp = zip(points)
    print(xp)
    print(yp)
    # load topo
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
       
        #plot SSC
        print(np.shape(var1))
        if stats == 'mean':
            vmin = 5.0
            vmax = 500.
        elif stats == 'P95':
            vmin = 0.0
            vmax = 1000.
        elif stats == 'P90':
            vmin = 0.0
            vmax = 1000.
        #     
        # sea x and y labels
        ax1.set_ylabel('Latitude $^{\circ}$S',fontsize = 14)
        ax2.set_ylabel('Latitude $^{\circ}$S',fontsize = 14)
        ax3.set_ylabel('Latitude $^{\circ}$S',fontsize = 14)
        ax3.set_xlabel('Longitude $^{\circ}$E',fontsize = 14)  
        # set x and y limits
        ax1.set_xlim(np.min(lon)+0.02,np.max(lon))
        ax1.set_ylim(np.min(lat), np.max(lat))
        ax2.set_xlim(np.min(lon)+0.02,np.max(lon))
        ax2.set_ylim(np.min(lat), np.max(lat))
        ax3.set_xlim(np.min(lon)+0.02,np.max(lon))
        ax3.set_ylim(np.min(lat), np.max(lat))
        # start plotting
        # plot1 
        var1 =var2nan(var1,vmin)
        c = ax1.pcolormesh(lon,lat,var1,vmin=vmin,vmax=vmax,cmap=cmo.matter)
        m1.plot(xp, yp, 'x',color='w')
        c.set_clim(vmin, vmax)
        # colorbar
        cbar = m1.colorbar(c, ax=ax1)
        cbar.ax.set_ylabel('SSC [mg/L]',size=11) 
        plotCommon(ax1)
        # plot2
        var2 = var2nan(var2,vmin)
        c = ax2.pcolormesh(lon,lat,var2,vmin=vmin,vmax=vmax,cmap=cmo.matter)
        m2.plot(xp, yp, 'x',color='w')
        c.set_clim(vmin, vmax)
        # colorbar
        cbar = m2.colorbar(c, ax=ax2)
        cbar.ax.set_ylabel('SSC [mg/L]',size=11) 
        plotCommon(ax2)
        # plot 3
        var3 = var2nan(var3,vmin)
        c = ax3.pcolormesh(lon,lat,var3,vmin=vmin,vmax=vmax,cmap=cmo.matter)
        m3.plot(xp, yp, 'x',color='w')
        c.set_clim(vmin, vmax)
        #colorbar
        cbar = m3.colorbar(c, ax=ax3)
        cbar.ax.set_ylabel('SSC [mg/L]',size=11) 
        plotCommon(ax3)

    elif method == 'DEP':
        #basemap
        map = Basemap(projection='cyl',llcrnrlat=ymin,urcrnrlat=ymax,
            llcrnrlon=xmin,urcrnrlon=xmax , resolution='f', ax=ax)
        # project in the original Basemap
        x,y = map(lon_source, lat_source)
        print("shape lon and lat_source: ", lon_source.shape, lat_source.shape,bnd1.T.shape)
        map.imshow(img, origin='upper', ax=ax)
        #
        ax.set_ylabel('Latitude $^{\circ}$S',fontsize = 14)
        ax.set_xlabel('Longitude $^{\circ}$E',fontsize = 14)  
        #
        ax.set_xlim(np.min(lon)+0.02,np.max(lon))
        ax.set_ylim(np.min(lat), np.max(lat))
        #
        var =var2nan(var,0.0)
        c = ax.pcolormesh(lon,lat,var,vmin=0.001,vmax=1.0,cmap=cmo.matter)
        c.set_clim(0.001, 1)
        cax = fig.add_axes([0.91, 0.145, 0.02, 0.70])
        cbar = fig.colorbar(c, cax=cax)
        cbar.ax.set_ylabel('Deposition Thickness [m] (bulking ratio %s) ' %(str(br)),size=11) 
        map.readshapefile(dredge_area, 'dredgeA')
        # plot 10 cm and 1 mm contours
        CS = ax.contour(lon,lat,var,levels=[0.005,0.1],colors=['grey','black'])
    
        map.plot(xp, yp, 'x',color='w')

        plotCommon(ax) 

    # display
    print('Saving figure ')
    plt.savefig(outname,bbox_inches = 'tight',
    pad_inches = 0.1)
