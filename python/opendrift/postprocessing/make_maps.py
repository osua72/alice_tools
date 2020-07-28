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
import cmocean.cm as cm


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