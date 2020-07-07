#!/usr/bin/env python2.7

import os,sys
import netCDF4
import numpy as np
import matplotlib.tri as mtri
from matplotlib.dates import YearLocator, DayLocator, DateFormatter, AutoDateLocator
import matplotlib.pyplot as plt
import pyproj
import rasterio
from rasterio.plot import show

import matplotlib.ticker as ticker

import shapefile as shp
from mpl_toolkits.basemap import Basemap
from osgeo import gdal

import datetime as dt

import json
from scipy.interpolate import griddata
import pylab

import cmocean
import cmocean.cm as cmo

def lat_lon_proportion(plt,ax):
	#Calculate the distances along the axis
	xlim=plt.xlim()
	ylim=plt.ylim()
	x_dist = np.diff(xlim);
	y_dist = np.diff(ylim);

	#Adjust the aspect ratio
	c_adj = np.cos(np.mean(np.deg2rad(xlim)));
	dar = [1,c_adj,1];
	pbar = [x_dist[0]*c_adj/y_dist[0],1,1 ];
	ax.set_aspect(abs(c_adj))
	#ax.set_aspect(abs(x_dist[0]*c_adj/y_dist[0]))

def near2d_selfe(x, y, x0, y0, nnodes=1, nreturn='multiple'):

    """
    Find the indexes of the grid point that is
    nearest a chosen (x0, y0).
    Usage: line, col = near2d(x, y, x0, y0)
    """   
    dx = np.abs(x - x0); dx = dx / dx.max()
    dy = np.abs(y - y0); dy = dy / dy.max()
    dn = dx + dy    

    if nnodes > 1:
        line = []
        for node in range(nnodes):
            fn = np.where(dn == dn.min())[0][0]
            line.append(fn)
            dn[fn] = 9999

    else:
        fn = np.where(dn == dn.min())[0][0]
        line = [int(fn)]  

    if nreturn == 'multiple':

        if nnodes > 1: return line
        else: return line[0]

    elif nreturn == 'single':
        return line[nnodes-1]


def plotpatch(bnd):

	data=np.loadtxt(bnd)
	x=data[:,0];
	y=data[:,1];

	splits=np.bitwise_or(np.isnan(x),x>99999999999,y==1).nonzero()[0]
	splits=np.insert(splits,0,-1)
	splits=np.insert(splits,len(splits),len(splits))
	pols=[]
	for ist in range(0,len(splits)-2):
		ind=np.arange(splits[ist]+1,splits[ist+1])
		plt.fill( x[ind],y[ind],'silver')

def get_min_max(dirin,params,Istart,Iend,level,quiver):
    Zmin=np.inf
    Zmax=np.inf*-1
    for k in range(Istart,Iend+1):
        for vars in params:
            fullfile=os.path.join(dirin,'schout_'+str(k)+'.nc')
            nc=netCDF4.Dataset(fullfile)
            nt=len(nc.variables['time'])
          
            for t in range(0,nt):
                if 'two'  in nc.variables[vars].dimensions:
                    tmp=nc.variables[vars][t]
                    u=tmp[...,0]
                    v=tmp[...,1]
                    Z=np.sqrt(u**2+v**2)
                    if quiver is not None:
                        tmp=nc.variables[vars][t]
                        u=tmp[...,0]
                        v=tmp[...,1]

                else:
                    Z=nc.variables[vars][t,:]
	    
                if len(Z.shape)>1:
                        Z=Z[level,:]

                if len(u.shape)>1:
                    u=u[level,:]
                    v=u[level,:]

        if Zmin==np.inf:
            z=Z

		
        Zmin=min(Zmin,min(Z))
        Zmax=max(Zmax,max(Z))




    return z,np.round(Zmax,2),np.round(Zmin,2),u,v

def extract_ts(Istart,Iend,node,dirin):
    E=[]
    T=[]

    for k in range(Istart,Iend+1):
        fullfile=os.path.join(dirin,'schout_'+str(k)+'.nc')
        nc=netCDF4.Dataset(fullfile)
        dtime = netCDF4.num2date(nc.variables['time'][:],nc.variables['time'].units)
        T=np.hstack((T,dtime))
        E=np.hstack((E,nc.variables['elev'][:,node]))
    return E,T

def get_data(dirin,varin,params):
    for vars in params:
        fullfile=os.path.join(dirin,varin+'.nc')
        nc=netCDF4.Dataset(fullfile)
        nt=len(nc.variables['time'])
          
        tmp_mean=nc.variables[vars+'_mean'][0,:,-1]
        tmp_max=nc.variables[vars+'_max'][0,:,-1]    
        tmp_min=nc.variables[vars+'_min'][0,:,-1]            

    return tmp_mean, tmp_max, tmp_min 

def plot_geo(dirin,topo,axtmp):

    def normalize(array):
        """Normalizes numpy arrays into scale 0.0 - 1.0"""
        array_min, array_max = array.min(), array.max()
        return ((array - array_min)/(array_max - array_min))

    # load topo
    print(topo)
    src = rasterio.open(os.path.join(dirin,topo))
    red = src.read(3)
    green = src.read(2)
    blue = src.read(1)
    print(src.bounds)
    #
    # Normalize the bands
    redn = normalize(red)
    greenn = normalize(green)
    bluen = normalize(blue)

    rgb = np.dstack((redn, greenn, bluen))
    #


    return src, rgb

def extract_timeslice(dirin,varin,Istart,Iend,params,part,lev):

    all_files_tmp = [y for x in os.walk(dirin) for y in glob(os.path.join(x[0], varin+'_*.nc'))]
    all_files=[]
    print(all_files_tmp)
    for file in all_files_tmp:
        [tmp,filenum]=file.replace('.nc','').split('_')
        if int(filenum)>=Istart and int(filenum)<=Iend:
            all_files.append(file)

    print("%.f files found" % len(all_files))

    for param in params:
        for i,file in enumerate(all_files):
            print('%s => %s' % (param,file))
            ncin = netCDF4.Dataset(file,  'r')
            print(param,part,lev)
            D=np.squeeze(ncin.variables[param][part,:,:,:])

            print(np.shape(D)) #hgrid,u/v
            ncin.close()

            if i==0:
                matrix[:]=np.squeeze(D[:])	

    return matrix

def read_ts(dirin,varin,istart,iend):
    infile = glob.glob(os.path.join(dirin,varin))
    print(infile)
    fin = open(infile[0],'r').readlines()
    fin=[x.split() for x in fin]
    time,depth,elev = [],[],[]
    for f in fin:
        if f[0] == 'Year':
            HEAD=f
        else:
            time = np.append(time,dt.strptime(','.join(f[0:6]),'%Y,%m,%d,%H,%M,%S'))
            depth = np.append(depth,float(f[6]))
            elev = np.append(elev,float(f[7])) 

    # find index of maximum elevation
    nlist = np.arange(istart,iend+1,1)
    nfiles = len(nlist)
    flen = len(elev)/nfiles
    print(flen)
    ti=int(np.where(elev==max(elev))[0])
    i = 0
    tidx = ti-13
    while tidx < ti+13:
        if tidx < (i+1)*flen and tidx >= i*flen:
            idx = nlist[i]
            Tidx = tidx
            if Tidx >= 48:
                print(i)
                Tidx = Tidx-(i*48)
            print('file idx ',idx, 'time idx', Tidx)
            # run variable extraction
            if tidx == ti-13:
                matrix = extract_timeslice(dirin,varin,idx,idx,'hvel',Tidx,-1)
            else: matrix = dstack([matrix,extract_timeslice('hvel',idx,Tidx,-1)])
            tidx +=1
        else:
            i+=1
            print(tidx,i,(i+1)*flen, i*flen)
            print(ti+13)
            continue

    print(np.shape(matrix))        
    print('Done')               
    return matrix

def process(figfile,dirin,varin,Istart,Iend,params,quiver,quiver_res,quiver_scale,zmin,zmax,level,lim,plot_elev,preview,point,topo,tsfile):
    
    figdir, file = os.path.split(figfile)
    if varin != 'schout':

        fig, (ax1,ax2,ax3) = plt.subplots(nrows=1,ncols=3,figsize=(11.69,8), sharex=True, sharey=True, gridspec_kw={'wspace': 0.05})

        ax1.set_aspect('equal')
        ax2.set_aspect('equal')
        ax3.set_aspect('equal')
        #
        #tt1 = ax1.text(.5, 1.05, '', transform = ax1.transAxes, va='center',fontsize = 11)
        ax1.tick_params(labelsize=10)
        ax1.xaxis.set_major_locator(plt.MaxNLocator(3))
        #
        #tt2 = ax2.text(.5, 1.05, '', transform = ax2.transAxes, va='center',fontsize = 11)
        ax2.tick_params(labelsize=10)
        ax2.xaxis.set_major_locator(plt.MaxNLocator(3))
        #
        #tt3 = ax3.text(.5, 1.05, '', transform = ax3.transAxes, va='center',fontsize = 11)
        ax3.tick_params(labelsize=10) 
        ax3.xaxis.set_major_locator(plt.MaxNLocator(3))
        #
    
    # get the static variable
    if varin == 'schout':
        staticfile=os.path.join(dirin,varin+'_'+str(Istart)+'.nc')
        if varin == 'schout' and tsfile is not None:
            hvel = read_ts(dirin,varin,istart,iend)

    else: staticfile=os.path.join(dirin,varin+'.nc')
    
    ncs=netCDF4.Dataset(staticfile)
    X=ncs.variables['SCHISM_hgrid_node_x'][:]
    Y=ncs.variables['SCHISM_hgrid_node_y'][:]
    Z=ncs.variables['depth'][:]
    ele=ncs.variables['SCHISM_hgrid_face_nodes'][...,:3]-1
    nt=len(ncs.variables['time'])
    
    if quiver is not None:
        if lim is not None:
            Xreg, Yreg = np.meshgrid(np.linspace(lim[0],lim[1], quiver_res), np.linspace(lim[2], lim[3], quiver_res))
        else:
            Xreg, Yreg = np.meshgrid(np.linspace(min(X[:]),max(X[:]), quiver_res), np.linspace(min(Y[:]), max(Y[:]), quiver_res))
        XY=np.array((X,Y)).T

    ncs.close()
    print(lim)
    if lim is not None:
        node   = near2d_selfe(X[:],Y[:],(lim[1]-lim[0])/2, (lim[3]-lim[2])/2,  nnodes=1, nreturn='single')
    else:
        node   = near2d_selfe(X[:],Y[:],np.mean(X[:]),np.mean(Y[:]),  nnodes=1, nreturn='single')
    
    if varin == 'schout':
        elev,time=extract_ts(Istart,Iend,node,dirin)
        ZZ,ZZmax,ZZmin,U,V=get_min_max(dirin,params,Istart,Iend,level,quiver)
    else: # assume plotting stats:
        print('Getting stats data')
        ZZ,ZZmax,ZZmin=get_data(dirin,varin,params) 
            
    def init_img():
            
            if varin == 'schout':
                ZZ[ZZ>Zmax]=Zmax
                ZZ[ZZ<Zmin]=Zmin
                levels = np.linspace(Zmin, Zmax, 60)
                F=plt.tricontourf(X,Y,ele,ZZ,vmin=Zmin,vmax=Zmax,cmap=plt.cm.Spectral_r,levels=levels)
                plt.clim(Zmin,Zmax)
            else:

                if topo is not None: # can't seem to get this to work, must be missing something with imshow and rgb..
                    print('plotting chart')
#                    src,rgb = plot_geo(dirin,topo,ax1)
#                    print(np.amax(rgb))
#                    print(np.amin(rgb))
#                    m1 = show(src,ax=ax1)
#                    src,rgb = plot_geo(dirin,topo,ax2)
#                    print(np.amax(rgb))
#                    print(np.amin(rgb))                    
#                    m2 = ax2.imshow(rgb)
#                    src,rgb = plot_geo(dirin,topo,ax3)
#                    m3 = ax3.imshow(rgb)            
#                    print(np.amax(rgb))
#                   print(np.amin(rgb))

                if point is not None:    

                    p1 = ax1.plot(point[0],point[1],color='k',marker = '.')
                    p2 = ax2.plot(point[0],point[1],color='k',marker = '.')
                    p3 = ax3.plot(point[0],point[1],color='k',marker = '.')

                print('Plotting: ',params[0])
                if params[0] == 'temp':
                    vmin = 12.# 12.
                    vmax = 30.# 17.
                    ylabel = r'Sea water temperature [$^{\circ}$C]'
                    cmap = cmo.thermal
                    levels = np.linspace(vmin,vmax,vmax-vmin)
                elif params[0] == 'salt':
                    vmin = 30.
                    vmax = 40.
                    ylabel = 'Sea water salinity [psu]'
                    cmap = cmo.haline
                    levels = np.linspace(vmin,vmax,vmax-vmin)

                print(np.amax(ZZmax))
                print(np.mean(ZZ))
                print(np.amin(ZZmin))

                F1=ax1.tricontourf(X,Y,ele,ZZ,vmin=vmin,vmax=vmax,cmap=cmap,levels=levels)
                F1.set_clim(vmin, vmax)
                #
                F2=ax2.tricontourf(X,Y,ele,ZZmin,vmin=vmin,vmax=vmax,cmap=cmap,levels=levels)
                F2.set_clim(vmin, vmax)
                #
                F3=ax3.tricontourf(X,Y,ele,ZZmax,vmin=vmin,vmax=vmax,cmap=cmap,levels=levels)
                F3.set_clim(vmin, vmax)                

            
                ax1.set_xlabel('Easting (meters)',fontsize = 12)
                ax2.set_xlabel('Easting (meters)',fontsize = 12)
                ax3.set_xlabel('Easting (meters)',fontsize = 12)
                ax1.set_ylabel('Northing (meters)',fontsize = 12)

                ax1.set_title('Mean', fontsize=12)
                ax2.set_title('Min', fontsize=12)
                ax3.set_title('Max', fontsize=12)


                # cbar = fig.colorbar(F1, ax=ax1)
                # cbar.ax.set_ylabel(ylabel,size=11)
                # cbar.ax.tick_params(labelsize=11) 
                # #
                # cbar = fig.colorbar(F2, ax=ax2)
                # cbar.ax.set_ylabel(ylabel,size=11)
                # cbar.ax.tick_params(labelsize=11) 
                # cax = fig.add_axes([0.91, 0.145, 0.02, 0.70])
                
                fig.subplots_adjust(right=0.8)
                cax = fig.add_axes([0.82, 0.25, 0.02, 0.49])
                cbar = fig.colorbar(F3,cax=cax)
                cbar.ax.set_ylabel(ylabel,size=11)
                #cbar.ax.tick_params(labelsize=11) 
            
                if lim is not None:
                    ax1.set_xlim([lim[0], lim[1]])
                    ax1.set_ylim([lim[2], lim[3]])
                    ax2.set_xlim([lim[0], lim[1]])
                    ax2.set_ylim([lim[2], lim[3]])
                    ax3.set_xlim([lim[0], lim[1]])
                    ax3.set_ylim([lim[2], lim[3]])                
                else:
                    ax1.set_xlim([X.min(), X.max()])
                    ax1.set_ylim([Y.min(), Y.max()])
                    ax2.set_xlim([X.min(), X.max()])
                    ax2.set_ylim([Y.min(), Y.max()])                
                    ax3.set_xlim([X.min(), X.max()])
                    ax3.set_ylim([Y.min(), Y.max()])                  

            if varin == 'schout':
                # ADD ELEVATION
                if plot_elev:
                    rect = [0.1,0.1,0.3,0.2] # l,b,w,h
                    ax4 = fig.add_axes(rect)
                    ax4.plot(time,elev,color='b', lw=2)
                    zeros = elev*0
                    tide=ax4.plot([time[0],time[0]],[-1,1], color='k')
                    ax4.set_ylim([-1,1])
                    ax4.set_ylabel('elevation [m]',fontsize = 30)
                    ax4.xaxis.set_major_locator(   DayLocator() )
                    ax4.xaxis.set_major_formatter( DateFormatter( '%d ' ) )
                    ax4.tick_params(labelsize=25)
  
            plt.savefig( os.path.join(figdir,varin+'_'+params[0]+'.png'),bbox_inches = 'tight',
                        pad_inches = 0.1)
                                
    init_img() 
            
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(prog='plot_schism.py', usage='%(prog)s fileout dirout varin INDstart INDend params [options]')
    ## main arguments
    parser.add_argument('fileout', type=str,help='name of the output PNG file')
    parser.add_argument('dirout', type=str,help='name of the output files where the data is')
    parser.add_argument('varin', type=str,help='output file name - w/o extension')
    parser.add_argument('INDstart', type=int,help='First file to take')
    parser.add_argument('INDend', type=int,help='Last file to take')
    parser.add_argument('params', type=str,nargs='+',help='name of the parameter to plot')
    ## options
    parser.add_argument('-lim', type=float,help='Xlim,Ylim',nargs='+')
    parser.add_argument('-zmin', type=float,help='minimum value')
    parser.add_argument('-zmax', type=float,help='maximum value')
    parser.add_argument('-quiver', type=str,help='name of the quiver variable to plot')
    parser.add_argument('-quiver_res', type=int,help='Quiver resolution (default 10)',default=10)
    parser.add_argument('-quiver_scale', type=float,help='Quiver scale (default 1)',default=1)
    parser.add_argument('-level',type=int,help='level to plot for 3D data',default=1)
    parser.add_argument('-plot_elev',type=int,help='Plot elevation graph at boundary',default=True)
    parser.add_argument('-preview',type=int,help='Plot elevation graph at boundary',default=False)
    parser.add_argument('-point',type=float,help='Plot data extraction location: X Y',nargs='+')
    parser.add_argument('-chart',type=str,help='Name of tiff file to plot')
    parser.add_argument('-tsfile',type=str,help='Name of ts file to read')
    args = parser.parse_args()

    ### PRINT IT ALL
    print('output name : %s' %(args.fileout))
    print('Direcory : %s' %(args.dirout))
    print('From file #%i and #%i' %(args.INDstart,args.INDend))
    print('Do parameters : %s' %(args.params))

    process(args.fileout,args.dirout,args.varin,args.INDstart,args.INDend,args.params,args.quiver,args.quiver_res,args.quiver_scale,args.zmin,args.zmax,args.level-1,args.lim,args.plot_elev,args.preview,args.point,args.chart,args.tsfile)
