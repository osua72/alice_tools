import sys, os
import netCDF4
from opendrift_calcs import *
import numpy as np

fdir = './outputs/wp2_nino_tshd/'
method =  'SSC' # DEP
stats = 'mean' # mean P95 P90
dredge = 'TSHD'

coords = [ [173.2531,173.2576,173.2624,173.2618,173.2633,173.2712,173.2734,173.277],
           [-41.2571,-41.2621,-41.2634,-41.2665,-41.2662,-41.2563,-41.2561,-41.2545] ]

# give sediment mass/volumes for mass per particle & thickness calculations
mass_sand = 1179887.01279352 # area 3 because other sites were normalised by this one
vol_sand = 737.429382995949 # area 3
#
mass_silt = 177939.305324063  # area 3 because other sites were normalised by this one
vol_silt = 355.878610648126 # area 3
#
mass_clay = 58690.3656624504 # area 3
vol_clay = 117.380731324901  # area 3
#
nb_parts = 1000
HOC = 2 # height of cell (m)
BR = 1.5 # Bulking ratio, for deposition calculation

topo = './nz-10m-satellite-imagery-2018-2019.tif'
dredge_area = './dredge_wgs'
var1,var2,var3 = [],[],[]
# calc ssc
if method == 'SSC':
    #if stats == 'mean':
    #    ssc_ts(fdir,mass_sand,mass_silt,mass_clay,nb_parts,HOC,dredge,coords,method,stats,0)

    [lon,lat,h,ssc_sand_s,ssc_sand_m,ssc_sand_b] = calc_ssc(fdir,mass_sand,nb_parts,method,HOC,BR,dredge,stats,'sand',1)

    [lon,lat,h,ssc_silt_s,ssc_silt_m,ssc_silt_b] = calc_ssc(fdir,mass_silt,nb_parts,method,HOC,BR,dredge,stats,'silt',1)

    [lon,lat,h,ssc_clay_s,ssc_clay_m,ssc_clay_b] = calc_ssc(fdir,mass_clay,nb_parts,method,HOC,BR,dredge,stats,'clay',1)

    ssc_tot_s = ssc_sand_s+ssc_silt_s+ssc_clay_s
    ssc_tot_m = ssc_sand_m+ssc_silt_m+ssc_clay_m
    ssc_tot_b = ssc_sand_b+ssc_silt_b+ssc_clay_b

    print(np.shape(ssc_tot_s))
    print(np.amin(ssc_tot_s),np.amax(ssc_tot_s))
    # save in case of plotting faff
    np.savez('%s_vars_%s_%s' %(method,dredge,stats),ssc_tot_s,ssc_tot_m,ssc_tot_b,lon,lat,h) 
    # map of mean values
    outname = '%s_%s_%s.png' %(method,dredge,stats)
    make_map(topo,coords,ssc_tot_s,ssc_tot_m,ssc_tot_b,lon,lat,h,method,dredge_area,outname,stats,1)

elif method == 'DEP':
    stats = 'mean'
    outname = '%s_%s.png' %(method,dredge)

    [lon,lat,h,thickness_sand] = calc_ssc(fdir,vol_sand,nb_parts,method,HOC,BR,dredge,stats,'sand',0)
    [lon,lat,h,thickness_silt] = calc_ssc(fdir,vol_silt,nb_parts,method,HOC,BR,dredge,stats,'silt',0)
    [lon,lat,h,thickness_clay] = calc_ssc(fdir,vol_clay,nb_parts,method,HOC,BR,dredge,stats,'clay',0)

    thickness_tot = thickness_sand + thickness_silt + thickness_clay
    print(np.amin(thickness_tot),np.amax(thickness_tot))
    # save in case of plotting faff
    np.savez('%s_vars_%s_%s' %(method,dredge,stats),thickness_tot,lon,lat,h)
    # make topo map of nelson port
    print('making thickness map')
    make_map(topo,coords,thickness_tot,var2,var3,lon,lat,h,method,dredge_area,outname,stats,BR)
  
