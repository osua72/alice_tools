###################################################################
# Run using official opendrift docker + xarray
# command to start below:
# 
# docker run -it  --rm -v /data2:/data2 -v /data2/simon/github/opendrift:/code/ -v /net/datastor1:/net/datastor1 opendrift/opendrift:py3-v1.0.7
# cd /data2/simon/0491_SydneyHarbour_Dredging
# pip install xarray # need to post-process results
# python process_opendrift_sydney_dredging.py
# 
####################################################################
import sys,os
import glob
import netCDF4
import numpy as np
import xarray

# append also folder with opendrift_postprocess.py to access all tools
sys.path.append('/data2/alice/alice_tools/python/opendrift/postprocessing')
#
import opendrift
from opendrift_postprocess_disposal import *
##########################################################
# set file path
run_folder_list =['/data2/alice/Projects/0435_Nelson/']
# save these options to config_obj, they will be used in compute_density_array()
pixelsize_m =20.0
frame = [173.210,173.240,-41.260,-41.240]#[lon1,lon2,lat1,lat2]
center_point = None
pdf_options = { 'pdf_method' : 'numpy.histogram2d',
                'vertical_levels' : [['all'],[0,-3],['mid+1','mid-1'],['seafloor+2','seafloor']]}

year = ['Nino','Nina']

#############################################################################################
# End of Loop processing into pdfs - Next weight the PDF's and do analysis
################################################################################################  

nb_parts = 1000.0
HOC = 2 # height of cell (m)

dredge = 'TSHD' # TSHD BHD 
#fdir = './outputs_wp3_long/Nina/%s/processed_pdfs/' %(dredge)
fdir = './wp3_test/Nino/%s/' %(dredge)
# shapefiles for plotting
topo = './nz-10m-satellite-imagery-2018-2019.tif'
dredge_area = './dredge_wgs'

# dredging coordinates
coords = [ [173.2216,173.2216,173.2294,173.22,173.2144],
           [-41.2483,-41.2426,-41.2481,-41.2539,-41.2483] ]

sed_type = ['sand','silt','clay']
vol2rem = [189400.,94100.,225200.] # m3 to dredge per area
sites = [3,2,3]
# fraction of sediment per area
sed_frac_sand = [44.11/100,68.63/100,45.18/100]
sed_frac_silt = [36.70/100,11.89/100,21.80/100]
sed_frac_clay = [9.40/100,4.66/100,7.19/100]
sed_frac_gravel = [9.79/100,14.82/100,25.82/100] # deposition only
#vol_gravel = [18547,13941,58149]
# Use this to establish the amount of sediment to deposit under the dredge
# m3 / num cells in area below dredge [m2] - amount to inrease thickness by in those areas
print(vol2rem)
hop_vol = [np.array(vol2rem)*.73] #TSHD / BHD 73% accounted for in deposition
print(hop_vol)

# temporary vars for plotting code
var1,var2,var3 = [],[],[]
mass_sand,mass_silt,mass_clay=[],[],[]
# give sediment mass/volumes for mass per particle & thickness calculations
for i in range(0,len(vol2rem)):
    sand = calc_vol(2,0,dredge,nb_parts,vol2rem[i]/sites[i],sites[i],'sand',sed_frac_sand[i])
    silt = calc_vol(2,0,dredge,nb_parts,vol2rem[i]/sites[i],sites[i],'silt',sed_frac_silt[i])
    clay = calc_vol(2,0,dredge,nb_parts,vol2rem[i]/sites[i],sites[i],'clay',sed_frac_clay[i])
    # for each of the dredging areas generate a list of the mass of sediment to deposit
    mass_sand.append(sand) #kg
    mass_silt.append(silt)
    mass_clay.append(clay)
    print(mass_sand,mass_silt,mass_clay)
    # convert into volume
    vol_sand = np.array([mass_sand]) / 1600.0 # dry sediment density (kg/m3)
    vol_silt = np.array([mass_silt]) / 500.0
    vol_clay = np.array([mass_clay]) / 500.0

################################################################################################
# Now do stats
################################################################################################
method = 'SSC'
BR = 1.5 # Bulking ratio, for deposition calculation
# for deposition, do statistics depending on the area dredged
area = 3 # 2 3 

if method == 'SSC':
    statslist = ['mean','P90','P95','P99','exceed']
elif method == 'DEP':
    statslist = ['mean']
switch = 0 # change how analysis is done

if method == 'SSC':
    BR = 0
#if method == 'SSC': # DEP SSC
# for each area dredged, calculate the surface, middle and bottom ssc per deposition site.
# ssc_tot_x = [ssc,disp_area]
#import pdb;pdb.set_trace()
just_plot = 0
if just_plot != 1:
    for r in range(1,len(coords[0])+1): 
        fname_s = os.path.join(fdir,'opendrift_nelson_%s_*_%s_*.nc' % (dredge,r))
        flist = sorted(glob.glob(fname_s))
        print(fname_s)
        
        if os.path.isfile(fdir+'Combined_%s_%s_%s_%s_%s.nc'%(method,dredge,coords[0][r-1],coords[1][r-1],area)):
            print('File exists! Loading: ',fdir+'Combined_%s_%s_%s_%s_%s.nc'%(method,dredge,coords[0][r-1],coords[1][r-1],area))
        else:
            if method == 'SSC':
                # return ssc [kg/m3] for each sediment type for the disposal site
                [lon,lat,h,ssc_sand_s,ssc_sand_m,ssc_sand_b] = calc_ssc(flist,mass_sand[area-1],nb_parts,method,dredge,'sand',area)
                [lon,lat,h,ssc_silt_s,ssc_silt_m,ssc_silt_b] = calc_ssc(flist,mass_silt[area-1],nb_parts,method,dredge,'silt',area)
                [lon,lat,h,ssc_clay_s,ssc_clay_m,ssc_clay_b] = calc_ssc(flist,mass_clay[area-1],nb_parts,method,dredge,'clay',area)

                # calculate total - account for differing lengths of pdf
                ntime = np.max((ssc_sand_s.shape[0],ssc_silt_s.shape[0],ssc_clay_s.shape[0]))
                ssc_tot_s = np.zeros((ntime,lat.shape[0],lon.shape[0]))
                ssc_tot_m = np.zeros((ntime,lat.shape[0],lon.shape[0]))
                ssc_tot_b = np.zeros((ntime,lat.shape[0],lon.shape[0]))
                #
                ntime = ssc_sand_s.shape[0] #get length of variable
                ssc_tot_s[0:ntime,:,:] = ssc_tot_s[0:ntime,:,:] + ssc_sand_s
                ssc_tot_m[0:ntime,:,:] = ssc_tot_m[0:ntime,:,:] + ssc_sand_m
                ssc_tot_b[0:ntime,:,:] = ssc_tot_b[0:ntime,:,:] + ssc_sand_b
                del(ssc_sand_s,ssc_sand_b,ssc_sand_m)
                #
                ntime = ssc_silt_s.shape[0]
                ssc_tot_s[0:ntime,:,:] = ssc_tot_s[0:ntime,:,:] + ssc_silt_s
                ssc_tot_m[0:ntime,:,:] = ssc_tot_m[0:ntime,:,:] + ssc_silt_m
                ssc_tot_b[0:ntime,:,:] = ssc_tot_b[0:ntime,:,:] + ssc_silt_b
                del(ssc_silt_s,ssc_silt_m,ssc_silt_b)
                #
                ntime = ssc_clay_s.shape[0]
                ssc_tot_s[0:ntime,:,:] = ssc_tot_s[0:ntime,:,:] + ssc_clay_s
                ssc_tot_m[0:ntime,:,:] = ssc_tot_m[0:ntime,:,:] + ssc_clay_m
                ssc_tot_b[0:ntime,:,:] = ssc_tot_b[0:ntime,:,:] + ssc_clay_b
                del(ssc_clay_s,ssc_clay_m,ssc_clay_b)
                print(np.shape(ssc_tot_s))
                print(np.amax(ssc_tot_s))          
                print(np.amax(ssc_tot_m)) 
                print(np.amax(ssc_tot_b))  
                #np.savez(fdir+'%s_%s_%s_%s' %(method,dredge,r,area),pdf_active_tot_s=ssc_tot_s,pdf_active_tot_m=ssc_tot_m,pdf_active_tot_b=ssc_tot_b,lon=lon,lat=lat,h=h) 
                # make netcdf
                # create a xarray dataset where stats will be saved
                T = 0
                for file in flist:
                    print(file)
                    ds = xarray.open_dataset(file)
                    if len(ds['time']) > T:
                        T = len(ds['time'])
                        fname = file
                ntime = ssc_tot_s.shape[0]
                #import pdb; pdb.set_trace()
                ds = xarray.open_dataset(fname)
                ds_stats = xarray.Dataset(coords={'lon': (['x'], ds['lon'].data),
                    'lat': (['y'], ds['lat'].data),
                    'lon_corner': (['x_corner'], ds['lon_corner'].data), # saved for plotting purposes only
                    'lat_corner': (['y_corner'], ds['lat_corner'].data), # saved for plotting purposes only
                    'time': range(0,ntime)})  # time is actually event number here
                ds_stats.attrs = ds.attrs # save attributes
                ds_stats.attrs['center_point'] = coords[0][r-1],coords[1][r-1]
                ds_stats.attrs['file_list'] = flist
                ds_stats.attrs['vertical_levels'] = '['+ds.attrs['vertical_levels'][10:-1]+']'

                ds_stats['ssc_tot_[0_-3]'] = (['time','y','x'],  np.zeros( (ntime,ds['lat'].shape[0],ds['lon'].shape[0]) ) )
                ds_stats['ssc_tot_[0_-3]'].attrs = {'units' : 'combined kg/m3'}
                ds_stats['ssc_tot_[mid+1_mid-1]'] = (['time','y','x'],  np.zeros( (ntime,ds['lat'].shape[0],ds['lon'].shape[0]) ) )
                ds_stats['ssc_tot_[mid+1_mid-1]'].attrs = {'units' : 'combined kg/m3'}
                ds_stats['ssc_tot_[seafloor+2_seafloor]'] = (['time','y','x'],  np.zeros( (ntime,ds['lat'].shape[0],ds['lon'].shape[0]) ) )
                ds_stats['ssc_tot_[seafloor+2_seafloor]'].attrs = {'units' : 'combined kg/m3'}
                #
                ds_stats['ssc_tot_[0_-3]'][:,:,:] = ssc_tot_s
                ds_stats['ssc_tot_[mid+1_mid-1]'][:,:,:] = ssc_tot_m
                ds_stats['ssc_tot_[seafloor+2_seafloor]'][:,:,:] = ssc_tot_b
                ds_stats.to_netcdf(path = fdir+'Combined_%s_%s_%s_%s_%s.nc'%(method,dredge,coords[0][r-1],coords[1][r-1],area))

            elif method == 'DEP':
                # return cumulative thickness of particles at the last timestep in m 
                [lon,lat,h,ssc_sand_settled] = calc_ssc(flist,mass_sand[area-1],nb_parts,method,dredge,'sand',area)
                [lon,lat,h,ssc_silt_settled] = calc_ssc(flist,mass_silt[area-1],nb_parts,method,dredge,'silt',area)
                [lon,lat,h,ssc_clay_settled] = calc_ssc(flist,mass_clay[area-1],nb_parts,method,dredge,'clay',area)
                # calculate total - account for differing lengths of pdf
                ntime = np.max((ssc_sand_settled.shape[-1],ssc_silt_settled.shape[-1],ssc_clay_settled.shape[-1]))
                ssc_tot_settled = np.zeros((lat.shape[0],lon.shape[0],ntime))
                #
                ssc_tot_settled = ssc_tot_settled + ssc_sand_settled
                del(ssc_sand_settled)
                #
                ssc_tot_settled = ssc_tot_settled + ssc_silt_settled
                del(ssc_silt_settled)
                #
                ssc_tot_settled = ssc_tot_settled + ssc_clay_settled
                del(ssc_clay_settled)
                print(np.shape(ssc_tot_settled))
                print(np.amax(ssc_tot_settled))           
                #np.savez(fdir+'Combined_%s_%s_%s_%s' %(method,dredge,coords[0][r-1],coords[1][r-1]),pdf_tot_settled=ssc_tot_settled,lon=lon,lat=lat,h=h) 

                # make netcdf?
                # create a xarray dataset where stats will be saved
                ntime = ssc_tot_settled.shape[-1] # cumulative total monthly sediment accumulation [m] per month
                ds = xarray.open_dataset(flist[0])
                ds_stats = xarray.Dataset(coords={'lon': (['x'], ds['lon'].data),
                    'lat': (['y'], ds['lat'].data),
                    'lon_corner': (['x_corner'], ds['lon_corner'].data), # saved for plotting purposes only
                    'lat_corner': (['y_corner'], ds['lat_corner'].data), # saved for plotting purposes only
                    'time': range(0,ntime)})  # time is actually event number here
                ds_stats.attrs = ds.attrs # save attributes
                ds_stats.attrs['center_point'] = coords[0][r-1],coords[1][r-1]
                ds_stats.attrs['file_list'] = flist
                ds_stats.attrs['vertical_levels'] = '['+ds.attrs['vertical_levels'][10:-1]+']'

                ds_stats['ssc_tot_settled'] = (['y','x','time'],  np.zeros( (ds['lat'].shape[0],ds['lon'].shape[0],ntime) ) )
                ds_stats['ssc_tot_settled'].attrs = {'units' : 'combined m'}
                #
                ds_stats['ssc_tot_settled'][:,:,:] = ssc_tot_settled
                ds_stats.to_netcdf(path = fdir+'Combined_%s_%s_%s_%s_%s.nc'%(method,dredge,coords[0][r-1],coords[1][r-1],area))
        # do stats
        if os.path.isfile(fdir+'%s_%s_%s_%s_stats_%s.nc'%(method,dredge,coords[0][r-1],coords[1][r-1],area)) is False:
            stats_dataset(fdir='./',stats=statslist,method=method,
                          file_list=fdir+'Combined_%s_%s_%s_%s_%s.nc'%(method,dredge,coords[0][r-1],coords[1][r-1],area), # must be surface, middle, bottom
                          output_filename=os.path.join(fdir,'%s_%s_%s_%s_stats_%s.nc' % (method,dredge,coords[0][r-1],coords[1][r-1],area)),
                          coords=[coords[0][r-1],coords[1][r-1]])
        else:                      
            print('check for : ',os.path.join(fdir,'Combined_%s_%s_%s_%s_stats_%s.nc' % (method,dredge,coords[0][r-1],coords[1][r-1],area)))

else:
    # map of mean values
    fname_s = os.path.join(fdir,'Combined_%s_%s_*_%s.nc' % (method,dredge,area))
    flist = sorted(glob.glob(fname_s))
    if method == 'SSC':
        make_map(flist=flist,statslist=statslist,method=method,topo=topo,dredge_area=disp_area,bulk_vol=None)
    elif method == 'DEP':    
        print(hop_vol[area])
        make_map(flist=flist,statslist=statslist,method=method,topo=topo,dredge_area=disp_area,bulk_vol=hop_vol[area])
        
        #outname = '%s_%s_%s.png' %(method,dredge,stats)
        #if statslist[0] == 'Ptime':
        #    make_map(fdir,topo,coords,ssc_tot_s,ssc_tot_m,ssc_tot_b,lon,lat,h,method,spoil_area,outname,stats,BR)
        #else:
        #make_map_ssc_dredge(fdir,topo,coords,ssc_tot_s,ssc_tot_m,ssc_tot_b,lon,lat,h,method,dredge_area,outname,stats,BR)