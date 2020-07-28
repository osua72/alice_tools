import sys, os
import glob
import netCDF4
import numpy as np
import argparse
import importlib
import gc
#
from datetime import datetime, timedelta
import opendrift
import logging
import xarray

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from mpl_toolkits.basemap import Basemap
from osgeo import gdal
import utils
import cmocean 
import cmocean.cm as cmo

class OpenDriftPostProcess(object):
    """Wraps all plotting operations on OpenDrift output files.
    
    Arguments:
    None. All variables generated during post processing to be parsed into routine.
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

    def load_xarray_dataset(self):
        # load as xarray dataset for future operations
        self.opendrift_xr = xarray.open_dataset(self.opendrift_output_file)

        super(OpenDriftPostProcess, self).__init__()

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
    
#########################################################################################################################
# Start of tools for plotting
#########################################################################################################################    

    def plotCommon(self,ax,tickMinutes):
            ax.xaxis.set_major_locator(ticker.MultipleLocator(tickMinutes/60.0))
            ax.yaxis.set_major_locator(ticker.MultipleLocator(tickMinutes/60.0))
            ax.xaxis.set_major_formatter(ticker.FuncFormatter(utils.formatDegreesE))
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(utils.formatDegreesS))

            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(11)
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(11)

    def var2nan(self,values,var,key=None):
        values=np.array(values)
        """Replace every 0 with 'nan' and return a copy."""
        if key is not None:
            print(key)    
            if key == 'upper':
                print('removing values above/= ', var)
                values[values>=var]=np.nan
            elif key == 'lower':
                print('removing values below/= ',var)
                values[values<=var]=np.nan
            else: 
                print('Need more information in kwargs')
                sys.exit()
        else:
            if np.round(var)==0:
                print('removing values below / = ',var)
                values[values<=var]=np.nan
            else:
                print('removing values above/= ', var)
                values[values>=var]=np.nan

        return values    

    def normalize_pdf(self,ds,nb):
        """Normalise time varying pdf with # of particles"""
        # number of particles per grid cell normalised by the number of particles in the domain at that timestep
        # to obtain the fraction of partciles in the grid cell per timestep gives value between 1 and 0.
        #pdf_norm = pdf/nb[:,None,None] # particles released:# particles per cell
        vol_per_part = nb
        for var in ds.data_vars:
            if 'pdf' in var: # pdf fields   
                if 'retired' in var: 
                    return
                if 'active' in var: 
                    print(var)
                    print('Normalising data by ', vol_per_part)
                    ds[var][:,:,:] = ds[var][:,:,:] / vol_per_part
                    # update units - nb part per m2
                    ds[var].attrs['units'] = ds[var].units.replace('_per_cell','_per_cell_norm_by_release')
                        
        print('Normalized partciles by the total number of particles released / len time')
           
        #try:
        #    if len(nb)==np.shape(pdf)[0] and len(np.shape(pdf))==3:
        #        pdf_norm = pdf/nb[:,None,None]
        #    else: print('ERROR!')
        #except:
        #    print('Normalising data by: ',nb)
        #    pdf_norm = pdf/nb
    
    
    def convert_to_m3(self,ds=None,var=None):
        # convert pdf from [kg/cell] to [kg/m2]
        cell_area = ds.attrs['pixelsize_m'] ** 2 # use last open dataset (same cell area for all)
        for var in ds.data_vars:
            if 'pdf' in var: # pdf fields
                if 'retired' in var: 
                    return
                print(var)
                ## normalize by cell area to get [particle/m2]
                ds[var][:,:,:] = ds[var][:,:,:] / cell_area
                # update units - nb part per m2
                ds[var].attrs['units'] = ds[var].units.replace('_per_cell','_per_m2')
                # For active/suspended particles, we need to normalize by depth to get concentrations in [nb_part/m3]
                if 'active' in var: 
                    print(var)
                    # scan the vertical depth band considered
                    vert_band = var[1+var.find('['):var.find(']')]
                    vert_band = vert_band.replace('_',',')
                    vert_band = vert_band.replace('seafloor','0')
                    vert_band = vert_band.replace('mid','0')
                    if vert_band != 'all' : # cover all cases except 'all'
                        # define depth_band_thickness
                        vert_band = eval(vert_band)
                        depth_band_thickness  = vert_band[0] - vert_band[1] # it is expected that depth is given as [shallower,deeper]
                        ds[var].attrs['depth_band_thickness'] = depth_band_thickness # add to dataset for sanity checks
                    else:
                        if 'water_depth' in ds.data_vars:
                            #tile water depth matrix in time-dimension to fit with xr_dataset[var][:,:,:] dimensions
                            water_depth = self.var2nan(ds['water_depth'],1.0,'lower') # remove values < 1 m
                            depth_band_thickness = np.tile(water_depth,[ds['time'].shape[0],1,1])
                            ds[var].attrs['depth_band_thickness'] = 'using variable: water_depth_center'
                        else :
                            print('Water depth information for processing level %s' % (vert_band))
                            print('Add water depth info using method : add_waterdepth_to_dataset()  ')
                            print('e.g. ds_combined = add_waterdepth_to_dataset (xr_dataset = ds_combined,config_obj = config_obj) ')
                            return

                    # normalize by depth band thickness to get to kg/m3
                    ds[var][:,:,:] = ds[var][:,:,:] / depth_band_thickness 
                    # update units - nb part per m3
                    ds[var].attrs['units'] = ds[var].units.replace('_per_m2','_per_m3')
                
        print('Normalized pdfs by cell surface (and depth for active particles)')

    def plot_topo(self,topo):
        print(topo)
        # load topo
        datafile = gdal.Open(topo)
        bnd1 = datafile.GetRasterBand(1).ReadAsArray()
        bnd2 = datafile.GetRasterBand(2).ReadAsArray()
        bnd3 = datafile.GetRasterBand(4).ReadAsArray()
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

        return img,ymin,ymax,xmin,xmax

    def get_points(self,ds=None,coords=None,vol_per_part=None):
        print('extracting data at POI')
        # find where lon | lat == coords
        #while True:
        try:
                pdf_A_points = np.zeros((len(ds['time']),len(coords[0])))
                site_depth = np.zeros((len(coords[0])))
                if ds.attrs['vertical_levels'] == "[['all']]":
                    for i in range(0,np.shape(coords)[1]):
                        print(i)
                        lonidx = int([ti for ti, x in enumerate(np.abs(ds['lon'] - coords[0][i])) if x == np.min(np.abs(ds['lon'] - coords[0][i]))][0])
                        latidx = int([ti for ti, x in enumerate(np.abs(ds['lat'] - coords[1][i])) if x == np.min(np.abs(ds['lat'] - coords[1][i]))][0])
                        # get PDF
                        #pdf_A_points[:,i] = ds['pdf_active[all]'][:,int(latidx),int(lonidx)]/(np.abs(np.diff(np.append(0,ds['nb_part_active']))))
                        pdf_A_points[:,i] = ds['pdf_active[all]'][:,int(latidx),int(lonidx)]
                        site_depth[i] = ds['water_depth'][int(latidx),int(lonidx)]
        except:
                lonidx = int([ti for ti, x in enumerate(np.abs(ds['lon'] - coords[0])) if x == np.min(np.abs(ds['lon'] - coords[0]))][0])
                latidx = int([ti for ti, x in enumerate(np.abs(ds['lat'] - coords[1])) if x == np.min(np.abs(ds['lat'] - coords[1]))][0])
                # get PDF
                pdf_A_points = ds['pdf_active[all]'][:,int(latidx),int(lonidx)]
                site_depth = ds['water_depth'][int(latidx),int(lonidx)]             
        if vol_per_part is not None:
            print('Normalising coords by ', vol_per_part)
            return (pdf_A_points/vol_per_part)
        else: return pdf_A_points

    def make_ts_plot(self,
                    coords=None,
                    fdir=None,
                    fout=None,
                    vol_per_part=None):

        ds = xarray.open_dataset(self.opendrift_output_file)
        print(ds)
         # normalise by m3
        self.convert_to_m3(ds)
        print(ds)
        self.normalize_pdf(ds,vol_per_part) #mg/m3

        print(len(ds['time']))
        pdf_a_m3 = self.get_points(ds,coords)

        print('Plotting TS!') 
 
        t = ds['time']
        tmax = t[-1]
        tmin = t[0]
        ymin = 0.0
        ymax = 3e-3
        # organise figure
        fig, ax = plt.subplots(figsize=(13,8))

        print(np.shape(t),np.shape(pdf_a_m3))
        print(np.nanmax(pdf_a_m3))
        print(np.nanmin(pdf_a_m3))

        c1=ax.plot(t,pdf_a_m3[:])
        ax.legend((c1), ('Black Reef', 'Clifton Shellfish','Ngaruroro','Short Outfall','Te Awanga','Te Awanga CR','Site1','Site2'),
        loc='upper left',ncol=4,fancybox=True)
        ax.set_ylabel('Concentration [mg/L]',fontsize = 14)
        ax.set_xlim(tmin,tmax)
        ax.set_ylim(ymin,ymax)
        ax.grid(True)

        # display
        print('Saving figure ')
        outname = '%s.png' %(fout)
        plt.savefig(fdir+outname,bbox_inches = 'tight',
        pad_inches = 0.1)
    
    def make_map(self,
                topo=None,
                fdir=None,
                fout=None,
                clabel=None,
                stats_list=None,
                tickMinutes=1.,
                Err=1e-3,
                coords=None,
                vol_per_part=None):
        print(self.__dict__)
        # load variables
        ds = xarray.open_dataset(self.opendrift_output_file)
        self.convert_to_m3(ds)
        self.normalize_pdf(ds,vol_per_part) #mg/m3
        
        
        if clabel == 'Dilution':
            log_switch = 0
        else: log_switch = 0        
        
        for stats in stats_list:

            if stats == 'mean':          
                var = ds['pdf_active[all]'].mean(dim = 'time')   
                upper = 1.e-2
                lower = 1.e-4
            elif stats == 'max':
                var=ds['pdf_active[all]'].max(dim='time')
                upper = 1.e-2
                lower = 1.e-4
            elif stats == 'p80':
                var = ds['pdf_active[all]'].quantile(q = 0.8, dim = 'time', interpolation = 'nearest')
                upper =1.e-2
                lower = 1.e-4
            elif stats == 'p90':            
                var = ds['pdf_active[all]'].quantile(q = 0.9, dim = 'time', interpolation = 'nearest')
                upper = 1.e-2
                lower = 1.e-4
            elif stats == 'p95': 
                var = ds['pdf_active[all]'].quantile(q = 0.95, dim = 'time', interpolation = 'nearest')
                upper = 1.e-2
                lower = 1.e-4
            elif stats == 'median':  
                var = ds['pdf_active[all]'].quantile(q = 0.5, dim = 'time', interpolation = 'nearest')
                upper = 1.e-2
                lower = 1.e-4            
          
            print(upper,lower)
            print(np.nanmax(var))    
            print(np.nanmin(var))    

            # remove values below lower threshold
            if lower > 0.0:
                var = self.var2nan(var,lower,'lower') # mask
            else:
                var = self.var2nan(var,0.0,'lower') # mask

            print(np.nanmax(var))
            print(np.nanmin(var))
            print(upper,lower)
            #levels = np.linspace(lower,upper)
            levels = np.arange(lower,upper,np.mean((lower,upper))/2)
            if clabel == 'Dilution':
                if log_switch == 1:
                    var = np.log10(1./var)

            print(np.nanmax(var))
            print(np.nanmin(var))
            print(1./lower, 1./upper)

            print(levels)
            # plot data
            fig,ax = plt.subplots(figsize=(11.69,8.27))  # a new figure window (landscape A4)

            if topo is not None:
                img,ymin,ymax,xmin,xmax = self.plot_topo(topo)
                print(ymin,ymax,xmin,xmax)
                #basemap
                map = Basemap(projection='cyl',llcrnrlat=ymin,urcrnrlat=ymax,
                        llcrnrlon=xmin,urcrnrlon=xmax , resolution='f', ax=ax)
                # project in the original Basemap
                map.imshow(img, origin='upper', ax=ax,alpha=0.7)
                #
            #
            ax.set_ylabel('Latitude $^{\circ}$S',fontsize = 14)
            ax.set_xlabel('Longitude $^{\circ}$E',fontsize = 14)  
            #
            ax.set_xlim(np.min(ds['lon'])+0.02,np.max(ds['lon']))
            ax.set_ylim(np.min(ds['lat']), np.max(ds['lat']))
            
            figout = '%s_%s_%s.png' %(clabel,stats,fout)

            if clabel == 'conc':
                im = ax.pcolormesh(ds['lon_corner'],ds['lat_corner'], var,
                               cmap='viridis') 
                im.set_clim([0.0, upper])
                label = 'Concentration [mg/L]'
                cbar = fig.colorbar(im, ax=ax)#,ticks=levels)
                cbar.locator = ticker.MaxNLocator(nbins=11)
                cbar.update_ticks()
                cbar.ax.set_ylabel(label,size=11)
                

            elif clabel == 'Dilution':
                #levels = np.round(np.sort(1./np.linspace(lower,upper,11,endpoint=True)))
                im = ax.pcolormesh(ds['lon_corner'],ds['lat_corner'], -var,
                               cmap='viridis_r')
                print(levels)
                im.set_clim([-upper,0.0])
                cbar = fig.colorbar(im, ax=ax)
                cbar.locator = ticker.MaxNLocator(nbins=11)
                cbar.update_ticks()
                cbar.ax.set_ylabel(clabel,size=11)
                cbar.ax.set_yticklabels(['100','111','125','142','167','200','250','333','500','1000','10000'])
        
            if log_switch == 1:                 
                im.set_clim(np.log10(1.), np.log10(Err))
                cbar.ax.set_yticklabels(10**(levels))

            #if clabel == 'conc':
            #    ax.contour(ds['lon'],ds['lat'], var,
            #                   levels=[1e-4,1e-3],colors='white')
            #else:
            #    ax.contour(ds['lon'],ds['lat'], var,
            #                   levels=[1e3],colors='white')

            # display
            print('Saving figure ')
            plt.savefig(fdir+figout,bbox_inches = 'tight',
            pad_inches = 0.1) 
            #ax.pcolormesh.remove()
