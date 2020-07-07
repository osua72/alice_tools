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
# sys.path.append('../../../toolbox_simon/python_tools')
# import download_metocean_uds
# import download_cmems
from datetime import datetime,timedelta
import opendrift
import logging
import xarray

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
# enf of OpenDriftPostProcess class methods 
#################################################################

def load_drillcutting_config(config_file = None):
    # load a drill cutting config file (used previously to run simulation) 
    sys.path.append('./run_drillcut_configs') 
    config_file = os.path.basename(config_file) # keep only filename
    # remove the .py extension if present
    if config_file[-3:] == '.py':
        config_file = config_file[:-3]
    config = importlib.import_module(config_file)
    config = config.config # keep only the config dictionary
    print('Loaded configuration file %s ...' % (config_file))
    # add filename for reference
    config['filename'] = config_file

    return config

def initialize_combined_drillcut_dataset(config_obj, folder):
    # The function initializes the xarray dataset to be used to save "combined" pdfs including contributions 
    # from all different sediment releases & classes
    # 
    # This is done by looping through all different processed discrete files (one per release, per class)
    # and taking the one with the longest time vector, all values are then set to zero
    # We force the dataset to include both pdf_settled, and pdf_stranded

    # Find the longest simulation, and use that as template file
    print('Initialize combined drillcut_dataset')
    ntime = 0
    for i_release,sed_release in enumerate(config_obj['sediment_release']['names']):
        for i_class,sed_class in enumerate(config_obj['sediment_classes']['names']): 
            fname = '/processed_pdfs/' + config_obj['site_name'] + '_' + sed_class + '_'+ sed_release + '_*_processed_dx' + str(int(config_obj['pixelsize_m'])) + '.nc' 
            out_fname = glob.glob(folder + fname)
            ds = xarray.open_dataset(out_fname[0])
            if ds['time'].shape[0] > ntime:
                fname_longest = out_fname[0]
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

def populate_combined_drillcut_dataset(xr_dataset= None, config_obj = None,folder = None):
    # populates the "combined" dataset by adding contributions from each sediment release/class found in "folder"
    # 
    # The function takes will convert the pdf in [particle/cell] into [kg/m3] for the active/suspended particles, 
    # and thickness in [m] for the settled particles. 
    #
    # These conversion take into account, for each sediment release depth and class :
    #       - actual amount of sediment released, per release, per class, 
    #       - volumic mass, 
    #       - bulking ratio 
    # 
    # Note provided thickness are the "cumulative" thickness over time (i.e. not the thickness deposited at each step)
    # 
    # 
    print('Populate combined drillcut_dataset')

    ################################################################################################
    # Loop through releases and sediment classes, and add contributions to total TSS and deposition
    ################################################################################################      
    for i_release,sed_release in enumerate(config_obj['sediment_release']['names']):
        for i_class,sed_class in enumerate(config_obj['sediment_classes']['names']):
            # load processed pdf for that one release/class
            fname = '/processed_pdfs/' + config_obj['site_name'] + '_' + sed_class + '_'+ sed_release + '_*_processed_dx' + str(int(config_obj['pixelsize_m'])) + '.nc'
            out_fname = glob.glob(folder + fname) 
            if len(out_fname)>1 :
                print('Several files for that wildcard : %s' % fname )
                import pdb;pdb.set_trace()
            # load that discrete file
            ds = xarray.open_dataset(out_fname[0])
            ntime = ds['time'].shape[0]
            # work out particle load (i.e. mass/particle) for that discrete release/sed class
            nb_part_per_run = ds.attrs['nb_part_total'] # sum of active, settled, stranded, etc... particles at end of simulation

            # determine mass/particle for that release, that sediment class
            # mass per particle(release,class) = (total_volume_relased(release) * fraction_of_sed_class(class) * volumic_mass(class) ) / total_nb_particles_released
            mass_per_part = (config_obj['sediment_release']['volume_total_m3'][i_release] * 
                            (config_obj['sediment_classes']['distribution_percent'][i_class]/100) * \
                            config_obj['sediment_classes']['wet_volumic_mass'][i_class]) / nb_part_per_run
            # equivalent volume per particle, accounting for volumic mass and bulking factor
            volume_per_part = (mass_per_part / config_obj['sediment_classes']['wet_volumic_mass'][i_class]) * \
                               config_obj['sediment_classes']['bulking_ratio'][i_class]
            
            # now loop through all different variables inside the xarray dataset (i.e. various processed pdfs),  
            # and add the contribution of that particular source(release,sediment class) to the combined dataset 
            # (for the TSS, deposition, or total mass)
            for var in ds.data_vars :
                print('Adding contribution of %s - %s to variable : %s' % (sed_class,sed_release,var) )
                if 'pdf' in var: # pdf fields
                    if 'active' in var:
                        # add contribution for the correct time period 
                        # [kg/cell] 
                        xr_dataset[var][0:ntime,:,:] = xr_dataset[var][0:ntime,:,:] + (ds[var].copy(deep=True)* mass_per_part)
                    else:
                        # For cumulative variables such as settled/stranded/died etc...we need to expand the datacube ds[var] 
                        # to xr_dataset[var] shape. The new array is backfilled using the last good values before deactivation 
                        last_step = ds[var][-1,:,:].copy(deep=True)
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
                        last_step = ds[var][-1].copy(deep=True)
                        last_step_tiled = np.tile(last_step,[xr_dataset[var].shape[0]-ds[var].shape[0]]) 
                        backfilled_vector = np.append(ds[var].copy(deep=True),last_step_tiled, axis = 0) # axis 0 is time dimension
                        xr_dataset[var] = xr_dataset[var] + (backfilled_vector * mass_per_part)
            del(ds)
            gc.collect()

    ################################################################################################
    # End of Loop through releases and sediment classes
    ################################################################################################  

    # convert pdf from [kg/cell] to [kg/m2]
    cell_area = xr_dataset.attrs['pixelsize_m'] ** 2 # use last open dataset (same cell area for all)

    for var in xr_dataset.data_vars:
        if 'pdf' in var: # pdf fields
            # normalize by cell area to get [kg/m2] for TSS, or [m3/m2]=[m] for deposition
            xr_dataset[var][:,:,:] = xr_dataset[var][:,:,:] / cell_area
            
            # For active/suspended particles, we need to normalize by depth to get concentrations in [kg/m3]
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
                    xr_dataset[var].attrs['depth_band_thickness'] = depth_band_thickness # add to dataset for sanity checks
                else:
                    if 'water_depth_center' in xr_dataset:
                        # tile water depth matrix in time-dimension to fit with xr_dataset[var][:,:,:] dimensions
                        depth_band_thickness = np.tile(xr_dataset['water_depth_center'],[xr_dataset['time'].shape[0],1,1])
                        xr_dataset[var].attrs['depth_band_thickness'] = 'using variable: water_depth_center'
                    else :
                        print('Water depth information for processing level %s' % (vert_band))
                        print('Add water depth info using method : add_waterdepth_to_dataset()  ')
                        print('e.g. ds_combined = add_waterdepth_to_dataset (xr_dataset = ds_combined,config_obj = config_obj) ')
                        return
                # normalize by depth band thickness to get to kg/m3
                xr_dataset[var][:,:,:] = xr_dataset[var][:,:,:] / depth_band_thickness       
            elif 'settled' in var or 'stranded' in var :
                # no normalization by vertical depth band required          
                pass
    return xr_dataset


def add_waterdepth_to_dataset(xr_dataset = None,config_obj = None,filename = None):  
    # the function add depth to a give xarray dataset using the source netcdf data used to run 
    # simulations. Netcdf filename with depth information is obtained from config_obj
    print('Add water depth to dataset')

    if config_obj is not None: # using config_object 
        nc_file_forcing = glob.glob(os.path.join(config_obj['ocean_forcing']['path'],config_obj['ocean_forcing']['filename']))[0]
    elif filename is not None: # when using a direct filename input # when using a direct filename input
        nc_file_forcing = filename
    # open forcing dataset
    ds = xarray.open_dataset(nc_file_forcing)
    if 'h' in ds.data_vars: # ROMS native file
        water_depth_varname = 'h'
    elif 'depth' in ds.data_vars:
        water_depth_varname = 'dep'
    else:
        print('no water depth information, or unrecognized variable name')
    # interpolate water depth to grid corners and centers
    from scipy import interpolate

    x = ds[water_depth_varname].coords['lon_rho'].data.flatten()
    y =ds[water_depth_varname].coords['lat_rho'].data.flatten()
    grid_nodes_native = np.vstack((x,y)).T
    wdep = ds[water_depth_varname].data.flatten()
    print('Building LinearNDInterpolator for water depth using variable %s from %s' % (water_depth_varname,nc_file_forcing))
    # build interpolator
    f = interpolate.LinearNDInterpolator(grid_nodes_native, wdep, fill_value=np.nan, rescale=False)
    # interpolate to grid corners
    xx,yy = np.meshgrid(xr_dataset['lon_corner'].data,xr_dataset['lat_corner'].data)
    xy_corner_processed = np.vstack((xx.flatten(),yy.flatten())).T
    water_depth_corner =f(xy_corner_processed) # interpolate to grid nodes
    water_depth_corner = np.reshape(water_depth_corner,xx.shape) # reshape to matrix
    # interpolate to grid centers
    xx,yy = np.meshgrid(xr_dataset['lon'].data,xr_dataset['lat'].data)
    xy_center_processed = np.vstack((xx.flatten(),yy.flatten())).T
    water_depth_center =f(xy_center_processed) # interpolate to grid nodes
    water_depth_center = np.reshape(water_depth_center,xx.shape) # reshape to matrix
    # now water depth to xr_dataset
    xr_dataset['water_depth_center'] = (['y','x'],water_depth_center)
    xr_dataset['water_depth_center'].attrs = ds[water_depth_varname].attrs # save native attributes
    xr_dataset['water_depth_corner'] = (['y_corner','x_corner'],water_depth_corner)
    xr_dataset['water_depth_corner'].attrs = ds[water_depth_varname].attrs # save native attributes
    # save native file from which bathymetry was interpolated
    xr_dataset.attrs['bathymetry_dataset_used_for_grid_bathy_interpolation'] = nc_file_forcing

    if False: #check_plot
        import matplotlib.pyplot as plt
        plt.ion()
        plt.figure(1)
        plt.pcolor(ds[water_depth_varname].coords['lon_rho'].data,ds[water_depth_varname].coords['lat_rho'].data,ds[water_depth_varname].data)
        plt.xlim(xr_dataset['lon_corner'].data[0],xr_dataset['lon_corner'].data[-1])
        plt.ylim(xr_dataset['lat_corner'].data[0],xr_dataset['lat_corner'].data[-1])

        plt.figure(2)
        plt.ion()
        plt.pcolor(xr_dataset['lon_corner'].data,xr_dataset['lat_corner'].data,water_depth_corner)

    return xr_dataset # udpated xarray dataset

def export_combined_drillcut_dataset(xr_dataset = None, config_obj = None ) :
    # create the folder where all processed timeseries will be saved 
    folder_combined_pdf = os.path.join(config_obj['run_folder'],'combined_pdfs_' + config_obj['filename'])
    if not os.path.exists(folder_combined_pdf):
        os.mkdir(folder_combined_pdf)    
    # work out timestamp
    timestamp = str(xr_dataset['time'].data[0]).replace("-", "").replace("T", "_").replace(":", "")[:15]
    # export combined dataset to file
    fname_combined_pdf = config_obj['filename'] + '_combined_pdfs_' + timestamp + '_dx' + str(int(xr_dataset.attrs['pixelsize_m'])) + 'm.nc'
    xr_dataset.to_netcdf(os.path.join(folder_combined_pdf,fname_combined_pdf))
    print('Combined dataset exported to %s' % (os.path.join(folder_combined_pdf,fname_combined_pdf)) )

def clean_individual_processed_files(config_obj = None,folder = None):
    # remove all intermediate files (i.e. one file per release, per class) for the considered resolution
    for i_release,sed_release in enumerate(config_obj['sediment_release']['names']):
        for i_class,sed_class in enumerate(config_obj['sediment_classes']['names']): 
            fname = '/processed_pdfs/' + config_obj['site_name'] + '_' + sed_class + '_'+ sed_release + '_*_processed_dx' + str(int(config_obj['pixelsize_m'])) + '.nc'
            out_fname = glob.glob(folder + fname)
            try:
                os.remove(out_fname[0])
            except:
                pass
    print('Removed files :  /processed_pdfs/' + config_obj['site_name'] + '_*_processed_dx' + str(int(config_obj['pixelsize_m'])) + '.nc')


def process_combined_drillcut_dataset(config_obj,folder,remove_intermediates_files = False):
    # top level function that reads a given config, and folder for a give event (which can include many individual outputs files,
    # for each sediment release/class). The function calls many sub-functions to initialize, populate, and export "combined" 
    # drill cutting dataset that includes Total Suspended Solid concentration and Deposition
    # 
    # Note it is expected that the config_obj will include the pdf_options required to call compute_density_array()
    # see example in __main__()  ...we could have them as explicit parameters ...

    # first check if that file has already been processed ?
    timestamp = folder.split('/')[-2] # folder name == timestamp string
    folder_combined_pdf = os.path.join(config_obj['run_folder'],'combined_pdfs_' + config_obj['filename'])
    fname_combined_pdf = config_obj['filename'] + '_combined_pdfs_' + timestamp + '_dx' + str(int(config_obj['pixelsize_m'])) + 'm.nc'
    if os.path.exists(os.path.join(folder_combined_pdf,fname_combined_pdf)):
        print( 'Combined file %s already exists in %s - skipping '  % (fname_combined_pdf,folder_combined_pdf ) ) 
        return

    file_list = glob.glob(os.path.join(folder,config_obj['site_name'] + '*.nc')) # files in that folder to be processed
    # post process time-varying particle density fields for each sediment release/class simulated
    for f in file_list:
        o = OpenDriftPostProcess(opendrift_output_file = f)
        ds = o.compute_density_array(pixelsize_m = config_obj['pixelsize_m'],
                                frame = config_obj['frame'],
                                center_point = config_obj['center_point'], # center grid on release point
                                export_to_netcdf = True,
                                pdf_options = config_obj['pdf_options'])
        # clear temporary variables
        del(ds)
        del(o)
        gc.collect()
    # 5) Now for each folder, combine the output files according to the config to produce the time-varying 
    # timeseries of Total Suspended Solids (TSS) and Deposition
    
    # Initialize the "combined" dataset for that simulation/folder (one folder per simulation))
    ds_combined = initialize_combined_drillcut_dataset(config_obj = config_obj,folder = folder)
    # get water depth info from external file required to get depth-averaged concentrations
    ds_combined = add_waterdepth_to_dataset(xr_dataset = ds_combined,config_obj = config_obj) 
    # populate the "combined" dataset by adding contributions from each sediemnt release/class
    ds_combined = populate_combined_drillcut_dataset(xr_dataset = ds_combined, config_obj = config_obj,folder = folder)
    export_combined_drillcut_dataset(xr_dataset = ds_combined, config_obj = config_obj)
    # clean the individual processed files (i.e. one per class, per release) to save space
    if remove_intermediates_files:
        clean_individual_processed_files(config_obj = config_obj,folder = folder)

    del(ds_combined) 
    gc.collect()



def stats_combined_drillcut_dataset(file_list = None,variable_to_process = None, output_filename = None, process_ssc = False):
    # compute a range of stats from the combined_drillcut_dataset file(s) in file_list.  
    # The combined_drillcut_dataset are produced by process_combined_drillcut_dataset()
     
    # thresholds and options - hardcoded for now
    deposition_thick_threshold_mm = [0.0025,0.05,1] # in mm
    ssc_stats_to_process = ['max','P99'] # ['mean','median','P99','max']
    ssc_threshold_mg_per_liters = [1,10,100] # in mg/L

    # create a xarray dataset where stats will be saved
    ds = xarray.open_dataset(file_list[0])
    ds_stats = xarray.Dataset(coords={'lon': (['x'], ds['lon'].data),
                    'lat': (['y'], ds['lat'].data),
                    'lon_corner': (['x_corner'], ds['lon_corner'].data), # saved for plotting purposes only
                    'lat_corner': (['y_corner'], ds['lat_corner'].data), # saved for plotting purposes only
                    'time': range(0,len(file_list))})  # time is actually event number here
    ds_stats.attrs = ds.attrs # save attributes
    ds_stats.attrs['file_list'] = file_list
    # make up a grid for distance computations
    xx,yy = np.meshgrid(ds_stats['lon'].data,ds_stats['lat'].data)
    xx = xx.ravel()
    yy = yy.ravel()

    # create variables and fill with zeros
    if variable_to_process is None:variable_to_process = ds.data_vars # process all variables by default
    for var in variable_to_process:
        if 'pdf' in var :
            if 'active' in var and process_ssc : 
                for stat_param in ssc_stats_to_process:
                    ds_stats[var + '_' + stat_param] = (['time','y','x'],  np.zeros( (len(file_list),ds['lat'].shape[0],ds['lon'].shape[0]) ) )
                    ds_stats[var + '_' + stat_param].attrs = ds[var].attrs
            else:
                    ds_stats[var + '_final'] = (['time','y','x'],  np.zeros( (len(file_list),ds['lat'].shape[0],ds['lon'].shape[0]) ) )
                    ds_stats[var + '_final'].attrs = ds[var].attrs
    # add threshold-dependent variables
    # deposition
    for thresh in deposition_thick_threshold_mm :
        thresh_str = str(thresh)
        ds_stats['pdf_settled_m2_above_' + thresh_str + 'mm'] = (['time'],  np.zeros( (len(file_list)))) # 10g/m2 ~ 0.0025 mm
        ds_stats['pdf_settled_m2_above_' + thresh_str + 'mm'].attrs = {'units' : 'm2'}
        ds_stats['pdf_settled_max_distance_' + thresh_str + 'mm'] = (['time'],  np.zeros( (len(file_list))))
        ds_stats['pdf_settled_max_distance_' + thresh_str + 'mm'].attrs = {'units' : 'm'}
    # add usual stats on deposition thickness e.g. RPS                
    ds_stats['m3_settled_in_frame'] = (['time'],  np.zeros( (len(file_list))))
    ds_stats['m3_settled_in_frame'].attrs = {'units' : 'm3'}
    ds_stats['pdf_settled_absolute_max'] = (['time'],  np.zeros( (len(file_list))))
    ds_stats['pdf_settled_absolute_max'].attrs = {'units' : 'm'}
    # ssc
    for thresh in ssc_threshold_mg_per_liters :
        for var in variable_to_process:
            if 'pdf_active' in var:
                thresh_str = str(thresh)
                ds_stats[var + '_m2_above_' + thresh_str + 'mg_per_liter'] = (['time'],  np.zeros( (len(file_list)))) # 10g/m2 ~ 0.0025 mm
                ds_stats[var + '_m2_above_' + thresh_str + 'mg_per_liter'].attrs = {'units' : 'm2'}

                ds_stats[var + '_max_distance_' + thresh_str + 'mg_per_liter'] = (['time'],  np.zeros( (len(file_list))))
                ds_stats[var + '_max_distance_' + thresh_str + 'mg_per_liter'].attrs = {'units' : 'm'}

    # loop through all files, and then variables to process stats 
    for ifile,fname in enumerate(file_list):
        ds = xarray.open_dataset(fname)
        print('loading %s' % (fname) )
        for var in variable_to_process:
            # print(var)
            if 'pdf' in var :
                if 'active' in var : 
                    if process_ssc:  # this is the time-consuming part of the processing, consider using dask ?? 
                        for stat_param in ssc_stats_to_process:
                            print(var + '_' + stat_param)
                            if stat_param == 'mean':
                                ds_stats[var + '_' + stat_param][ifile,:,:] = ds[var].mean(dim = 'time')
                            elif stat_param == 'median':
                                ds_stats[var + '_' + stat_param][ifile,:,:] = ds[var].median(dim = 'time')
                            elif stat_param == 'P99':
                                ds_stats[var + '_' + stat_param][ifile,:,:] = ds[var].quantile(q=0.99,dim = 'time',interpolation = 'higher')
                            elif stat_param == 'max':
                                ds_stats[var + '_' + stat_param][ifile,:,:] = np.max(ds[var].data,axis = 0) #ds[var].max(dim = 'time') 
                            else:
                                pass
                else : # stranded or settled, in that case want the last time step
                    ds_stats[var + '_final'][ifile,:,:] = ds[var].sel(time=ds['time'][-1])
        
        # add some diagnostic metrics
        # 
        # deposition
        ds_stats['m3_settled_in_frame'][ifile] = ds['pdf_settled'].sel(time=ds['time'][-1]).sum(dim='x').sum(dim='y')*ds.attrs['pixelsize_m']**2 # including bulking
        ds_stats['pdf_settled_absolute_max'][ifile] = ds['pdf_settled'].sel(time=ds['time'][-1]).max()
      
        for thresh in deposition_thick_threshold_mm :
            thresh_str = str(thresh)
            cell_above = np.where(ds_stats['pdf_settled_final'][ifile,:,:].data.ravel()*1e3 >= thresh) # find cells where thickness > thresh
            ref_site_lon = np.tile(ds_stats.attrs['center_point'][0],len(yy[cell_above]))
            ref_site_lat= np.tile(ds_stats.attrs['center_point'][1],len(yy[cell_above]))
            dist = distance(s_lat=yy[cell_above], s_lng=xx[cell_above], e_lat=ref_site_lat, e_lng=ref_site_lon)*1000 # compute distance from release
            if len(dist) == 0:dist = np.array(1.e36)
            ds_stats['pdf_settled_m2_above_' + thresh_str + 'mm'][ifile] = cell_above[0].shape[0] * ds.attrs['pixelsize_m']**2 # compute associated area using grid cell surface
            ds_stats['pdf_settled_max_distance_' + thresh_str + 'mm'][ifile]  = dist.max()
        
        if process_ssc:
            # ssc
            for thresh in ssc_threshold_mg_per_liters :
                for var in variable_to_process:
                    if 'pdf_active' in var:
                        for stat_param in ssc_stats_to_process:
                            thresh_str = str(thresh)
                            cell_above = np.where(ds_stats[var + '_' + stat_param][ifile,:,:].data.ravel()*1e3 >= thresh) # find cells where thickness > thresh
                            ref_site_lon = np.tile(ds_stats.attrs['center_point'][0],len(yy[cell_above]))
                            ref_site_lat= np.tile(ds_stats.attrs['center_point'][1],len(yy[cell_above]))
                            dist = distance(s_lat=yy[cell_above], s_lng=xx[cell_above], e_lat=ref_site_lat, e_lng=ref_site_lon)*1000 # compute distance from release
                            if len(dist) == 0:dist = np.array(1.e36)
                            ds_stats[var + '_m2_above_' + thresh_str + 'mg_per_liter'][ifile] = cell_above[0].shape[0] * ds.attrs['pixelsize_m']**2 # compute associated area using grid cell surface
                            ds_stats[var + '_max_distance_' + thresh_str + 'mg_per_liter'][ifile] =  dist.max()
        # import pdb;pdb.set_trace()

        # mass_suspended = ds['pdf_active[all]'].sum(dim='x').sum(dim='y')*ds['pdf_active[all]'].attrs['depth_band_thickness']*ds.attrs['pixelsize_m']**2
    
    if output_filename is not None:
        ds_stats.to_netcdf(path = output_filename)

    check_plot = False
    if check_plot:
        #####################################################################################################
        # DEPOSITION PLOTS
        #####################################################################################################
        # import matplotlib.pyplot as plt
        # plt.ion()
        # fig, ax = plt.subplots()
        # im = ax.pcolormesh(ds_stats['lon'],ds_stats['lat'],ds_stats['pdf_settled_final'][ifile,:,:])
        # ax.plot(xx[cell_above],yy[cell_above],'k.')
        # ax.plot(ds_stats.attrs['center_point'][0],ds_stats.attrs['center_point'][1],'r.')
        # cbar = fig.colorbar(im, ax=ax)

        import matplotlib.pyplot as plt
        plt.ion()
        fig, ax = plt.subplots()
        im = ax.pcolormesh(ds_stats['lon'],ds_stats['lat'],ds_stats['pdf_settled_final'].max(dim = 'time'))
        ax.contour(ds_stats['lon'],ds_stats['lat'],ds_stats['pdf_settled_final'].max(dim = 'time'),levels=[1*1e-3],colors='k')
        ax.contour(ds_stats['lon'],ds_stats['lat'],ds_stats['pdf_settled_final'].max(dim = 'time'),levels=[.1*1e-3],colors='r')
        ax.contour(ds_stats['lon'],ds_stats['lat'],ds_stats['pdf_settled_final'].max(dim = 'time'),levels=[0.05*1e-3],colors='g')
        ax.plot(ds_stats.attrs['center_point'][0],ds_stats.attrs['center_point'][1],'r.')
        cbar = fig.colorbar(im, ax=ax)

        import matplotlib.pyplot as plt
        plt.ion()
        fig, ax = plt.subplots()
        im = ax.pcolormesh(ds_stats['lon'],ds_stats['lat'],ds_stats['pdf_settled_final'].quantile(q=0.95,dim = 'time'))
        ax.plot(ds_stats.attrs['center_point'][0],ds_stats.attrs['center_point'][1],'r.')
        cbar = fig.colorbar(im, ax=ax)

        import matplotlib.pyplot as plt
        plt.ion()
        fig, ax = plt.subplots(4)
        ax[0].plot(ds_stats['m3_settled_in_frame']);ax[0].set_ylabel('m3_settled_in_frame')
        ax[1].plot(ds_stats['pdf_settled_m2_above_1mm']);ax[1].set_ylabel('pdf_settled_m2_above_1mm')
        ax[2].plot(ds_stats['pdf_settled_max_distance_1mm']);ax[2].set_ylabel('pdf_settled_max_distance_1mm')
        ax[3].plot(ds_stats['pdf_settled_absolute_max']);ax[3].set_ylabel('pdf_settled_absolute_max')

        #####################################################################################################
        # SSC PLOTS
        #####################################################################################################
        if process_ssc:
            ifile = 1
            import matplotlib.pyplot as plt
            plt.ion()
            fig, ax = plt.subplots()
            im = ax.pcolormesh(ds_stats['lon'],ds_stats['lat'],ds_stats['pdf_active[-25_-35]_max'][ifile,:,:])
            cbar = fig.colorbar(im, ax=ax)

            fig, ax = plt.subplots()
            im = ax.pcolormesh(ds_stats['lon'],ds_stats['lat'],ds_stats['pdf_active[mid+5_mid-5]_max'][ifile,:,:])
            cbar = fig.colorbar(im, ax=ax)

            fig, ax = plt.subplots()
            im = ax.pcolormesh(ds_stats['lon'],ds_stats['lat'],ds_stats['pdf_active[seafloor+10_seafloor]_max'][ifile,:,:])
            cbar = fig.colorbar(im, ax=ax)

            # max of all events
            fig, ax = plt.subplots()
            im = ax.pcolormesh(ds_stats['lon'],ds_stats['lat'],ds_stats['pdf_active[-25_-35]_max'].max(dim = 'time'))
            ax.contour(ds_stats['lon'],ds_stats['lat'],ds_stats['pdf_active[-25_-35]_max'].max(dim = 'time'),levels=[1*1e-3],colors='k')
            cbar = fig.colorbar(im, ax=ax)

            fig, ax = plt.subplots()
            im = ax.pcolormesh(ds_stats['lon'],ds_stats['lat'],ds_stats['pdf_active[mid+5_mid-5]_max'].max(dim = 'time'))
            cbar = fig.colorbar(im, ax=ax)

            fig, ax = plt.subplots()
            im = ax.pcolormesh(ds_stats['lon'],ds_stats['lat'],ds_stats['pdf_active[seafloor+10_seafloor]_max'].max(dim = 'time'))
            cbar = fig.colorbar(im, ax=ax)

    if False:
        # GEOTIFF
        # da = xarray.open_rasterio('/media/simon/Seagate Backup Plus Drive/metocean/aerials_geotiffs/NZ_500m_wgs84.tif')
        # img = plt.imshow(da[:3, :, :].data.transpose((1, 2, 0)))
        # img = plt.imshow(da.coords['x'].data,da.coords['y'].data,da[:3, :, :].data.transpose((1, 2, 0)))
        # https://rasterio.readthedocs.io/en/latest/topics/plotting.html
        import rasterio
        from rasterio.plot import show
        src = rasterio.open('/media/simon/Seagate Backup Plus Drive/metocean/aerials_geotiffs/NZ_500m_wgs84.tif')
        import matplotlib.pyplot as plt
        plt.ion()
        show(src.read(), transform=src.transform)
        im = plt.pcolormesh(ds['lon'],ds['lat'],ds['pdf_settled'][-1,:,:])
        plt.xlim(ds['lon_corner'].data[0]-1,ds['lon_corner'].data[-1]+1)
        plt.ylim(ds['lat_corner'].data[0]-1,ds['lat_corner'].data[-1]+1)


def plot_site_location(deotiff = None,site_xy = None):
    # 
    # To Tidy up
    # 
    # GEOTIFF
    # da = xarray.open_rasterio('/media/simon/Seagate Backup Plus Drive/metocean/aerials_geotiffs/NZ_500m_wgs84.tif')
    # img = plt.imshow(da[:3, :, :].data.transpose((1, 2, 0)))
    # img = plt.imshow(da.coords['x'].data,da.coords['y'].data,da[:3, :, :].data.transpose((1, 2, 0)))
    # https://rasterio.readthedocs.io/en/latest/topics/plotting.html
    import rasterio
    from rasterio.plot import show
    src = rasterio.open('/media/simon/Seagate Backup Plus Drive/metocean/aerials_geotiffs/NZ_500m_wgs84.tif')
    import matplotlib.pyplot as plt
    plt.ion()
    show(src.read(), transform=src.transform)
    im = plt.pcolormesh(ds['lon'],ds['lat'],ds['pdf_settled'][-1,:,:])
    plt.xlim(ds['lon_corner'].data[0]-1,ds['lon_corner'].data[-1]+1)
    plt.ylim(ds['lat_corner'].data[0]-1,ds['lat_corner'].data[-1]+1)
    

def distance(s_lat, s_lng, e_lat, e_lng):
    # https://gist.github.com/rochacbruno/2883505
    # approximate radius of earth in km
    R = 6373.0
    
    s_lat = s_lat*np.pi/180.0                      
    s_lng = np.deg2rad(s_lng)     
    e_lat = np.deg2rad(e_lat)                       
    e_lng = np.deg2rad(e_lng)  
    
    d = np.sin((e_lat - s_lat)/2)**2 + np.cos(s_lat)*np.cos(e_lat) * np.sin((e_lng - s_lng)/2)**2
    
    return 2 * R * np.arcsin(np.sqrt(d)) 

def process_oilspill_statistics(self):
    pass
    # see some of the updated matlab scripts as well as rps/apasa report for a list of useful metrics
    # threshold concentration on surface (for cleaning perspective) 10g/m2
    # ideally track down to 0.5g/m2 for deterministic simulations
    # see threshold in LeBreton
    # -discontinue weathering when oil reach shore
    # run until less than 5% is within model boundary ..at least duration+30days
    #
    # METRICS:  
    # - total length of shoreline oiled (above 10g/m2)
    # - total volume reaching the shore 
    # - locations of impact
    # - more metrics in RPS/APASA report
    # mini time for shoreline contact above threshold
    # mini time for surface concentration above threshold
    # proba of surface oiling above threshold
    # stats on time before oiling based on stochastic runs 
    # proba of oil beaching above threshold
    # cumulative oil beached etc..
    # RPS:
    # Surface oil average thickness >0.04 μm - available from opendrift
    # Shore oil average thickness >1.0 μm
    # Subsurface (within the water column) dissolved aromatic concentrations >1ppb

    # 
    # >> align with RPS/APASA
    # >> address all comments in review
    # 
    # FrenchMcCay 2016 Potential Effects Thresholds for Oil Spill Risk Assessments
    # 
    # ***consider not taking into account weathering for surface slicks and beaching for conservatism
    # but include weathering for the concentration of dissolved, dispersed etc...


#################################################################
# For testing
#################################################################
if __name__ == '__main__':

    pth = '/media/simon/Seagate Backup Plus Drive/metocean/0495_BeachEnergy_drillcut_oilspill/20050528_220020/'
    
    if False:
        ##############################################################
        # Single File Processing & Plotting                          #
        ##############################################################
        #create the OpenDriftPostProcess object
        o = OpenDriftPostProcess(opendrift_output_file = pth + 'wherry1_claysilt_fasterset_surface_20050528_220020.nc')
        # o.opendrift_object.plot() # using the Opendrift_object function
        # o.plot_simple(frame = [173.010988638889-2,173.010988638889+2,-45.0711106388889-2,-45.0711106388889+2])

        ds = o.compute_density_array(pixelsize_m = 1000,
                                frame = [173.010988638889-2,173.010988638889+2,-45.0711106388889-2,-45.0711106388889+2],
                                center_point = [173.010988638889,-45.0711106388889],
                                export_to_netcdf = True)

        # import pdb;pdb.set_trace()

        o.animate_density_array()
        o.plot_density_array()
        ##############################################################
    if False:
        ##############################################################
        # Single File Processing with pdf_options                    #
        ##############################################################
        pixelsize_m = 1000
        frame = [173.010988638889-1,173.010988638889+1,-45.0711106388889-1,-45.0711106388889+1]
        center_point = [173.010988638889,-45.0711106388889]
        pdf_options = { 'pdf_method' : 'numpy.histogram2d',
                        'vertical_levels' : [['all'],[-25,-50],['mid+10','mid-10'],['seafloor+10','seafloor']]}

        #create the OpenDriftPostProcess object
        o = OpenDriftPostProcess(opendrift_output_file = pth + 'wherry1_claysilt_fasterset_surface_20050528_220020.nc')                        
        ds = o.compute_density_array(pixelsize_m = 1000,
                                frame = [173.010988638889-1,173.010988638889+1,-45.0711106388889-1,-45.0711106388889+1],
                                center_point = [173.010988638889,-45.0711106388889],
                                export_to_netcdf = True,
                                pdf_options = pdf_options)
        del o
        del ds

        o = OpenDriftPostProcess(opendrift_output_file = pth + 'wherry1_claysilt_fasterset_nearbed_20050528_220020.nc')                                      
        ds = o.compute_density_array(pixelsize_m = 1000,
                                frame = [173.010988638889-1,173.010988638889+1,-45.0711106388889-1,-45.0711106388889+1],
                                center_point = [173.010988638889,-45.0711106388889],
                                export_to_netcdf = True,
                                pdf_options = pdf_options)

        import pdb;pdb.set_trace()

        o.animate_density_array()
        o.plot_density_array()
        ##############################################################
 
    if False :
        ##############################################################
        # Multi File Processing  : loop compute_density_array()                                   #
        ############################################################## 
        t0 = datetime.now()
        # import pdb;pdb.set_trace()

        fname_wildcard = 'wherry1_*.nc'
        file_list = glob.glob(pth+fname_wildcard)
        for f in file_list:
            o = OpenDriftPostProcess(opendrift_output_file = f)
            ds = o.compute_density_array(pixelsize_m = 1000,
                                    frame = [173.010988638889-2,173.010988638889+2,-45.0711106388889-2,-45.0711106388889+2],
                                    center_point = [173.010988638889,-45.0711106388889],
                                    export_to_netcdf = True)
            # clear variables for memory efficiency
            del ds
            del o
        t1 = datetime.now()
        print(t1-t0)
        #################################################################

    if False:
        # post process a single file for checking model for example
        fname = '/net/datastor1//data2/simon/0495_beachenergy_drillcutting/out_config_wherry_drillcutting/20000511_104655/wherry1_claysilt_fasterset_surface_20000511_104655.nc'
        fname = '/net/datastor1//data2/simon/0495_beachenergy_drillcutting/out_config_wherry_drillcutting/20000511_104655/wherry1_veryfinesand_surface_20000511_104655.nc'
        pixelsize_m = 200.0 # in meters
        square_frame_extent = 50000 # in meters = length/width of the square frame
        mean_lat = -45.0
        deg_lat = square_frame_extent/111000.0  # meters to degrees - length in meter of a degree latitude is constant
        deg_lon = deg_lat/np.cos(np.radians(mean_lat)) # length in meter of a degree longitude changes with latitude
        frame = [173.010988638889-deg_lon/2,173.010988638889+deg_lon/2,-45.0711106388889-deg_lat/2,-45.0711106388889+deg_lat/2]

        pdf_options = { 'pdf_method' : 'numpy.histogram2d',
                        'vertical_levels' : [['all'],[-25,-35],['mid+5','mid-5'],['seafloor+10','seafloor']]}

        o = OpenDriftPostProcess(opendrift_output_file = fname)
        ds = o.compute_density_array(pixelsize_m = pixelsize_m,
                                    frame = frame,
                                    center_point = [173.010988638889,-45.0711106388889],
                                    export_to_netcdf = True,
                                    pdf_options = pdf_options)

    if False :

        ####################################################################
        # Processing using a drill cutting config file               #
        # This will combine many different simulations into total 
        # suspended solid concentration [kg/m3] and deposition thickness [m]
        ####################################################################

        # Note all these tools could be moved to a Class inheriting from OpenDriftPostProcess()
        # e.g. OpenDriftPostProcess_DrillCuttings
        # with methods initialize_combined_drillcut_dataset(),add_waterdepth_to_dataset(),
        # populate_combined_drillcut_dataset(), load_drillcutting_config() etc... 
        # 
        # some of these functions may be useful for other cases e.g. oilspill so let's re-arrange further down the track


        # PROCESSING OPTIONS #####################################
        
        # load config file 
        config_file = './run_drillcut_configs/config_wherry_drillcutting.py'
        config_obj = load_drillcutting_config(config_file = config_file)
        
        # proba density function details
        # # test different frame/resolution
        # frame =  [173.010988638889-1,173.010988638889+1,-45.0711106388889-1,-45.0711106388889+1]
        # pixelsize_m = 1000.0
        # frame = [173.010988638889-.02,173.010988638889+.02,-45.0711106388889-.02,-45.0711106388889+.02]
        # pixelsize_m= 20.0
        pixelsize_m = 200 # in meters
        frame = [173.010988638889-.2,173.010988638889+.2,-45.0711106388889-.2,-45.0711106388889+.2]

        # to get a square grid, we need to work out frame extents based on latitude 
        pixelsize_m = 200 # in meters
        square_frame_extent = 40000 # in meters = length/width of the square frame
        mean_lat = -45.0
        # pixelsize_m = 20 # in meters
        # square_frame_extent = 4000 # in meters = length/width of the square frame
        DX_lat = square_frame_extent/111000.0  # meters to degrees - length in meter of a degree latitude is constant
        DX_lon = DX_lat/np.cos(np.radians(mean_lat)) # length in meter of a degree longitude changes with latitude
        frame = [173.010988638889-DX_lon/2,173.010988638889+DX_lon/2,-45.0711106388889-DX_lat/2,-45.0711106388889+DX_lat/2]
        import pdb;pdb.set_trace()

        center_point =  [config_obj['lon'],config_obj['lat']]
        pdf_options = { 'pdf_method' : 'numpy.histogram2d',
                        'vertical_levels' : [['all'],[-25,-35],['mid+5','mid-5'],['seafloor+10','seafloor']]}
                        
        # save these options to config_obj, they will be used in compute_density_array()
        config_obj['pixelsize_m'] = pixelsize_m
        config_obj['frame'] = frame
        config_obj['center_point'] = center_point
        config_obj['pdf_options'] = pdf_options
        ##########################################################

        # for local testing - only one folder 
        # config_obj['run_folder'] = '/media/simon/Seagate Backup Plus Drive/metocean/0495_BeachEnergy_drillcut_oilspill'
        # for local testing - add prefix to see data 
        config_obj['run_folder'] = '/net/datastor1/' + config_obj['run_folder']
        run_folder_list = glob.glob(os.path.join(config_obj['run_folder'],'out_' + config_obj['filename']) + '/*/')
        run_folder_list = glob.glob(os.path.join(config_obj['run_folder'],'out_' + config_obj['filename']) + '/20000426_123348/')
        # loop through folder list, and process combined dataset for each 
        for fold in run_folder_list:
            # clean_individual_processed_files(config_obj = config_obj,folder = fold) # to make sure we start fresh
            process_combined_drillcut_dataset(config_obj = config_obj,folder  = fold,remove_intermediates_files = True) # for each single event 

    if True :
        ####################################################################
        # Compute statistics from a list of "combined" datasets 
        # (from drill cutting simulation)
        ####################################################################
        fname = '/net/datastor1//data2/simon/0495_beachenergy_drillcutting/combined_pdfs_config_wherry_drillcutting/config_wherry_drillcutting_combined_pdfs_*_dx20m.nc'
        file_list = glob.glob(fname)
        stats_combined_drillcut_dataset(file_list = file_list,output_filename = '/media/simon/Seagate Backup Plus Drive/metocean/0495_BeachEnergy_drillcut_oilspill/wherry_20m_STATS_ssc.nc',process_ssc = True)

        # fname = '/net/datastor1//data2/simon/0495_beachenergy_drillcutting/combined_pdfs_config_wherry_drillcutting_respud/config_wherry_drillcutting_respud_combined_pdfs_*_dx20m.nc'
        # file_list = glob.glob(fname)
        # stats_combined_drillcut_dataset(file_list = file_list,output_filename = '/media/simon/Seagate Backup Plus Drive/metocean/0495_BeachEnergy_drillcut_oilspill/wherry_respud_20m_STATS_ssc.nc')