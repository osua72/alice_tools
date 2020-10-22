#!/usr/bin/env python
# docker run -it  --rm -v /data2:/data2 -v /data2/simon/github/opendrift:/code/ -v /net/datastor1:/net/datastor1 opendrift/opendrift:py3-v1.0.7
# cd /data2/alice/Projects/PXXX_Napier
import sys, os
import importlib

#make high resolution coastline 
if False:
    # shape file - coast
    # this shape file needs to contain closed polygons of the land and island regions
    # I tend to load in the hgrid.gr3 schism file into QGmesh and extract the coastline files from this,
    # and make sure they are all closed before exporting as a shapefile
    shapein = './coastline/hastings_shoreline_wgs.shp' # check this
    coastout = './coastline/napier_shoreline_opendrift.shore'     
    # packages required
    import os, sys
    import shapefile
    import numpy 
    import scipy
    def write_coast_opendrift(shapein,coastout):
        # filenames to use
        sf = shapefile.Reader(shapein)
        shapes = sf.shapes()
        with open(coastout,'w') as f:
            for i in range(0,len(shapes)):
                f.write("NaN"+' '+"NaN"+'\n')
                for point in shapes[i].points:
                    f.write(str(point[0])+' '+str(point[1])+'\n')

            f.write("NaN"+' '+"NaN"+'\n')
        f.close()
        print("Done! Check for file: ",coastout) 
    if __name__ == '__main__':
        write_coast_opendrift(shapein,coastout)

# opendrift model nomenclature: schism_northport20150101_01z_3D.nc etc
# need to pseudo link files to opendrift directory and rename:
if False: # 
    # packages required
    import os,sys,glob
    import netCDF4
    import numpy as np
    import datetime
    # variables required:
    schism_dir = '/data2/simon/hastings_flow_fields'
    mon_year = ['HB02']
    outname = 'schism_napier'
    # define local directory
    fdir = '/data2/alice/Projects/PXXX_Napier/' #Existing/schis
    # define opendrift input directories
    schism_fle = 'selfie_napier20020101_01z_3D.nc' # in this case it's just one file, sometimes we use multiple output files
    #  define opendrift output directories
    outdir = fdir+'%s/'%(mon_year[0])

    if os.path.isdir(outdir) is False:
        os.mkdir(outdir)

    def symlink_schism(mon_year,schism_dir,out_dir,outname,fdir):
        # for each folder, loop through and get each file to create a symbolic link with
        for i in range(0,len(mon_year)):
            out_dir = os.path.join(schism_dir,mon_year[i])
            if os.path.isdir(schism_dir) is True:
                if os.path.isdir(os.path.join(out_dir,'outputs')) is True:
                    flist=np.sort(os.listdir(os.path.join(out_dir,'outputs')))                
                    for file in flist:
                        print(file)
                        # get start time for each file and write into symbolic link file name:
                        nc=netCDF4.Dataset(os.path.join(out_dir,'outputs',file))
                        t0=netCDF4.num2date(nc.variables['time'][0],nc.variables['time'].units)
                        nc.close()
                        # make ready for output
                        fout = '%s%i%02i%02i_%02iz_3D.nc'%(outname,t0.year,t0.month,t0.day,t0.hour)
                        symout = os.path.join(fdir,'schism',fout)
                        schism_file = os.path.join(out_dir,'outputs',file)
                        # create symbolic link
                        os.system('ln -s %s %s' %(schism_file,symout))
                else:  
                    flist=np.sort(os.listdir(schism_dir)) 
                    idx = np.flatnonzero(np.core.defchararray.find(flist,mon_year[0])!=-1)[0]
                    file = flist[idx]
                    print(file)
                    # get start time for each file and write into symbolic link file name:
                    nc=netCDF4.Dataset(os.path.join(schism_dir,file))
                    try:
                        t0=netCDF4.num2date(nc.variables['time'][0],nc.variables['time'].units)
                    except:
                        year = '20'+mon_year[0][-2::]
                        month = '01'
                        day = '01'
                        start_date = datetime.datetime.strptime(year+month+day,'%Y%m%d')
                        t0 = netCDF4.num2date(nc.variables['time'][0],'seconds since %s'%(start_date))
                    nc.close()
                    # make ready for output
                    fout = '%s%i%02i%02i_%02iz_3D.nc'%(outname,t0.year,t0.month,t0.day,t0.hour)
                    symout = os.path.join(outdir,fout)
                    schism_file = os.path.join(schism_dir,file)
                    # create symbolic link
                    #import pdb; pdb.set_trace()
                    os.system('ln -s %s %s' %(schism_file,symout))                    

    if __name__ == '__main__':
        symlink_schism(mon_year,schism_dir,outdir,outname,fdir)                    

# run opendrift
if False:
    #!/usr/bin/env python
    import sys
    import datetime
    import numpy
    import os
    import importlib

    #################################################################
    # The script runs a wwtp event for a specified release location and duration
    # 
    #################################################################
    # Usage :
    # opendrift_napier_WWTP.py <config_filename>
    # 
    # config_filename is a <.py> file that includes a dictionary 
    # with all simulation detail. That file is expected to be found 
    # in subfolder ./run_wwtp_configs
    # 
    # e.g. 
    # opendrift_napier_WWTP.py ./run_wwtp_configs/config_napier_wwtp.py
    # 
    # AMMEND the below to suit different run locations? i/e Area, LAT LON
    # To able to split runs across several nodes, using the same config file,
    # the config parameters 'year_start_stochastic' and 'year_end_stochastic'
    # can be input at the function as follow (this will overrid parameters in config file):
    # 
    # opendrift_napier_WWTP.py ./run_wwtp_configs/config_napier_wwtp.py 1 -39 173.5
    # 
    # 
    #################################################################

    def run_opendrift_wwtp(config = None):
        import sys,os
        sys.path.append('/data2/alice/alice_tools/python/opendrift/postprocessing')
        from opendrift.readers import reader_netCDF_CF_unstructured_selfe,reader_global_landmask
        from opendrift.models.oceandrift3D import OceanDrift3D
        from opendrift.readers import reader_landmask_custom
        import numpy as np
        from datetime import datetime, timedelta
        import time
        import glob

        ###############################
        # READERS
        ############################### 
        # OCEAN FORCING
        reader_fname = config['ocean_forcing']['path'] + config['ocean_forcing']['filename']
        reader_ocean = reader_netCDF_CF_unstructured_selfe.Reader(reader_fname, proj4 = config['proj_epsg_2135'], use_3d = True, name='reader_ocean')
        #
        # MASKING
        # check if we use the native landmask from hydro files, or reader_global_landmask from opendrift
        if config['use_custom_landmask'] is True:
            use_global_landmask = False
            reader_custom = reader_landmask_custom.Reader(polygon_file = config['custom_coastfile'])
        else:
            use_global_landmask = True # use by default
            if 'use_native_landmask' in config['ocean_forcing']:
                use_global_landmask = not config['ocean_forcing']['use_native_landmask']
            if use_global_landmask:
                # Landmask - make it fit with largest reader
                reader_landmask = reader_global_landmask.Reader(
                                llcrnrlon=corners[0][0], llcrnrlat=corners[1][0],
                                urcrnrlon=corners[0][1], urcrnrlat=corners[1][1])

        # check if we use multiprocessing or not
        if 'use_multiprocessing' not in config['ocean_forcing']:
            config['ocean_forcing']['use_multiprocessing'] = True # use multiprocessing by default ( multiprocessing is used in lonlat2xy() )
        if not config['ocean_forcing']['use_multiprocessing']:
            reader_ocean.multiptocessing_fail = True # this reader attribute will prevent use of multiprocessing - see basereader.py line796

        ###############################
        # MODEL
        ###############################
        simulation_name = config['site_name']+'_'+config['start_datetime_filename']
        o = OceanDrift3D(loglevel = 0,logfile = simulation_name + '.log')   # Set loglevel to 0 for debug information 
        o.set_projection(proj4='+proj=latlong') # coordinate system used to run simulations
        ###############################
        # READERS
        ############################### 
        if use_global_landmask: 
            o.add_reader([reader_landmask,reader_ocean,reader_wave1,reader_wave2,reader_wind1,reader_wind2])
        else:
            o.add_reader([reader_custom,reader_ocean])
            o.set_config('general:use_auto_landmask', False) # prevent opendrift from make a new dynamical landmask
            o.fallback_values['land_binary_mask'] = 0 # set anything beyond the mask to ocean 

        ###################################
        # Define Particle release locations
        ###################################
        release_depth = o.get_environment(['sea_floor_depth_below_sea_level'],time=reader_ocean.start_time, lon=np.array((config['lon'],)),lat=np.array((config['lat'],)),z=np.zeros(1), profiles=None)
        site_release_depth = np.zeros(1)
        for i in range(0,len(site_release_depth)): # leave this in, incase future scenarios want a release at multiple sites along the pipe
            site_release_depth[[i]] = [release_depth[0][i][0]]

        print(config['lon'],' ',config['lat'],' ','SR depth: ',site_release_depth)

        # for release site and set day increments
        parts_per_rel = config['nb_parts']
        #
        lon_rel = config['lon']
        lat_rel = config['lat']
        if 'z_rel_watercolumn' in config['ptm_release']['rel_depth']:
            z_rel = np.random.uniform(-np.abs(site_release_depth)+1.0,0,size=int(parts_per_rel*1.0))
            
        ##################
        # Define timings
        ##################
        t0 = reader_ocean.start_time + timedelta(days = config['start_time']) # start of release sequence - 1/6/2002
        print('TO: ', t0)
        #
        t_end = t0 + timedelta(days = config['ptm_release']['duration_days'][0]) # month long release (2 x spring neap cycles, El Nino)
        ###############################
        # PARTICLE SEEDING
        ###############################
        # seeding the dredging part - bottom (due to drag head) / bucket
        o.seed_elements(lon = config['lon'],
                        lat = config['lat'],
                        number = config['nb_parts'],
                        z = z_rel,
                        radius=config['rad'] ,
                        time= [t0,t_end])     
        # release_start_i = config['start_time_stochastic'] + datetime.timedelta(days = days_since_t0)
        run_dir = os.path.join(config['run_folder'],'out_' + config['start_datetime_filename'])
        if not os.path.exists(run_dir): # create run folder if it doesnt exist yet
            print('Creating run folder : %s' % (run_dir))
            os.mkdir(run_dir) # one folder per event
        print('Running event %s' % (t0))
        # define name of that single simulation
        scenario = config['site_name'] + '_'+ str(config['area']) + '_' + str(t0.strftime('%Y%m%d_%H%M%S'))
        out_file = os.path.join(run_dir,scenario + '.nc')

        if os.path.exists(out_file):
            print('%s already exists  - skipping' %  (out_file) )
            return
        else:
            print('Running %s' % (out_file))

        ###############################
        # PHYSICS
        ###############################
        o.list_config()
        o.list_configspec()
        #
        o.fallback_values['x_sea_water_velocity'] = 0
        o.fallback_values['y_sea_water_velocity'] = 0
        o.fallback_values['x_wind'] = 0
        o.fallback_values['y_wind'] = 0
        o.fallback_values['sea_floor_depth_below_sea_level'] = 100.0

        dt_seconds = 900.0 # hard-coded for now
        # diffusion -
        o.fallback_values['ocean_vertical_diffusivity'] = config['ocean_vertical_diffusivity']
        # current uncertainty can be set to be equivalent to horizontal diffusion coefficient as follow:
        #o.set_config('drift:current_uncertainty', (2*config['ocean_horizontal_diffusivity']/dt_seconds)**0.5 )
        # drift
        o.set_config('drift:scheme','runge-kutta4') # or 'runge-kutta'
        # o.get_config('drift:max_age_seconds',config['max_run_length_days']*3600)
        o.set_config('drift:wind_uncertainty', 0.0)
        o.set_config('drift:stokes_drift', False) # not available in wave files
        o.set_config('drift:use_tabularised_stokes_drift', False) # not using estimations from wind for now.
        o.set_config('drift:tabularised_stokes_drift_fetch','25000')

        o.set_config('processes:verticaladvection' , False) # no vertical current available, so no vertical advection
        # vertical mixing options
        o.set_config('processes:turbulentmixing', True)
        o.set_config('turbulentmixing:timestep',  dt_seconds) # use same as simulation dt , unless we have diffusivity profiles - not yet functional
        o.set_config('turbulentmixing:diffusivitymodel', 'environment') # i.e. specified from model or constant
        o.set_config('turbulentmixing:TSprofiles',False)
        o.set_config('turbulentmixing:timestep', 90.0)  # if some ocean_vertical_diffusivity!=0, turbulentmixing:timestep should be less than 900 seconds (15min)    
        # 
        o.set_config('general:coastline_action','previous') # option('none', 'stranding', 'previous', default='stranding')
        ###############################
        # RUN 
        ############################### 
        o.run(stop_on_error = True,
            time_step = dt_seconds, 
            end_time = t0 + timedelta(days = config['max_run_length_days']), 
            outfile = out_file,
            time_step_output = 1800.0,
            export_variables = ['x_sea_water_velocity','y_sea_water_velocity','z','sea_floor_depth_below_sea_level'])

    if __name__ == '__main__':
        #################################################################
        # load the configuration file
        #################################################################
        config_file = sys.argv[1] # get configuration file as input
        # add the folder where all configs are expected to be found to path
        sys.path.append('./run_wwtp_configs')
        config_file = os.path.basename(config_file) # keep only filename
        # remove the .py extension if present
        if config_file[-3:] == '.py':
            config_file = config_file[:-3]
        config = importlib.import_module(config_file)
        config = config.config # keep only the config dictionary
        config['config_file'] = config_file # save config filename
        print('Loaded configuration file %s ...' % (config_file))
        if len(sys.argv) ==4: 
            # enable overriding of area, lon and lat so that several jobs can be submitted
            # This is useful to submit several jobs on many nodes, using command line but same config file
            # e.g. opendrift_run_drillcuttings_stochastic.py config_file 2001 2003
            config['area'] = int(sys.argv[2])
            config['lon'] = int(sys.argv[3])
            config['lat'] = int(sys.argv[4])

        #################################################################
        # create folder where outputs of each simulation will be saved
        #################################################################   
        out_dir = os.path.join(config['run_folder'],'out_' + config['start_datetime_filename'])
        if not os.path.exists(out_dir): 
            os.mkdir(out_dir)
        #################################################################
        # start the loop to run all events
        #################################################################
        print('Starting model run loop...')
        run_opendrift_wwtp(config = config)
        print('Model run loop completed')    

# post process
if False:
    # To run:
    # pip install xarray - opendrift doesn't run with xarray installed currently, bu tneeded for post-processing
    # pip install cmocean
    # python opendrift_napier_WWTP.py run_wwtp_configs/config_post_process.py
    import sys
    sys.path.append('/data2/alice/alice_tools/python/opendrift/postprocessing')
    import datetime
    import numpy
    import os
    import importlib    
    from process_opendrift_wwtp import *
    import glob
    def postprocess_opendrift(config=None):
        obj = OpenDriftPostProcess(opendrift_output_file = glob.glob(os.path.join(config['run_folder'],'out_'+config['start_datetime_filename'],config['site_name']+'_%s_*.nc'%(config['area'])))[0])

        if False:
            obj.plot_particle_clouds(frame = frame,filename = 'test.png')

        ds = obj.compute_density_array(pixelsize_m = config['pixelsize_m'], 
                                    frame = config['frame'], # None or [lon1,lon2,lat1,lat2] 
                                    weight_name=None,  # use a variable to weight the pdf, expects variable name e.g. 'oil_mass''
                                    center_point = config['center_point'], # None or [lon,lat] used to adjust grid so that middle point is right on the user-input location
                                    export_to_netcdf = True,
                                    pdf_options = config['pdf_options'])
                                    # normalized_pdf = True

    if __name__ == '__main__':
        #################################################################
        # load the configuration file
        #################################################################
        config_file = sys.argv[1] # get configuration file as input
        # add the folder where all configs are expected to be found to path
        sys.path.append('./run_wwtp_configs')
        config_file = os.path.basename(config_file) # keep only filename
        # remove the .py extension if present
        if config_file[-3:] == '.py':
            config_file = config_file[:-3]
        config = importlib.import_module(config_file)
        config = config.config # keep only the config dictionary
        config['config_file'] = config_file # save config filename
        print('Loaded configuration file %s ...' % (config_file))
        if len(sys.argv) ==3: 
            # enable overriding of area, lon and lat so that several jobs can be submitted
            # This is useful to submit several jobs on many nodes, using command line but same config file
            # e.g. opendrift_run_drillcuttings_stochastic.py config_file 2001 2003
            config['area'] = int(sys.argv[2])
            config['centre_point'] = [int(sys.argv[3]),int(sys.argv[4])]

        #################################################################
        # create folder where outputs of each simulation will be saved
        #################################################################   
        out_dir = os.path.join(config['run_folder'],'out_' + config['start_datetime_filename'],'processed')
        if not os.path.exists(out_dir): 
            os.mkdir(out_dir)
        #################################################################
        # start the loop to run all events
        #################################################################
        print('Starting model run loop...')
        postprocess_opendrift(config = config)
        print('Model run loop completed') 
