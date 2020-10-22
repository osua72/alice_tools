###############################
# PROJECT-SPECIFIC INFO
# 
# configuration is saved as a python dictionary, called config
# 
###############################
# OPENDRIFT CONFIGS
###############################
import numpy as np
config = {'site_name' : 'napier'}
config.update({'area' : 1})
config.update({'lon' : 176.965})
config.update({'lat' :-39.5778})
config.update({'z' : 'z_rel_watercolumn'}) 
config.update({'rad' : 1}) # radius of release

# file name and number of particles to release over the duration of the run
config.update({'start_datetime_filename' : '20020101_01'})
config.update({'start_time':151}) # days after input file start date
config.update({'nb_parts' : 288000}) # total number of run with the period [year_start_stochastic,year_end_stochastic] 
config.update({'max_run_length_days': 30.+14.}) # total run duration in days, approx.

# sort out projection:
config.update({'proj_epsg_2135': '+proj=utm +zone=60 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs '})
# land mask
config.update({'use_custom_landmask': True}) # use by default
config.update({'custom_coastfile': './coastline/napier_shoreline_opendrift.shore'})

# PARTICLE RELEASE
ptm_release = {
    'rel_depth': ['z_rel_watercolumn'],
    'duration_days' : [30.0]} # days to run model for
config.update({'ptm_release':ptm_release})    

config.update({'ocean_horizontal_diffusivity' :  0.01}) #Okubo for 30m = 0.01, for 60 0.02, for 100m 0.04 
config.update({'ocean_vertical_diffusivity' : 0.0001}) #>> to reduce too # m2/s-1
config.update({'drift:current_uncertainty': np.sqrt(2*0.01/600.) }) # 2*ohd/600.
config.update({'wind_drift_factor' : 0.03}) # ~3% as per NOAA's GNOME\
config.update({'drift:max_age_seconds': 30*24*3600}) # 30 days

# DATASETS TO BE USED AS FORCING
# Required variables
# 'x_sea_water_velocity', 'y_sea_water_velocity', 'sea_floor_depth_below_sea_level',
# 'sea_surface_wave_significant_height', 'sea_surface_wave_period_at_variance_spectral_density_maximum',
# 'x_wind', 'y_wind', 'land_binary_mask'
# 
# and ideally these too:
# 
# 'sea_water_temperature', 'sea_water_salinity'
# 'sea_surface_wave_stokes_drift_x_velocity', 'sea_surface_wave_stokes_drift_y_velocity'
# 
# required variables that are not available should be specifiec as :
# config.update({'sea_water_temperature' : 18.0}) 
# config.update({'sea_water_salinity' : 18.0}) 
ocean_forcing = {
    'path': '/data2/alice/Projects/PXXX_Napier/HB02/',
    'filename': 'schism_napier20020101_01z_3D.nc',
    'reader_type' : 'unstructured_selfe', # not currently used
}
config.update({'ocean_forcing':ocean_forcing})

# The general lanmask will be sorted out according to reader coverage (if not using the native_landmask)

# run folder used for simulation in servers, in /data2
config.update({'run_folder': '/data2/alice/Projects/PXXX_Napier/'})
