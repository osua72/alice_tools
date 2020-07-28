import sys, os
import netCDF4
import numpy as np 

sys.path.append('/data2/alice/alice_tools/python/opendrift/postprocessing')
import opendrift
from process_opendrift_wwtp import *
##########################################################
# set file paths
#config_file = '/data2/simon/github/toolbox_simon/python_tools/opendrift/run_drillcut_configs/config_wherry_newloc_drillcutting_respud.py'
#config_obj = load_drillcutting_config(config_file = config_file)
run_folder_list =['/data2/alice/Projects/0507_Hastings/opendrift/outputs/HB02']
area = int(sys.argv[1]) # which part of the pipe for the flow to be released from
# define points where outflow is estimated
coords = [[176.965,176.953,176.9467],[-39.5778,-39.582,-39.5842]] # release locations
# probability desnity function details
pdf_options = { 'pdf_method' : 'numpy.histogram2d',
                'vertical_levels' : [['all']]}
pixelsize_m = 20.0 # in meters
square_frame_extent = 10000 # in meters = length/width of the square frame
center_point = [coords[0][2],coords[1][2]]
print(center_point)
mean_lat = np.round(center_point[1])
print(mean_lat)
deg_lat = square_frame_extent/111000.0  # meters to degrees - length in meter of a degree latitude is constant
deg_lon = deg_lat/np.cos(np.radians(mean_lat)) # length in meter of a degree longitude changes with latitude
frame = [center_point[0]-deg_lon/2,center_point[0]+deg_lon/2,center_point[1]-deg_lat/2,center_point[1]+deg_lat/2]                               
##########################################################
# save these options to config_obj, they will be used in compute_density_array()
#config_obj['pixelsize_m'] = pixelsize_m
#config_obj['frame'] = frame
#config_obj['center_point'] = center_point
#config_obj['pdf_options'] = pdf_options
##########################################################
#
## change the ratios of different classes
# import pdb;pdb.set_trace()
#config_obj['filename'] =  config_obj['filename'] + '_new_ratio'
#
obj = OpenDriftPostProcess(opendrift_output_file = './outputs/HB02/opendrift_hastings_%s_20020601_010000.nc' % (area))

ds = obj.compute_density_array(pixelsize_m = pixelsize_m,
                               frame = frame, # None or [lon1,lon2,lat1,lat2] 
                               weight_name=None,  # use a variable to weight the pdf, expects variable name e.g. 'oil_mass''
                               center_point = None, # None or [lon,lat] used to adjust grid so that middle point is right on the user-input location
                               export_to_netcdf = True,
                               pdf_options = pdf_options)
                               #normalize_pdf = True)
##########################################################
# get volumes for the normalisation:

# calculate the dilution:
# 
fdir = './outputs/HB02/processed_pdfs/'#_largeframe/'
topo = 'chart-nz-561-approaches-to-napier_WGS.tif' #plot results over chart
### Dilution ###
# To get the dilution, the sum of the active and retired pdf array is 
# averaged and divided by the total number of particles released
### Concentration ###
# the sum of the active and retired particles is divided by the cell volume
# 
switch = 1
obj = OpenDriftPostProcess(opendrift_output_file = os.path.join(fdir,'opendrift_hastings_%s_20020601_010000_processed_dx%s.nc' %(area,int(pixelsize_m))))
# depths from log: 1: 10. 2: 8.83688 3. 5.0891
if area == 1:
        h = 10.
        r = 10.
elif area == 2:
        h = 8.8368
        r = 1.
elif area == 3:
        r = 1.
        h = 5.0891
vol_per_part = 10./h # normalise by nearfield concentration of partciles within the release location (nb_part/m3)
# volume per timestep = mass flowrate [m3/s] * timestep [s] 
if switch == 1:
        clabel = 'Dilution'#'conc' #'Dilution'
        stats_list=['p95','median','max']#['p95','median']
        print(stats_list)
        obj.make_map(topo=topo,
                        fdir='./',
                        fout=area,
                        clabel=clabel,
                        stats_list=stats_list,
                        tickMinutes=2.,
                        Err=None,
                        coords=center_point,
                        vol_per_part = vol_per_part)
##########################################################                    
#
# make TS plot
#
if switch == 0:
        POIS = [[177.071,177.023,176.931,176.943,176.991,176.984,176.960,176.963],
                [-39.636,-39.645,-39.567,-39.584,-39.628,-39.631,-39.579,-39.572]]
        fout = 'ts_%s'%(area)
        obj.make_ts_plot(coords=POIS,
                        fdir='./',
                        fout=fout,
                        vol_per_part=vol_per_part)

