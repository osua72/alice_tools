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
config.update({'area':1})
config.update({'pixelsize_m' : 30})
config.update({'center_point' : [176.965,-39.5778]})
config.update({'frame': [176.54,177.1140,-39.6823,-39.39]})
config.update({'start_datetime_filename' : '20020101_01'}
# set pdf options
pdf_options = { 'pdf_method' : 'numpy.histogram2d',
                'vertical_levels' : [['all']]}
config.update({'pdf_options': pdf_options})
