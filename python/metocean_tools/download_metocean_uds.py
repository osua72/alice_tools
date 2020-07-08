#!/usr/bin/env python
import urllib
import requests
import os
from datetime import datetime, timedelta
import yaml
import logging

# "light" UDS class object to handle queries
class UDS_obj(object):
    def __init__(self,
        udshost='http://uds1.rag.metocean.co.nz:9191/uds',
        dset=['roms_cnz_surf'],
        var=['um','vm'],
        bnd=[173.51 ,174.5,-45.2,-45.4],
        time_start=datetime.now(),
        time_end=datetime.now()+ timedelta(days=1),
        fmt='nc',
        dt=3.0,
        dtype='hc',
        datum='msl',
        nomissing='True',
        cons = 'M2,S2,N2,K2,K1,O1,P1,Q1,MF,MM,M4,MS4,MN4',
        path_data = '/',
        output_filename=[]):
        super(UDS_obj, self).__init__()
        self.udshost=udshost
        self.dset=dset
        self.var=var
        self.bnd=bnd
        self.time_start=time_start
        self.time_end=time_end
        self.fmt=fmt
        self.dt=dt
        self.type=dtype
        self.datum=datum
        self.nomissing=nomissing
        self.cons = cons
        self.path_data = path_data
        self.output_filename=output_filename

    def query(self):
        # do the UDS query
        if not(os.path.exists(self.path_data)):
            os.mkdir(self.path_data)
        self.output_filename = [] # reset output_filename to [] 
        # 1) split time if required
        duration = self.time_end - self.time_start # timedelta
        # split in 15-day blocks 
        if duration.days > 32:
            print 'duration>32 days : split period in 15-day blocks '
            download_block = timedelta(days=15)
            t0 = self.time_start
            t1 = self.time_start+download_block
        else:
            t0 = self.time_start
            t1 = self.time_end
            download_block = timedelta(days=1000)
        # 2) convert list to string if required
        dset_str = list2string(self.dset)
        var_str = list2string(self.var)
        bnd_str = list2string(self.bnd)
        timestep_str = str(self.dt)
        fmt_str = self.fmt
        if hasattr(self,'cons'):
            cons_str = list2string(self.cons)

        while t1<=self.time_end : # time block loop
            time_str = t0.strftime('%Y%m%d.%H%M%Sz') + ',' + t1.strftime('%Y%m%d.%H%M%Sz')
            
            if 'tide' not in self.type: 
                syntax ='%s?&dset=%s&var=%s&bnd=%s&time=%s&fmt=%s&dt=%s&type=%s' % (self.udshost,dset_str,var_str,bnd_str,time_str,fmt_str,timestep_str,self.type)
            elif 'tide' in dset_str:
                syntax ='%s?&dset=%s&var=%s&bnd=%s&time=%s&fmt=%s&dt=%s&type=tide&cons=%s&datum=msl&nomissing=True' % (self.udshost,dset_str,var_str,bnd_str,time_str,fmt_str,timestep_str,cons_str)

            fname = '%s_' % (dset_str)  + t0.strftime('%Y%m%d') + '.nc'
            full_fname  = os.path.join(self.path_data,fname)
            logging.debug('WRAPPER:DOWNLOAD: Downloading ' + full_fname)
            # self.output_filename = os.path.join(self.path_data,fname)
            # check if file exists to avoid  re-downloading
            if os.path.exists(full_fname):
                logging.debug('WRAPPER:DOWNLOAD: already exists : skip to next file to download')
                self.output_filename.append(full_fname) # save name of output file
                pass
            else:
                self.output_filename.append(full_fname) # save name of output file
                logging.debug('WRAPPER:DOWNLOAD: UDS query : ' + syntax)
                print syntax
                urllib.urlretrieve(syntax, full_fname)
                logging.debug('WRAPPER:DOWNLOAD: Done')
            # move to next time block
            t0 = t0 + download_block
            t1 = t1 + download_block

def read_yaml_config(config_yaml_file):
    ''' read a YAML config file for OpenDrift - could be moved to opendrift/__init__.py eventually '''
    try:
        config = yaml.load(open(config_yaml_file).read())
        return config
    except Exception as ex: # handle arbitrary exception
        logging.error('WRAPPER:DOWNLOAD: ' + ex)
        logging.error('WRAPPER:DOWNLOAD: Cannot read ' + config_yaml_file)
        sys.exit('Cannot read ' + config_yaml_file)

def list2string(input_list):
    ''' turn string list into comma-separated string'''
    if isinstance(input_list,list):
        input_list_str=''
        for s in (input_list): 
            input_list_str = input_list_str+str(s)+','
        input_list_str=input_list_str[:-1]
        return input_list_str
    else:
        return input_list

def download(uds_host = 'http://uds1.rag.metocean.co.nz:9191/uds',
             path_data = './uds_data',
             datasets=['roms_nz_3D','nzra1_nz'], 
             variables=[['um','vm'],['ugrd10m','vgrd10m']],
             time_start=datetime.now(),
             time_end=datetime.now() + timedelta(days=1),
             timestep = 3.0,
             boundary=[173.51 ,174.5,-45.2,-45.4],
             datatype='fc,hc,tide',
             data_format= 'nc', 
             resolution = None,
             tide_cons='M2,S2,N2,K2,K1,O1,P1,Q1,MF,MM,M4,MS4,MN4',
             **kwargs):
             # datum='msl',
             # nomissing = True,
             # data_sorting = ['quality', '-res'],
             # spinup = '0.0',
             # stepback = 1,

    use_config = ('config_file_opendrift' in kwargs) or ('config_dict' in kwargs)

    if not use_config :
        # do the query based on function arguments 
        uds_obj = UDS_obj()
        uds_obj.time_start = time_start
        uds_obj.time_end = time_end
        uds_obj.path_data = path_data
        # for each dataset, update uds_obj, then query
        for ib,block in enumerate(datasets):
            uds_obj.udshost = uds_host
            uds_obj.dset=block # string 
            # check if variables is a list, or list of list (in case of multi dataset query)
            if any(isinstance(el, list) for el in variables) : 
                uds_obj.var=variables[ib] # read ib-th element of list variables
            else: 
                uds_obj.var=variables # # keep whole one-element list variables
            uds_obj.bnd=boundary
            uds_obj.dt =timestep
            uds_obj.type=datatype
            if uds_obj.type is 'tide':
                uds_obj.cons = tide_cons
            uds_obj.query() # query data for that UDS_obj
    # use of config file, or config dictionary to define download specs
    elif use_config: 
        if 'config_file_opendrift' in kwargs :
            # setup download using specifications in YAML configuration file
            config = read_yaml_config(kwargs['config_file_opendrift'])
        elif 'config_dict' in kwargs :
            # input was already a dictionary
            config = kwargs['config_dict']

        uds_obj = UDS_obj()
        uds_obj.time_start=datetime.strptime(config['start_time'],'%d-%m-%Y %H:%M')
        uds_obj.time_end=datetime.strptime(config['end_time_run'],'%d-%m-%Y %H:%M')
        # add some time before/after start
        uds_obj.time_start = uds_obj.time_start - timedelta(hours=6.0) 
        uds_obj.time_end = uds_obj.time_end + timedelta(hours=6.0) 
        uds_obj.path_data = os.path.join(config['rootdir'],'uds_data')
        # for each reader block - update uds_obj, then query
        for ib,block in enumerate(config['readers']):
            reader_block = config['readers'][block]
            if reader_block is not None and 'udshost' in config['readers'][block].keys():
                uds_obj.udshost = reader_block['udshost']    
                uds_obj.dset=reader_block['dset']
                uds_obj.var=reader_block['vars']
                uds_obj.bnd=reader_block['boundary']
                uds_obj.dt =reader_block['timestep']
                uds_obj.type=reader_block['datatype']
                # data_format= 'nc', 
                # resolution = None,
                # tide_cons='M2,S2,N2,K2,K1,O1,P1,Q1,MF,MM,M4,MS4,MN4',
                if 'constituents' in reader_block.keys():
                    uds_obj.cons = reader_block['constituents']
                 # query data for that UDS_obj
                uds_obj.query()
                # add output_filename to config['readers'][block]
                if len(uds_obj.output_filename)>1:# convert into a wildcard name
                    pth,f = os.path.split(uds_obj.output_filename[0])                  
                    f_wildcard = reader_block['dset'][0] + '_*.nc'#replace f by dataset_name_*.nc
                    config['readers'][block]['filename'] = [os.path.join(pth,f_wildcard)]
                else:
                    config['readers'][block]['filename'] = uds_obj.output_filename
        return config # return updated config object with readers filenames


#############################################################################
# Some Examples and tests
#############################################################################
if __name__ == '__main__':
 
    # using function as standalone
    if False: 
        # one dataset, one set of variables - hindcast
        download(uds_host = 'http://uds1.rag.metocean.co.nz:9191/uds',
                 datasets=['roms_cnz_surf'],
                 variables=['um','vm'],
                 time_start=datetime(2014,1,1),
                 time_end=datetime(2014,1,1)+timedelta(days=1),
                 timestep = 3.0,
                 boundary=[173.0,174.0,-41.0,-40.0],
                 datatype='hc')

        # 2 datasets, 2 sets of variables - hindcast
        download(uds_host = 'http://uds1.rag.metocean.co.nz:9191/uds',
                 datasets=[['roms_nz_surf'],['nzra1_nz']],
                 variables=[['umo','vmo'],['ugrd10m','vgrd10m']],
                 time_start=datetime(2014,1,1),
                 time_end=datetime(2014,1,1)+timedelta(days=1),
                 timestep = 3.0,
                 boundary=[173.0,174.0,-41.0,-40.0],
                 datatype='hc')

        # one dataset, one set of variable - tide
        download(uds_host = 'http://uds1.rag.metocean.co.nz:9191/uds',
                 datasets=['nz_tide'],
                 variables=['et','ut','vt'],
                 time_start=datetime.now(),
                 time_end=datetime.now() + timedelta(days=1),
                 timestep = 3.0,
                 boundary=[173.0,174.0,-41.0,-40.0],
                 datatype='tide',
                 tide_cons='M2,S2,N2,K2,K1,O1,P1,Q1,MF,MM,M4,MS4,MN4')

        # # one dataset, one set of variable - forecast
        # download(uds_host = 'http://uds-ops.service.consul:9191/uds',
        #          datasets=[['roms_cnz_surf'],['nzra1_nz']],
        #          variables=[['um','vm'],['ugrd10m','vgrd10m']],
        #          time_start=datetime(2014,1,1),
        #          time_end=datetime(2014,1,1)+timedelta(days=1),
        #          timestep = 3.0,
        #          boundary=[173.0,174.0,-41.0,-40.0],
        #          datatype='fc')
    
    # using a config file (to be used with OpenDrift )
    config_file_opendrift = 'C:\github\opendrift\examples_msl\config\opendrift.config_oceandrift3d_uds_download.yml'
    download(config_file_opendrift=config_file_opendrift)


#################################################################################
# Example script to use these functions
#################################################################################

# import sys
# sys.path.append('/home/metocean/github/toolbox_simon/python_tools')
# import download_metocean_uds
# from datetime import datetime,timedelta

# year = 2010
# month=range(1,13)

# for im,mth in enumerate(month):
#     download_metocean_uds.download(uds_host = 'http://uds1.rag.metocean.co.nz:9191/uds',
#              datasets=['roms_nz_3D'],
#              variables=['uo','vo','dep'],
#              time_start=datetime(year,mth,1),
#              time_end=datetime(year,mth+1,1),
#              timestep = 3.0,
#              boundary=[165.1,179.9,-47.9,-33.1],
#              datatype='hc',
#              path_data = './flow_fields')
#################################################################################