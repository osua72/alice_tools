import os,sys
import glob
import numpy as np

import xarray
from fiona import collection
#https://gis.stackexchange.com/questions/97545/using-fiona-to-write-a-new-shapefile-from-scratch/97563
from shapely.geometry import MultiPoint, Point, mapping
#https://shapely.readthedocs.io/en/latest/manual.html
# A nice little write up here:
#https://macwright.com/2012/10/31/gis-with-python-shapely-fiona.html
#
def var2zero(values,var):
        """Replace every 0 with 'zero' and return a copy."""
        values[ values<=var ]=0.0
        return values
        
method='DEP'
fdir = '/home/agowardbrown/Documents/Projects/0435_Nelson/opendrfit_wp3/Nino_%s_BHD/'%(method)
#import pdb; pdb.set_trace()
flist = np.sort(glob.glob(os.path.join(fdir,'*.nc')))
# fiona
#driver = 'ESRI_ShapeFile'
#source_crs = {'no_defs': True, 'ellps': 'WGS84', 'datum': 'WGS84', 'proj': 'longlat'}

schema = { 'geometry': 'Point', 'properties': {'units': 'str' , 'data':'float'} }
for file in flist:
    print(file)
    ds = xarray.open_dataset(file)
    #import pdb; pdb.set_trace()
    outdir = os.path.join(str.split(fdir,'/')[-2],str.split(str.split(file,'/')[-1],'.nc')[0])    
    for var in ds.data_vars:        
        outname = var
        # make directory if it doesn't exist already
        if os.path.isdir(outdir) is False:
            os.mkdir(outdir)

        data_len = int(ds['lon'].shape[0])*int(ds['lat'].shape[0])
        # sort data
        if method == 'SSC':
            if 'ssc_mean' in var:
                #import pdb; pdb.set_trace()
                vmin = 1.0
                vout = np.squeeze(ds[var].data)
                vtmp = var2zero(vout,vmin)
            if 'ssc_P90' in var:
                vmin = 1.0
                vout = np.squeeze(ds[var].data)
                vtmp = var2zero(vout,vmin)
            if 'ssc_P95' in var:
                vmin = 0.0
                vout = np.squeeze(ds[var].data)
                vtmp = var2zero(vout,vmin)
            if 'ssc_P99' in var:
                vmin = 0.0
                vout = np.squeeze(ds[var].data)
                vtmp = var2zero(vout,vmin)
            if 'ssc_exceed_10' in var:
                vmin = 0.0
                vout = np.squeeze(ds[var].data)
                vtmp = var2zero(vout,vmin)
            if 'ssc_exceed_50' in var:
                vmin = 0.0
                vout = np.squeeze(ds[var].data)
                vtmp = var2zero(vout,vmin)
            if 'ssc_exceed_100' in var:
                vmin = 0.0
                vout = np.squeeze(ds[var].data)
                vtmp = var2zero(vout,vmin)
            else: continue

        if method == 'DEP':
            if 'ssc_tot_settled_mean' in var:
                vmin = 0.0001
                vout = np.squeeze(ds[var].data)
                vout = var2zero(vout,vmin).reshape(data_len)
            elif 'ssc_tot_settled_P90' in var:
                vmin = 0.0001
                vout = np.squeeze(ds[var].data)
                vout = var2zero(vout,vmin).reshape(data_len)   
            elif 'ssc_tot_settled_P95' in var:
                vmin = 0.0001
                vout = np.squeeze(ds[var].data)
                vout = var2zero(vout,vmin).reshape(data_len)    
            elif 'ssc_tot_settled_P99' in var:
                vmin = 0.0001
                vout = np.squeeze(ds[var].data)
                vout = var2zero(vout,vmin).reshape(data_len)                                                
            else: continue
        # write to shapefile, in this case data needs to be written point by point, so the X,Y and Z data is reshaped 
        if len(vout.shape) == 3: # data is split into surface, middle and bottom arrays
            for i in range(0,len(vout.shape)):
                vlev = str.split(ds.attrs['vertical_levels'],'], [')[i]
                if vlev.startswith('[['): vlev = str.split(vlev,'[[')[-1]
                if vlev.endswith(']]'): vlev = str.split(vlev,']]')[0]
                vlev = str.replace(vlev,', ','_')
                vlev = str.replace(vlev,"'",'')
                shapename = ('%s/'+'%s_%s.shp')%(outdir,outname,vlev)
                try:
                    vout = np.squeeze(vtmp[i,:,:]).reshape(data_len)
                except: import pdb; pdb.set_trace()
                # one shapefile per variable
                with collection(
                    shapename, "w", "ESRI Shapefile", schema) as output:
                    #import pdb; pdb.set_trace()
                    [X,Y]=np.meshgrid(ds['lon'],ds['lat'])
                    for it in range(0,data_len):
                        #import pdb; pdb.set_trace()
                        if vout[it] > 0.0:
                            point = Point(float(X.reshape(data_len)[it]), float(Y.reshape(data_len)[it]))                    
                            output.write({
                                'properties': {'units': ds[var].attrs['units'],
                                               'data': vout.reshape(data_len)[it]},
                                'geometry': mapping(point)
                            })                
        else:
            shapename = ("%s"+'/'+"%s.shp")%(outdir,outname)
        # one shapefile per variable
        with collection(
            shapename, "w", "ESRI Shapefile", schema) as output:
                #import pdb; pdb.set_trace()
                [X,Y]=np.meshgrid(ds['lon'],ds['lat'])
                for it in range(0,data_len):
                    #import pdb; pdb.set_trace()
                    if vout[it] > 0.0:
                        point = Point(float(X.reshape(data_len)[it]), float(Y.reshape(data_len)[it]))                    
                        output.write({
                            'properties': {'units': ds[var].attrs['units'],
                                           'data': vout.reshape(data_len)[it]},
                            'geometry': mapping(point)
                        })


