#!/usr/bin/env python2.7

import os,sys
import netCDF4
import numpy
import scipy.interpolate as spint
import scipy.spatial.qhull as qhull
import matplotlib.tri as mtri
import json


PREF='schout_'

def interp3D(zdata,udata,wanted_depth):
	m,n=zdata.shape
	zz=zdata.flatten(1)
	yinterp=numpy.zeros(shape=(n,1))
	yv=udata.flatten(1)
	for i in range(0,n):
		if((wanted_depth[i] < zz[((i*m)+0)]) or (wanted_depth[i] > zz[(i*m)+(m-1)])):
			yinterp[i] = numpy.NaN
		else:
			for j in range(0,m):
				if(wanted_depth[i]<=zz[(i*m)+j]):
					yinterp[i] = (wanted_depth[i]-zz[(i*m)+(j-1)]) / (zz[(i*m)+j]-zz[(i*m)+(j-1)]) * (yv[(i*m)+j]-yv[(i*m)+(j-1)]) + yv[(i*m)+(j-1)]; 
					           									
					break
				else:
					yinterp[i] = numpy.NaN
	return yinterp

def create_output_file(filout,params,depth,dtype):
	f = file(filout, "w")
	# header
	f.write('Year\tMonth\tDay\tHour\tMin\tSec\tDepth\t')
	col_tot=0
	val=-1

	for vars in params:      
		val=val+1
		if dtype[val][0]==2: # 3D variable 
			if dtype[val][1] ==2: # U and V data		
				for n in range(0,len(depth)):
					f.write('%s\t' % ('u_'+vars+'lev_'+str(depth[n])))
					col_tot=col_tot+1

				for n in range(0,len(depth)):
					f.write('%s\t' % ('v_'+vars+'_lev_'+str(depth[n])))
					col_tot=col_tot+1
			else:
				for n in range(0,len(depth)):
					f.write('%s\t' % (vars+'_lev_'+str(depth[n])))
					col_tot=col_tot+1

		else: # 2D variable
			if dtype[val][1] ==2: # U and V data
				f.write('%s\t%s\t' % ('u_'+vars,'v_'+vars))
				col_tot=col_tot+2
			else:
				f.write('%s\t' % (vars))
				col_tot=col_tot+1

	f.write('\n')
	return f,col_tot
 
def pair(arg):
    # For simplity, assume arg is a pair of integers
    # separated by a comma. If you want to do more
    # validation, raise argparse.ArgumentError if you
    # encounter a problem.
    return [float(x) for x in arg.split(',')]

def check_3d(params,dirin,Istart):
	# Quick check for 3D data
	dtype=[]
	for vars in params:
		ncfile=os.path.join(dirin,PREF+str(Istart)+'.nc')
		print ncfile
		ncs=netCDF4.Dataset(ncfile)   
		if 'two'  in ncs.variables[vars].dimensions: 
			B=2
		else:
			B=1

		if 'nSCHISM_vgrid_layers'  in ncs.variables[vars].dimensions:
			A=2
		else:
			A=1

		dtype.append([A,B])

		ncs.close()


	## extract static vars
	ncfile=os.path.join(dirin,PREF+str(Istart)+'.nc')
	ncs=netCDF4.Dataset(ncfile)
	X=ncs.variables['SCHISM_hgrid_node_x'][:]
	Y=ncs.variables['SCHISM_hgrid_node_y'][:]

	ele=ncs.variables['SCHISM_hgrid_face_nodes'][...,:3]-1
	triang = mtri.Triangulation(X, Y, ele)
	nt=len(ncs.variables['time'])
	ncs.close()

	return dtype,nt,X,Y,triang

def get_depth(x,y,triang,ncs):

    Z=ncs.variables['depth'][:]
    interp_lin = mtri.LinearTriInterpolator(triang, Z)
    return interp_lin(x,y).data

def process(POS,fileout,dirin,Istart,Iend,params,depth):

	dtype,nt,X,Y,triang=check_3d(params,dirin,Istart)


	for k in range(0,len(POS[0]),2):
		x=POS[0][k]
		y=POS[0][k+1]
		# create the output file (one for each point)
		fout=fileout+'_'+str(x)+'_'+str(y)+'.txt'    
		print 'Extract file : '+fout
		f,col_tot=create_output_file(fout,params,depth,dtype)

#### For each of the file
		for nfile in range(Istart,Iend+1):
		
			Total= numpy.zeros(shape=(nt,col_tot+7))
			print ' read file: '+str(nfile)
#### For each variable
			val=0

			for i,vars in enumerate(params):
				print '		Read %s' % (vars)
				ncfile=os.path.join(dirin,PREF+str(nfile)+'.nc')
				ncs=netCDF4.Dataset(ncfile)

			
				if val==0:
					dtime = netCDF4.num2date(ncs.variables['time'][:],ncs.variables['time'].units)
					Total[:,0]=[dt.year for dt in dtime.astype(object)]
					Total[:,1]=[dt.month for dt in dtime.astype(object)]
					Total[:,2]=[dt.day for dt in dtime.astype(object)]
					Total[:,3]=[dt.hour for dt in dtime.astype(object)]
					Total[:,4]=[dt.minute for dt in dtime.astype(object)]
					Total[:,5]=[dt.second for dt in dtime.astype(object)]
					Total[:,6]=get_depth(x,y,triang,ncs)
					val=val+6


				data=ncs.variables[vars][:]
				data = numpy.ma.masked_where(data==-9999,data)
### For each timestep
				if dtype[i][0]==2: # 3D data
					Zdata=ncs.variables['zcor'][:]
					Zdata = numpy.ma.masked_where(Zdata==-9999,Zdata)			

 					for Nt in range(0,nt):	
 						if Nt==0:
							val=val+1
							
 						zdata=Zdata[Nt,:]
 						Data=data[Nt,:]
						
 						Vdata=numpy.zeros(shape=(Data.shape[1],1))
						if 'two'  in ncs.variables[vars].dimensions:
							Udata=numpy.zeros(shape=(Data.shape[1],1))
	
 						ZZdata=numpy.zeros(shape=(Data.shape[1],1))
 						for lev in range(0,Data.shape[1]):
 												 						
							if 'two'  in ncs.variables[vars].dimensions:
								interp_lin = mtri.LinearTriInterpolator(triang,Data[:,lev,0])
                                                        	Vdata[lev]=interp_lin(x,y).data
								interp_lin = mtri.LinearTriInterpolator(triang,Data[:,lev,1])
                                                        	Udata[lev]=interp_lin(x,y).data
							else:
								interp_lin = mtri.LinearTriInterpolator(triang,Data[:,lev])
                                        		        Vdata[lev]=interp_lin(x,y).data
								
 							interp_lin = mtri.LinearTriInterpolator(triang,zdata[:,lev])
 							ZZdata[lev]=interp_lin(x,y).data
						
 						ZZdata = numpy.ma.masked_where(ZZdata==-9999,ZZdata)
						ZZdata = numpy.ma.masked_where(ZZdata>1e36,ZZdata)
 						Vdata = numpy.ma.masked_where(Vdata==-9999,Vdata)
						Vdata = numpy.ma.masked_where(Vdata>1e36,Vdata)
						if 'two'  in ncs.variables[vars].dimensions:
							Udata = numpy.ma.masked_where(Udata==-9999,Udata)
							Udata = numpy.ma.masked_where(Udata>1e36,Udata)
 						nn=-1
						for n in range(0,len(depth)):
							nn=nn+1	
											
 							if depth[n]>0: # above sea bed
 							      wanted_depth=numpy.min(ZZdata,0)+depth[n];
 							if depth[n]<=0: # below sea surface
 							      wanted_depth=numpy.max(ZZdata,0)+depth[n];
						
 							Total[Nt,val+nn]=interp3D(ZZdata,Vdata,wanted_depth)
							if 'two'  in ncs.variables[vars].dimensions:
								Total[Nt,val+len(depth)+nn]=interp3D(ZZdata,Udata,wanted_depth)
								
                                                        	
                                                                	

				else:

					for Nt in range(0,nt):
						if Nt==0:
							val=val+1

						Data=data[Nt,:]
						if 'two'  in ncs.variables[vars].dimensions:
							Uinterp = mtri.LinearTriInterpolator(triang,Data[...,0])
							Vinterp = mtri.LinearTriInterpolator(triang,Data[...,1])
							Total[Nt,val]=Uinterp(x,y).data
							Total[Nt,val+1]=Vinterp(x,y).data
							if Nt==nt-1:
								val=val+1
						else:
							Uinterp = mtri.LinearTriInterpolator(triang,Data)
							Total[Nt,val]=Uinterp(x,y).data


					
			ncs.close()			



			numpy.savetxt(f,Total,fmt='%g',delimiter='\t')
	f.close()


     
  
    
    
if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(prog='extract_point.py', usage='%(prog)s fileout dirout INDstart INDend params (i.e python extract_point.py caca /home/remy/Buisness/Ciberon/ 15 16 ''{"elev":"elev","temp":"temp"}'' -POS 256842,9291527 -depth -5 -4')
	## main arguments
	
	parser.add_argument('fileout', type=str,help='name of the output file (without the extension)')
	parser.add_argument('dirout', type=str,help='name of the output where are the SELFE files')
	parser.add_argument('INDstart', type=int,help='First file to take')
	parser.add_argument('INDend', type=int,help='Last file to take')
	parser.add_argument('params', type=str,nargs='+',help='name of the parameter to plot')
	parser.add_argument('-POS', type=pair,nargs='+',help='XY position to extract')
        parser.add_argument('-depth', type=float,nargs='+',help='Z position to extract must be negative for bsl or pos for asb')
	args = parser.parse_args()
	

	### PRINT IT ALL
	print 'output name : %s' % (args.fileout)
	print 'Direcory : %s' % (args.dirout)
	print 'From file #%i and #%i' % (args.INDstart,args.INDend)
	print 'Do parameters : %s' % (args.params)

	
	
	process(args.POS,args.fileout,args.dirout,args.INDstart,args.INDend,args.params,args.depth)

