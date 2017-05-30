import numpy as np
import netCDF4 as nc
from os import listdir
import os
import argparse


parser = argparse.ArgumentParser(prog='PyCLES')
parser.add_argument("path")
args = parser.parse_args()
path = args.path
# path = '/Users/yaircohen/Documents/PyCLES_out/Output.TRMM_LBA.Tidke1/fields1/'
files_ = os.listdir(path)
print('')
print files_
files = []
print('')
print files
for i in range(len(files_)):
    if files_[i][-4:] == '0.nc':
        files = np.append(files, files_[i])
#files = files[:]
print('')
print files
outfile = os.path.join(path, 'to_paraview.nc')
# outfile = '/Users/yaircohen/Documents/PyCLES_out/Output.TRMM_LBA.Tidke1/fields1/to_paraview.nc'
nfiles = len(files)
print nfiles
print('')

t = np.arange(nfiles)

vars = [ 'u', 'v', 'w', 'qt',  'ql', 'temperature', 's'] # qr
units = ['m/s', 'm/s', 'm/s', 'kg/kg', 'kg/kg', 'K', 'K'] # kg/kg
#Get dimension information from one of the files
if files[0][-3:] == '.nc':
    print(os.path.join(path,files[0]))
    ingrp = nc.Dataset(os.path.join(path,files[0]),'r')
else:
    print(os.path.join(path,files[1]))
    ingrp = nc.Dataset(os.path.join(path,files[1]), 'r')
fields = ingrp['fields']
v = fields[vars[0]][:,:,:]
x = np.arange(v.shape[0]) * 200.0 
y = np.arange(v.shape[1]) * 200.0
z = np.arange(v.shape[2]) * 100.0
print('')
rootgrp = nc.Dataset(outfile,'w',format='NETCDF3_64BIT')
xd = rootgrp.createDimension('x', x.shape[0])
yd = rootgrp.createDimension('y', y.shape[0])
zd = rootgrp.createDimension( 'z', z.shape[0])
td = rootgrp.createDimension('Times', nfiles)
print('dims', xd, yd, zd, td)
xv = rootgrp.createVariable('x','f8',('x',))
xv[:] = x
yv = rootgrp.createVariable('y','f8',('y',))
yv[:] = y
zv = rootgrp.createVariable('z','f8',('z',))
zv[:] = z
tv = rootgrp.createVariable('Times','f8',('Times',))
tv.units = "<time length> since <date>"
tv[:] = t

ingrp.close()
rootgrp.sync()
vcount = 0
print('closed')

for v in vars:
    print('sync var', v)
    rootgrp.createVariable(v, 'f8', ('Times', 'z', 'y', 'x',))
    rootgrp[v].units = units[vcount]
    rootgrp.sync()
print('synced')
for v in vars:
    print('var', v)
    count = 0
    for file in files:
        print file
        ingrp = nc.Dataset(os.path.join(path,file))
        fields = ingrp['fields']
        data = fields[v][:,:,:]
        dt = np.zeros((z.shape[0],y.shape[0],x.shape[0]),dtype=np.double)
        #Need to transpose data now
        for i in range(x.shape[0]):
            for j in range(y.shape[0]):
                for k in range(z.shape[0]):
                    dt[k,j,i] = data[i,j,k]

        print np.max(dt)
        rootgrp[v][count,:,:,:] = dt
        ingrp.close()
        count += 1
        rootgrp.sync()
    vcount += 1
rootgrp.close()

'''
var = 'w'

fin = '9000.nc'
ingrp = nc.Dataset(fin,'r')
fields = ingrp['fields']
w = fields[var][:,:,]
#w = np.concatenate([w,w],axis=0)
#w = np.concatenate([w,w],axis=1)
x = np.arange(w.shape[0]) * 50.0
y = np.arange(w.shape[1]) * 50.0
z = np.arange(w.shape[2]) * 50.0
ingrp.close()

#x = np.arange(10, dtype=np.double)
#y = np.arange(10, dtype=np.double)
#z = np.arange(10, dtype=np.double)

data = w

#data[:,:,:6] = 2.0

rootgrp = nc.Dataset(var + '.nc','w')
xd = rootgrp.createDimension('x', x.shape[0])
yd = rootgrp.createDimension('y', y.shape[0])
zd = rootgrp.createDimension( 'z', z.shape[0])

#xv = rootgrp.createVariable('x','f8',('x,'))

xv = rootgrp.createVariable('x','f8',('x',))
xv[:] = x
yv = rootgrp.createVariable('y','f8',('y',))
yv[:] = y
zv = rootgrp.createVariable('z','f8',('z',))
zv[:] = z


datav  = rootgrp.createVariable(var,'f8',('x','y','z'))
datav[:,:,:] = data

rootgrp.close()
'''
