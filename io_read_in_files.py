import netCDF4 as nc
import numpy as np
import json as  simplejson
import os

def read_in_nml(path, case_name):
    global nx, dx
    nml = simplejson.loads(open(os.path.join(path,case_name + '.in')).read())
    dx = [nml['grid']['dx'], nml['grid']['dy'], nml['grid']['dz']]
    nx = [nml['grid']['nx'], nml['grid']['ny'], nml['grid']['nz']]
    ntot = nx[0]*nx[1]*nx[2]
    print('dx: ', dx)
    print('nx: ', nx)
    dt = nml['fields_io']['frequency']
    print('dt: ', dt)

    # global out_path
    # out_path = path
    # global in_path
    # in_path = path
    return dx, nx, dt


# ----------------------------------------------------------------------
def read_in_netcdf(variable_name, group_name, fullpath_in):
    print('io_read_in_files: read in netcdf', variable_name, group_name, fullpath_in)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    # grp = rootgrp.groups[group_name]
    if group_name == 'reference':
        var = rootgrp.groups[group_name].variables[variable_name][:]
        rootgrp.close()
    elif group_name == 'profiles':
        var = rootgrp.groups[group_name].variables[variable_name][:,:]
        rootgrp.close()
    elif group_name == 'fields':
        var = rootgrp.groups[group_name].variables[variable_name][:,:,:]
        rootgrp.close()
    return var
    # var = grp.variables[variable_name]
    # # var = rootgrp.groups[group_name].variables[variable_name]
    # data = np.ndarray(shape=var.shape)
    # data = var[:]
    # rootgrp.close()
    # return data


def read_in_netcdf_fields(variable_name, fullpath_in):
    print('io_read_in_files: read in netcdf fields', variable_name, fullpath_in)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    # data = np.ndarray(shape = var.shape)
    # data[:,:,:] = var[:,:,:]
    var = rootgrp.groups['fields'].variables[variable_name][:,:,:]
    rootgrp.close()
    return var


# def test():
#     rootgrp = Dataset('test.nc','w',format='NETCDF4')
#
#     '''create Groups & Subgroups'''
#     # (1) create Groups
#     grp_a = rootgrp.createGroup('test_data')
#     # (2) create Sub-Groups
#     grp_a_a = grp_a.createGroup('test_data_group')
#     print(rootgrp.groups)
#     print(grp_a.groups)
#     print(grp_a.groups.values)
#
#     ''' create Variables '''
#     # (1) define dimensions
#     time = grp_a_a.createDimension('time', size =None)
#     print('Dimensions:', grp_a_a.dimensions)
#     grp_a_a.createDimension('lat', 73)
#     print('Dimensions:', grp_a_a.dimensions)
#     print('Time dimension:', 'time'.__len__())
#     print('Time', len('time'))
#     print('Time', len(time))
#     print('Lat dimension:', 'lat'.__len__())
#     dim = len(grp_a_a.dimensions['time'])
#     print('dd', dim)
#     print(grp_a_a.dimensions.values())
#
#     # (2) create Variable
#     times = grp_a_a.createVariable('time','f8',('time',))
#     latitudes = grp_a_a.createVariable('latitude','f4',('lat',))
#
#     rootgrp.close()
#
#     # ___________
#
#
#
#     f = Dataset('100/0.nc')
#     for grp in walktree(f):
#         for a in grp:
#             print(a)
#             # print(f.groups)
#
#
# def walktree(top):
#     values = top.groups.values()
#     yield values
#     for value in top.groups.values():
#         for children in walktree(value):
#             yield children



