import netCDF4 as nc
import numpy as np
import json as  simplejson
import os

def read_in_nml(path, case_name):
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
def read_in_netcdf_fields(variable_name, fullpath_in):
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups['fields'].variables[variable_name]

    # print('shape:',var.shape)
    data = np.ndarray(shape=var.shape)
    data = var[:]
    rootgrp.close()
    return data

# ----------------------------------------------------------------------
def read_in_netcdf(variable_name, group_name, fullpath_in):
    print('read in netcdf', variable_name, group_name)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    grp = rootgrp.groups[group_name]
    var = grp.variables[variable_name]
    # var = rootgrp.groups[group_name].variables[variable_name]

    shape = var.shape
    # print('shape:',var.shape)
    data = np.ndarray(shape=var.shape)
    data = var[:]
    rootgrp.close()
    return data