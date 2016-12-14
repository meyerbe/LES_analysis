import netCDF4 as nc
import argparse
import os
import json as  simplejson
import pylab as plt
from matplotlib.colors import LogNorm

import numpy as np
import itertools
from scipy import linalg
from sklearn import mixture

from arfit_py import arfit


label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size


'''
Read in: parameters from Gaussian Mixture Model (means, stochastic, relative weights)

Calculation:
linear stochastic model for transition (auto-regression (AR) model)
(a) from one time-step to the next at a given point in space
(b) from one level to the next at a given time step
--> compute transformation matrix A and stochastic noise term

AR model: autoregressive model
v_n = w + (A_1*v_1 + ... + A_{n-1}*v_{n-1}) + eps_n
    v_n: state vector at time t_n
    w: vector of intercept terms (allows non-zero mean)
    A_i: coefficient methods (i=1, ..n-1)
    eps_n(C): noise term (uncorrelated random vectors)

Method for computing AR model:
(1) estimpate model parameters p, A_i, w, C
    (a) model order p: order selection criterion
    (b) A_i, w: least squares method
    (c) C = residual covariance matrix of the estimated model

- eigenmodes of estimated autoregressive (AR) models of first order
'''

def main():
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("path")
    args = parser.parse_args()

    case_name = 'Bomex'
    read_in_nml(args.path, case_name)

    # '''(1a) Read in Univariate'''
    # fullpath_in = os.path.join(in_path,'EM2_univar_3600.nc')
    # print('fullpath_in', fullpath_in)
    # var = 'w'
    # means = read_in_netcdf(var, 'means', fullpath_in)
    # covar = read_in_netcdf(var, 'covariances', fullpath_in)
    # print(means.shape, covar.shape)
    # nz_, ncomp, nvar = means.shape
    # print(nz_, ncomp, nvar)
    #
    # '''(2) AR(1) Model'''
    # # for uni-variate EM2 model: state-vector v = [means1,means2,covar1,covar2]
    # # transition from t1=3600s to t2=7200s at fixed height z0
    #
    # # v = np.ndarray(shape=(ncomp*(1+)))
    # # for i in range(ncomp):
    # i_z = 0
    # i_var = 0
    # # (a) state vector
    # v = [means[i_z,0,i_var],means[i_z,1,i_var],covar[i_z,0,i_var,i_var],covar[i_z,1,i_var,i_var]]

    # -----
    pmin = 1
    pmax = 2
    v = np.ones(shape=(10,2))
    v[:,1] = np.linspace(0,9,10)
    v[:,0] = np.linspace(0, 9, 10)
    print(v)
    arfit(v,pmin,pmax)

    # -----




    '''(3) Eigendecomposition of AR(1)'''
    # v_n = A*v_{n-1} + eps_n


    # global time
    # time = np.zeros((1))
    # # for d in files:
    # #     time[i] = d[0:-3]
    #
    # files = os.listdir(os.path.join(args.path,'fields'))
    # N = len(files)
    # print('Found the following directories', files, N)
    # for d in files:
    #     time = np.sort(np.append(time, np.int(d[0:-3])))
    # # print(time)
    #
    # # (5) IO
    # # (a) create file for eddy fields
    # var_list = ['w','s','qt']
    # for d in files:
    #     time = np.int(d[0:-3])
    #     nc_file_name = 'EM2_' + str(time)
    #     create_statistics_file(fullpath_out, nc_file_name)
    #
    # '''
    # (1) uni-variate PDF for single variable
    # '''
    # data = np.ndarray(shape=((nx*ny),1))
    # zrange= range(0,30,10)
    # # means = np.ndarray(shape=(np.len(varlist)))
    # for var in var_list:
    #     i = 0
    #     for d in files:
    #         t = np.int(d[0:-3])
    #         # print('t summing ' + d)
    #         fullpath_in = os.path.join(args.path, 'fields', d)
    #         print(fullpath_in)
    #         data_ = read_in_netcdf_fields(var,fullpath_in).reshape((nx*ny),nz)
    #         for i in zrange:
    #             data[:,0] = data_[:,i]
    #             means, covariance = Gaussian_mixture_univariate(data, var, t, i*dz)
    #
    #             # # (b) dump eddy fields
    #             # #    add_field(os.path.join(out_path,nc_file_name+'.nc'), var_name)
    #             # dump_variables(os.path.join(out_path, nc_file_name + '.nc'), 'means', u_eddy)
    #             # dump_variables(os.path.join(out_path, nc_file_name + '.nc'), 'v_eddy', v_eddy)
    #             # dump_variables(os.path.join(out_path, nc_file_name + '.nc'), 'w_eddy', w_eddy)
    #             # dump_variables(os.path.join(out_path, nc_file_name + '.nc'), 'phi_eddy', phi_eddy)
    #
    # '''
    # (2) multi-variate PDF for (s,qt,w)
    # '''
    # data = np.ndarray(shape=((nx*ny),2))
    # zrange= range(0,np.int(nz*dz), np.int(5*dz))
    # print('zrange', zrange)
    # zrange = np.asarray(zrange)/dz
    # print('zrange', zrange)
    # for d in files:
    #     fullpath_in = os.path.join(args.path, 'fields', d)
    #     print(fullpath_in)
    #     for var1 in ['w']:
    #         data1_ = read_in_netcdf_fields(var1,fullpath_in).reshape((nx*ny,nz))
    #         for var2 in ['s']:
    #             data2_ = read_in_netcdf_fields(var2, fullpath_in).reshape((nx*ny,nz))
    #             for i in range(0,nz,5):
    #                 data[:,0] = data1_[:,i]
    #                 data[:,1] = data2_[:, i]
    #
    #                 means, covariance = Gaussian_mixture_bivariate(data, var1, var2, np.int(d[0:-3]), i*dz)
    #
    #

    return



# ____________________






def create_statistics_file(path,file_name):
    print('create file:', path)
    rootgrp = nc.Dataset(path+file_name+'.nc', 'w', format='NETCDF4')
    dimgrp = rootgrp.createGroup('dims')
    fieldgrp = rootgrp.createGroup('means')
    fieldgrp = rootgrp.createGroup('covariances')
    fieldgrp.createDimension('n',nx*ny*nz)
    fieldgrp.createDimension('nx', nx)
    fieldgrp.createDimension('ny', ny)
    fieldgrp.createDimension('nz', nz)
    rootgrp.close()
    print('create file end')
    return

def dump_variables(path, var_name, var):
    print('dump variables', path, var_name, var.shape)
    data = np.empty((n[0],n[1],n[2]),dtype=np.double,order='c')
    #        double[:] data = np.empty((Gr.dims.npl,), dtype=np.double, order='c')
    add_field(path, var_name)
    for i in range(0, n[0]):
        for j in range(0, n[1]):
            for k in range(0, n[2]):
                data[i,j,k] = var[i,j,k]
    write_field(path,var_name, data)
    return

def add_field(path, var_name):
    print('add field: ', var_name)
#    rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
    rootgrp = nc.Dataset(path, 'r+')
    fieldgrp = rootgrp.groups['fields']
    #    fieldgrp.createVariable(var_name, 'f8', ('n'))
    var = fieldgrp.createVariable(var_name, 'f8', ('nx', 'ny', 'nz'))
    rootgrp.close()
    return

def write_field(path, var_name, data):
    print('write field:', path, var_name, data.shape)
    rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
    fieldgrp = rootgrp.groups['fields']
    var = fieldgrp.variables[var_name]
    #    var[:] = np.array(data)
    var[:, :, :] = data
    rootgrp.close()
    return


#----------------------------------------------------------------------
def read_in_netcdf_fields(variable_name, fullpath_in):
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups['fields'].variables[variable_name]
    
    shape = var.shape
    # print('shape:',var.shape)
    data = np.ndarray(shape = var.shape)
    data = var[:]
    rootgrp.close()
    return data


# ----------------------------------------------------------------------
def read_in_netcdf(variable_name, group_name, fullpath_in):
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups[group_name].variables[variable_name]

    shape = var.shape
    # print('shape:',var.shape)
    data = np.ndarray(shape=var.shape)
    data = var[:]
    rootgrp.close()
    return data


#----------------------------------------------------------------------
def read_in(variable_name, group_name, fullpath_in):
    f = File(fullpath_in)
    
    #Get access to the profiles group
    group = f[group_name]
    #Get access to the variable dataset
    variable_dataset = group[variable_name]
    #Get the current shape of the dataset
    variable_dataset_shape = variable_dataset.shape

    variable = np.ndarray(shape = variable_dataset_shape)
    for t in range(variable_dataset_shape[0]):
        if group_name == "timeseries":
            variable[t] = variable_dataset[t]
        elif group_name == "profiles":
            variable[t,:] = variable_dataset[t, :]
        elif group_name == "correlations":
            variable[t,:] = variable_dataset[t, :]
        elif group_name == "fields":
            variable[t] = variable_dataset[t]

    f.close()
    return variable
    return

#----------------------------------------------------------------------
def read_in_nml(path, case_name):
    nml = simplejson.loads(open(path + case_name + '.in').read())
    global dz
    dx = nml['grid']['dx']
    dy = nml['grid']['dy']
    dz = nml['grid']['dz']
    global nx, ny, nz, ntot
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    ntot = nx*ny*nz
    print('nx,ny,nz; ntot:', nx, ny, nz, ntot)
    global dt
    # dt = np.int(args.time2) - np.int(args.time1)
    # print('dt', dt)
    global out_path
    out_path = path
    global in_path
    in_path = path


#----------------------------------------------------------------------


if __name__ == "__main__":
    main()