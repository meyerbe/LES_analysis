import netCDF4 as nc
import argparse
import os
import json as  simplejson
import pylab as plt
from matplotlib.colors import LogNorm
import pickle


import numpy as np
import itertools
# from scipy import linalg
# from sklearn import mixture
from Stats.VAR_LES import var_modeling
import statsmodels as sm

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
    # ______________________
    case_name = 'Bomex'
    read_in_nml(args.path, case_name)
    # ______________________
    in_path = args.path
    files = os.listdir(os.path.join(args.path, 'fields'))
    N = len(files)
    print('Found the following directories', files, N)
    # ______________________


    # ------------------------------------------------------------------------------------------
    '''Using Example Data'''
    mdata = sm.datasets.macrodata.load().data

    # ------------------------------------------------------------------------------------------



    # ------------------------------------------------------------------------------------------
    '''
    Using statistics:
    means = nz x ncomp x nvar
    covariance = nz x ncomp x nvar x nvar
    nobs:     number of observations (time steps or levels)
    nphi:     number of parameters of PDFs
    '''
    # ------------------------------------------------------------------------------------------


    '''(1) Read in EM2, Univariate'''
    # ----------------------
    print '............'
    # get shape
    nc_file_name = 'EM2_univar_' + str(files[0])
    print 'nc file-name: ' + str(nc_file_name)
    means = read_in_netcdf('w', 'means', os.path.join(in_path, nc_file_name))
    nz_, ncomp, nvar = means.shape
    print('nz, ncomp, nvar', nz_, ncomp, nvar)
    print '............'
    # ----------------------

    '''(a) using time-steps as observables'''
    # do separately for each variable
    nobs = len(files)  # here number of timesteps
    nphi = 2*ncomp      # number of parameters for univariate PDF: nphi=2*ncomp
    # nphi = means.shape[1] + covar.shape[1]  # for univariate PDF nvar=2*ncomp
    arr = np.ndarray(shape=(nobs, nphi))
    phi = np.zeros(shape=(len(files),2*nvar*ncomp))
    # print('phi', phi.shape)

    '''read in nc-files'''
    var = 's'
    z0 = 0
    count = 0
    for d in files:
        nc_file_name = 'EM2_univar_' + str(d)
        fullpath_in = os.path.join(in_path, nc_file_name)
        print('fullpath_in', fullpath_in)
        means = read_in_netcdf(var, 'means', fullpath_in)
        covar = read_in_netcdf(var, 'covariances', fullpath_in)
        phi1 = means[z0, :].reshape(means[z0, :].size)
        phi2 = covar[z0, :, :, :].reshape(covar[z0, :, :, :].size)
        phi[count,:] = np.concatenate((phi1, phi2))
        # phi = np.concatenate((means[z0, :].reshape(means[z0, :].size), covar[z0, :, :, :].reshape(covar[z0, :, :, :].size)))
        # print means.shape
        # print covar.shape
        # print phi.shape
        print(means.shape, covar.shape)
        print(phi1.shape, phi2.shape)
        print(phi.shape)
        print '............'

        arr[count, 0:2] = means[z0, :, 0]
        arr[count, 2:4] = covar[z0, :, 0,0]
        count += 1
    print arr.shape, type(arr)
    # var_modeling(arr)


    '''read in pickles-files'''
    count = 0
    arr = np.ndarray(shape=(nobs,nphi))
    for d in files:
        pkl_file_name = 'EM2_univar_' + str(d[0:-3]) + '.pkl'
        fullpath_in = os.path.join(in_path, pkl_file_name)
        print('fullpath_in', fullpath_in)
        print '............'
        data_var = read_in_pickle(var, fullpath_in)
        means = data_var['means']
        covar = data_var['covars']
        # means = read_in_pickle(var, 'means', fullpath_in)
        # covar = read_in_pickle(var, 'covars', fullpath_in)
        phi1 = means[z0, :].reshape(means[z0, :].size)
        phi2 = covar[z0, :, :, :].reshape(covar[z0, :, :, :].size)
        phi[count, :] = np.concatenate((phi1, phi2))
        # # phi = np.concatenate((means[z0, :].reshape(means[z0, :].size), covar[z0, :, :, :].reshape(covar[z0, :, :, :].size)))
        # print covar.shape
        # print phi.shape
        # print phi
        # print(arr.shape)
        # print(arr[count, 0:2].shape)
        # print means.shape
        # print(means[z0, :, 0].shape)
        arr[count, 0:2] = means[z0, :, 0]
        arr[count, 2:4] = covar[z0, :, 0, 0]
        count += 1

    # var_modeling(arr)
    print '............'

    # for d in files:
    #     pkl_file_name = 'EM2_univar_' + str(d[0:-3]) + '.pkl'
    #     fullpath_in = os.path.join(in_path, pkl_file_name)
    #     print 'xxxxxxxxxxx'
    #     # statsmodels.iolib.smpickle.load_pickle('test/EM2_univar_7200.pkl')
    #     f = sm.iolib.smpickle.load_pickle(fullpath_in)
    #     # df = sm.datasets.get_rdataset("Guerry", "HistData").data
    #     print f
    # print 'xxxxxxxxxxx'



    '''(b) using levels as observables'''
    nobs = means.shape[0]                   # here number of z
    nphi = means.shape[1] + covar.shape[1]  # for univariate PDF nvar=2*ncomp
    arr = np.ndarray(shape=(nobs, nphi))
    phi = np.zeros(shape=(len(files), 2 * nvar * ncomp))
    print('phi', phi.shape)

    '''read in pickles-files'''
    count = 0
    for d in files:
        pkl_file_name = 'EM2_univar_' + str(d[0:-3]) + '.pkl'
        fullpath_in = os.path.join(in_path, pkl_file_name)
        print('fullpath_in', fullpath_in)
        print
        '............'
        data_var = read_in_pickle(var, fullpath_in)
        means = data_var['means']
        covar = data_var['covars']
        print('shape', means.shape)

        # do separately for each variable
        arr[:, 0:2] = means[:, :, 0]
        arr[:, 2:4] = covar[:, :, 0, 0]
        print('arr: ', arr.shape)

        var_modeling(arr,2)
    print '............'



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

    # '''(3) Eigendecomposition of AR(1)'''
    # v_n = A*v_{n-1} + eps_n



    # ------------------------------
    # ------------------------------
    pmin = 1
    pmax = 2
    v = np.ones(shape=(10,2))
    v[:,1] = np.linspace(0,9,10)
    v[:,0] = np.linspace(0, 9, 10)
    # print(v)
    # arfit(v,pmin,pmax)

    # ------------------------------
    # ------------------------------


    return


# ____________________
def dump_pickle(data,out_path,file_name):
    data_ = (1.4,42)
    # output = open(os.path.join(out_path,'data.pkl'), 'w')
    output = open(os.path.join(out_path, file_name), 'w')
    pickle.dump(data, output)
    output.close()
    return

# def read_in_pickle(var_name, type, fullpath_in):
def read_in_pickle(var_name, fullpath_in):
    print ''
    print '------- read in pickle: '+var_name+' ------'
    # fullpath_in = os.path.join(in_path,file_name)
    f = open(fullpath_in)
    data = pickle.load(f)
    # print(data)
    # print ''
    var = data[var_name]
    # print(var)
    # var_data = var[type]
    # print(var_data)
    print '-------------------------'
    print ''
    # return var_data
    return var
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
    print '------- read in nc-file: ' + group_name + ', ' + variable_name + ' ------'
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
    nml = simplejson.loads(open(os.path.join(path,case_name + '.in')).read())
    global dx, dy, dz
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