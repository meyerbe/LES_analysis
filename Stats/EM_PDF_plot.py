import netCDF4 as nc
import argparse
import os
import json as  simplejson
import pylab as plt
from matplotlib.colors import LogNorm

import pickle

import numpy as np
import itertools
from scipy import linalg
from sklearn import mixture


label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size




def main():
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("path")
    parser.add_argument("casename")
    args = parser.parse_args()
    # ______________________
    case_name = args.casename
    read_in_nml(args.path, case_name)
    print('nx,ny,nz; ntot:', nx, ny, nz, ntot)
    global fullpath_out
    fullpath_out = args.path
    print('fullpath_out:', fullpath_out)
    # ______________________
    files = os.listdir(os.path.join(args.path,'fields'))
    N = len(files)
    print('Found the following directories', files, N)
    # ______________________
    global time
    time = np.zeros((1))
    for d in files:
        time = np.sort(np.append(time, np.int(d[0:-3])))
    # ______________________
    # ______________________
    '''
    zrange:     z-values for which the PDF is fitted
    var_list:   list of variables that are included in (multi-variate) PDF
    '''
    # zrange = map(int,np.linspace(0, 24, 13))
    # print('zrange', zrange)
    var_list = ['w','s','qt']
    # var_list = ['w','s']
    # var_list = ['s']


    '''read in nc-files'''
    for var in var_list:
        for d in files:
            nc_file_name = 'EM2_univar_' + str(d)
            fullpath_in = os.path.join(in_path, nc_file_name)
            print('fullpath_in', fullpath_in)
            means = read_in_netcdf(var, 'means', fullpath_in)
            covar = read_in_netcdf(var, 'covariances', fullpath_in)
            weights = read_in_netcdf(var, 'weights', fullpath_in)

            if var == 'w':
                max = 0.8
                min = -max
            elif var == 's':
                min = 6950
                max = 6970
            elif var == 'qt':
                max = 1.0
                min = -max
            plot_PDFs_levels(var, means, covar, weights, d[0:-3], min, max)
            print('...............varvarvar...............', var)
            print('')




    return



#----------------------------------------------------------------------
#----------------------------------------------------------------------
def plot_PDFs_levels(var_name, means,covars,weights,t,min,max):
    import matplotlib.mlab as mlab
    import matplotlib.cm as cm

    cmap1 = cm.get_cmap('winter')
    cmap2 = cm.get_cmap('spring')
    cmap3 = cm.get_cmap('gray')
    cmap4 = cm.get_cmap('jet')

    global dz
    print('plot PDFs')

    # means = data_['means']
    # covar = data_['covars']
    # weights = data_['weights']

    plt.figure(figsize=(12,12))
    x = np.linspace(min, max, 300)
    # mu = 0
    # sigma = 1
    # plt.plot(x,mlab.normpdf(x,mu,sigma))
    nz = means.shape[0]
    nz = 10
    nzi = 1. / nz
    plt.subplot(2,2,1)
    for i in range(nz):
        print('i',i,i/nz)
        # color = cm.jet(i)
        plt.plot(x, mlab.normpdf(x, means[i,0,0], covars[i,0,0,0]),color=cmap4(float(i)*nzi),linewidth=2,label='z='+str(i*dz))
    plt.legend(fontsize=10)
    plt.title(var_name + ': PDF 1')

    plt.subplot(2,2, 2)
    for i in range(nz):
        print('i', i, i / nz)
        plt.plot(x, mlab.normpdf(x, means[i, 1, 0], covars[i, 1, 0, 0]), '-', color=cmap4(float(i)*nzi), linewidth=2,label='z='+str(i*dz))
    plt.legend(fontsize=10)
    plt.title(var_name+': PDF 2')

    plt.subplot(2,2, 3)
    for i in range(nz):
        print('i', i, i / nz)
        # color = cm.jet(i)
        plt.plot(x, mlab.normpdf(x, means[i, 0, 0], covars[i, 0, 0, 0]), color=cmap4(float(i)*nzi), linewidth=2,
                 label='z=' + str(i*dz))
        plt.plot(x, mlab.normpdf(x, means[i, 1, 0], covars[i, 1, 0, 0]), '--', color=cmap4(float(i)*nzi), linewidth=2)
    plt.legend(fontsize=10)
    plt.title(var_name+': PDF 1+2')

    plt.subplot(2,2, 4)
    for i in range(nz):
        print('i', i, i / nz)
        plt.plot(x, mlab.normpdf(x, means[i, 0, 0], covars[i, 0, 0, 0]) + mlab.normpdf(x, means[i, 1, 0],
                                                                                       covars[i, 1, 0, 0]),
                 color=cmap4(float(i)*nzi), linewidth=2,
                 label='z=' + str(i*dz))
    plt.legend(fontsize=10)
    plt.title(var_name+': PDF 1+2')

    plt.savefig('../figures_EM/EM2_PDF_univar_levels_' + var_name + '_' + str(t) + '.png')
    plt.close()
    return


# ----------------------------------------------------------------------
def read_in_netcdf(variable_name, group_name, fullpath_in):
    print('read in netcdf', variable_name, group_name)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups[group_name].variables[variable_name]

    shape = var.shape
    # print('shape:',var.shape)
    data = np.ndarray(shape=var.shape)
    data = var[:]
    rootgrp.close()
    return data

#----------------------------------------------------------------------
# def read_in_netcdf_fields(variable_name, fullpath_in):
#     rootgrp = nc.Dataset(fullpath_in, 'r')
#     var = rootgrp.groups['fields'].variables[variable_name]
#
#     shape = var.shape
#     # print('shape:',var.shape)
#     data = np.ndarray(shape = var.shape)
#     data = var[:]
#     rootgrp.close()
#     return data

#----------------------------------------------------------------------
# def read_in(variable_name, group_name, fullpath_in):
#     f = File(fullpath_in)
#
#     #Get access to the profiles group
#     profiles_group = f[group_name]
#     #Get access to the variable dataset
#     variable_dataset = profiles_group[variable_name]
#     #Get the current shape of the dataset
#     variable_dataset_shape = variable_dataset.shape
#
#     variable = np.ndarray(shape = variable_dataset_shape)
#     for t in range(variable_dataset_shape[0]):
#         if group_name == "timeseries":
#             variable[t] = variable_dataset[t]
#         elif group_name == "profiles":
#             variable[t,:] = variable_dataset[t, :]
#         elif group_name == "correlations":
#             variable[t,:] = variable_dataset[t, :]
#         elif group_name == "fields":
#             variable[t] = variable_dataset[t]
#
#     f.close()
#     return variable

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
    dz = nml['grid']['dz']
    print('dz: ', dz)
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