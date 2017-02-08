import netCDF4 as nc
import argparse
import os
import sys
import json as  simplejson
import pylab as plt
from matplotlib.colors import LogNorm
import matplotlib.cm as cm
from cycler import cycler
import pickle

import numpy as np
import itertools
from scipy import linalg
from sklearn import mixture

sys.path.append("..")
from io_read_in_files import read_in_netcdf_fields, read_in_netcdf

label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['legend.fontsize'] = 10

'''
means(z) = [m1(z),m2(z)]                        --> shape = nz x ncomp x nvar
covars(z) = [[c11(z),c12(z)],[c21(z),c22(z)]]   --> shape = nz x ncomp x nvar x nvar
weight(z)                                       --> shape = ncomp
'''


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
    '''
    zrange:     z-values for which the PDF is fitted
    var_list:   list of variables that are included in (multi-variate) PDF
    '''
    global zrange
    # zrange = np.arange(10, 31, 10)
    # zrange = np.asarray([5,10,20,30,50,70,80,100])
    zrange = np.asarray([5, 10, 20, 30, 50, 70])
    print('zrange', zrange*dz)
    print('_______________________')
    if case_name == 'DCBLSoares':
        var_list = ['w', 'u', 's']
    else:
        var_list = ['w', 's', 'qt']
    var_corr_list = ['wtemperatureqt']


    # Test File
    var = var_corr_list[0]
    nc_file_name = 'EM2_trivar_alltimes.nc'
    fullpath_in = os.path.join(in_path, 'EM2_trivar_alltimes', nc_file_name)
    print(fullpath_in)
    means = read_in_netcdf(var, 'means', fullpath_in)
    time_ = read_in_netcdf('t', 'time', fullpath_in)
    print('time', time, time_)
    # global z_max, zrange_, ncomp, nvar
    # zrange_ = read_in_netcdf('height', 'z-profile', fullpath_in)
    # print('zrange from data: ', type(zrange_))
    z_max = means.shape[0]
    ncomp = means.shape[1]
    nvar = means.shape[2]
    print('ncomp, nvar: ', ncomp, nvar)

    # # data = np.ndarray(shape=((nx * ny), nz, nvar))
    data = np.ndarray(shape=((nx * ny), nvar))
    data_all = np.ndarray(shape=(0, nvar))
    '''(1) read in nc   -files variables'''
    fig = plt.figure(figsize=(15,10))
    colmap = cm.get_cmap('viridis')
    for i in range(len(zrange)):
        iz = zrange[i]
        print('i = ' + np.str(iz) + ': ' + np.str(data_all.shape))
        for d in files:
            fullpath_in = os.path.join(args.path, 'fields', d)
            print(fullpath_in)
            var1 = 'w'
            var2 = 'temperature'
            var3 = 'qt'
            data1_ = read_in_netcdf_fields(var1, fullpath_in).reshape((nx * ny, nz))
            data2_ = read_in_netcdf_fields(var2, fullpath_in).reshape((nx * ny, nz))
            data3_ = read_in_netcdf_fields(var3, fullpath_in).reshape((nx * ny, nz))
            print('----', var1, 'T', var3)
            data[:, 0] = data1_[:, iz]
            data[:, 1] = data2_[:, iz]
            data[:, 2] = data3_[:, iz]
            data_all = np.append(data_all, data, axis=0)

        amp_w = 1
        amp_qt = 1e2
        data_aux = np.array(data, copy=True)
        if var1 == 'w':
            print('normalising w')
            data_aux[:, 0] = data[:, 0] * amp_w
        if var2 == 's':
            s_mean = np.average(data[:, 1])
            print('normalising s', s_mean)
            data_aux[:, 1] = data[:, 1] - s_mean
        if var3 == 'qt':
            print('normalising qt')
            data_aux[:, 2] = data[:, 2] * amp_qt

        plt.subplot(1, 3, 1)
        # plt.rc('lines', linewidth=4)
        # plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y']) +
        #                            cycler('linestyle', ['-', '--', ':', '-.'])))
        # plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y'])))
        # plt.rc('axes', prop_cycle=(cycler('color', colmap)))
        plt.scatter(data[:, 0], data[:, 1], s=2, alpha=0.3, color = colmap(0.1*i/z_max), label='z='+str(iz*dz)+'m')
        axis_label(var1, var2, 1, amp_w)
        plt.subplot(1, 3, 2)
        plt.scatter(data[:, 0], data[:, 2], s=2, alpha=0.3, color = colmap(0.1*i/z_max), label='z='+str(iz*dz)+'m')
        axis_label(var1, var3, 1, amp_w)
        plt.subplot(1, 3, 3)
        plt.scatter(data[:, 1], data[:, 2], s=2, alpha=0.3, color = colmap(0.1*i/z_max), label='z='+str(iz*dz)+'m')
        axis_label(var2, var3, 1, amp_w)
        plt.legend(loc=4)

    fig.suptitle('LES Data: ' + var1 + ' ' + var2 + ' ' + var3, fontsize=20)
    plt.savefig(os.path.join(
        fullpath_out, 'EM2_trivar_alltimes_figures',
        'Data_' + var1 + var2+ var3+ '.png')
    )
    # plt.show()


#     '''(2) read in trivar EM2 PDF parameters'''
#     for var in var_corr_list:
#         count_t = 0
#         nc_file_name = 'EM2_trivar_alltimes.nc'
#         fullpath_in = os.path.join(in_path, 'EM2_trivar_alltimes', nc_file_name)
#         print('fullpath_in', fullpath_in)
#         means = read_in_netcdf(var, 'means', fullpath_in)
#         covars = read_in_netcdf(var, 'covariances', fullpath_in)
#         weights = read_in_netcdf(var, 'weights', fullpath_in)
#
#         mean_tot_time[count_t,:,:], covar_tot_time[count_t,:,:,:] = covariance_estimate_from_multicomp_pdf(means, covars, weights)
#         count_t += 1
#
#     '''(2b) sort PDF components according to their mean'''
#         for i1 in range(2):  # loop over variables
#             for n in range(time.size):
#                 for k in range(z_max):
#                     if means_time[n, k, 0, i1] < means_time[n, k, 1, i1]:
#                         aux = means_time[n, k, 1, i1]
#                         means_time[n, k, 1, i1] = means_time[n, k, 0, i1]
#                         means_time[n, k, 0, i1] = aux
#                         for i2 in range(2):
#                             aux = covariance_time[n, k, 1, i1, i2]
#                             covariance_time[n, k, 1, i1, i2] = covariance_time[n, k, 0, i1, i2]
#                             covariance_time[n, k, 0, i1, i2] = aux
#
#     '''(3) plot Data and Histogram'''
#         for i in range(len(zrange)):
#             iz = zrange[i]
#             # plot_PDF_samples(data, var_name1, var_name2, clf, time, z)
#             plot_PDF_samples(var, np.int(d[0:-3]), iz * dz)
#
#
#     '''(4) plot Means and Covariance'''
#     #     print('calling plot covars: ')
#     #     print(zrange_, zrange_.shape, means_time.shape)
#     #     print(time)
#     #     for z0 in range(zrange_.shape[0]):
#     #         # for t0 in [1,6,12]:
#     #         # for t0 in [0,6,12,18]:
#     #         for t0 in [0,1]:
#     #             bivar_plot_means(var, means_time, covariance_time, mean_tot_time, covar_tot_time, z0, t0, time_)
#     #             bivar_plot_covars(var, means_time, covariance_time, mean_tot_time, covar_tot_time, z0, t0, time_)
    return

# #----------------------------------------------------------------------
#----------------------------------------------------------------------
def axis_label(var_name1, var_name2, amp_qt, amp_w):
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    if var_name1 == 'qt':
        plt.xlabel(var_name1 + '  (* ' + np.str(amp_qt) + ')')
        plt.ylabel(var_name2)
    elif var_name1 == 'w':
        plt.xlabel(var_name1 + '  (* ' + np.str(amp_w) + ')')
        plt.ylabel(var_name2)
    if var_name2 == 'qt':
        plt.xlabel(var_name1)
        plt.ylabel(var_name2 + '  (* ' + np.str(amp_qt) + ')')
    elif var_name2 == 'w':
        plt.xlabel(var_name1)
        plt.ylabel(var_name2 + '  (* ' + np.str(amp_w) + ')')

    return


# ----------------------------------------------------------------------
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
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def read_in_netcdf(variable_name, group_name, fullpath_in):
    print('-------', variable_name, group_name, fullpath_in)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups[group_name].variables[variable_name]

    shape = var.shape
    # print('shape:',var.shape)
    data = np.ndarray(shape=var.shape)
    data = var[:]
    rootgrp.close()
    return data


#----------------------------------------------------------------------
#----------------------------------------------------------------------


if __name__ == "__main__":
    main()