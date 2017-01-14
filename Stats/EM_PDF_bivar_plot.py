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

from read_in_files import read_in_netcdf

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
    time = np.zeros(len(files))
    i = 0
    for d in files:
        time[i] = np.int(d[0:-3])
        i+=1
    time = np.sort(time)
    # print('time', time)
    print('')
    # ______________________
    # ______________________
    '''
    zrange:     z-values for which the PDF is fitted
    var_list:   list of variables that are included in (multi-variate) PDF
    '''
    var_list = ['ws','wqt','sqt']

    '''read in nc-files - bivar'''
    count_t = 0

    d = files[0]
    var = 'ws'
    nc_file_name = 'EM2_bivar_' + str(d)
    fullpath_in = os.path.join(in_path, nc_file_name)
    means = read_in_netcdf(var, 'means', fullpath_in)
    time_ = read_in_netcdf('t', 'time', fullpath_in)
    print('time', time, time_)
    global z_max, zrange_, ncomp, nvar
    zrange_ = read_in_netcdf('height', 'z-profile', fullpath_in)
    # zrange = map(int,np.linspace(0, 24, 13))
    print('zrange from data: ', zrange_)
    z_max = means.shape[0]
    ncomp = means.shape[1]
    nvar = 2

    means_time = np.ndarray(shape=(len(files), z_max, ncomp, nvar))
    covariance_time = np.zeros(shape=(len(files), z_max, ncomp, nvar, nvar))
    weights_time = np.zeros(shape=(len(files), z_max, ncomp))
    mean_tot_time = np.ndarray(shape=(len(files), z_max, nvar))
    covar_tot_time = np.zeros(shape=(len(files), z_max, nvar, nvar))

    for var in var_list:
        for d in files:
            nc_file_name = 'EM2_bivar_' + str(d)
            fullpath_in = os.path.join(in_path, nc_file_name)
            print('fullpath_in', fullpath_in)
            means = read_in_netcdf(var, 'means', fullpath_in)
            covars = read_in_netcdf(var, 'covariances', fullpath_in)
            weights = read_in_netcdf(var, 'weights', fullpath_in)

            means_time[count_t,:,:,:] = means[:,:,:]
            covariance_time[count_t, :, :, :] = covars[:, :, :]
            weights_time[count_t,:,:] = weights[:,:]

            mean_tot_time[count_t,:,:], covar_tot_time[count_t,:,:,:] = covariance_estimate_from_multicomp_pdf(means, covars, weights)
            count_t += 1

        print('calling plot covars: ')
        print(zrange_, zrange_.shape, means_time.shape)
        print(time)
        for z0 in range(zrange_.shape[0]):
            # for t0 in [1,6,12]:
            # for t0 in [0,6,12,18]:
            for t0 in [0,1]:
                bivar_plot_means(var, means_time, covariance_time, mean_tot_time, covar_tot_time, z0,t0,time_)
                bivar_plot_covars(var, means_time, covariance_time, mean_tot_time, covar_tot_time, z0,t0,time_)
    return

#----------------------------------------------------------------------
#----------------------------------------------------------------------
def covariance_estimate_from_multicomp_pdf(means, covars, weights):
    '''
    Input: clf - Expectation Maximization (EM) Gaussian Mixture Model (GMM) with weights, and parameters of all PDF components
    Output: parameters of total PDF
    '''
    z_max, ncomp, nvar = means.shape
    mean_tot = np.zeros((z_max, nvar))
    covar_tot = np.zeros((z_max, nvar,nvar))
    for i in range(ncomp):
        for z in range(z_max):
            mean_tot[z,:] += weights[z,i]*means[z,i,:]
            covar_tot[z,:] += weights[z,i]*covars[z,i,:,:]

    return mean_tot, covar_tot
#----------------------------------------------------------------------
def bivar_plot_means(var_name, means_, covars_, mean_tot_, covar_tot_, z0,t0, time_):
    print('plotting means')
    print(os.path.join(fullpath_out, 'figures_EM2_bivar'))
    global time, ncomp, zrange_
    nt = time.size
    colors = ['b', 'g']

    # over time at given level
    means = means_[:,z0,:,:]
    covars = covars_[:,z0,:,:]
    means_tot = mean_tot_[:,z0,:]
    covar_tot = covar_tot_[:, z0, :, :]
    print(time.shape, means.shape)

    fig = plt.figure(figsize=(10,8))
    fig.suptitle(var_name+ r': mean values of $f_1$, $f_2$ (z=' + np.str(zrange_[z0] * dz) + 'm)', fontsize=20)
    for j in range(2):      # loop over variables
        plt.subplot(2,2,j+1)
        for i in range(nt):
            if means[i,0,j] < means[i,1,j]:
                aux = means[i,1,j]
                means[i,1,j] = means[i,0,j]
                means[i, 0, j] = aux
        plt.plot(time[:], means[:, 0, j], 'o-')
        plt.plot(time[:], means[:, 1, j], 'o-')
        plt.plot(time[:], means_tot[:, j], 'ro-')
        # plt.title('means '+var_name[j]+' (z=' + np.str(zrange_[z0] * dz) + 'm)')
        plt.xlabel('time')
        plt.ylabel(r'$<$'+var_name[j]+r'$>$')

        plt.subplot(2, 2, 2+j+1)
        plt.plot(time[:],means[:,0,j],'o',color=colors[0])
        plt.plot(time[:], means[:, 1, j], 'o', color=colors[1])
        plt.plot(time[:], means_tot[:, j], 'ro')
        for comp in range(ncomp):
            for n in range(nt):
                bar = 0.5*np.sqrt(covars[n,comp,j,j])
                plt.plot([time[n],time[n]],[means[n,comp,j]-bar, means[n,comp,j]+bar],color=colors[comp])
                bar = 0.5 * np.sqrt(covar_tot[n, j, j])
                plt.plot([time[n], time[n]], [means_tot[n, j] - bar, means_tot[n, j] + bar], color='r')
        # plt.title('means ' + var_name[j] + ' (z=' + np.str(zrange_[z0] * dz) + 'm)')
        plt.xlabel('time')
        plt.ylabel(r'$<$' + var_name[j] + r'$>$')
    plt.savefig(os.path.join(fullpath_out,'figures_EM2_bivar/') + 'means_time_a_' + var_name + '_z' + str(np.int(zrange_[z0]*dz)) + 'm.png')
    plt.close()

    # over levels, at given time
    means = means_[t0, :, :, :]
    covars = covars_[t0, :, :, :]
    means_tot = mean_tot_[t0, :, :]
    covar_tot = covar_tot_[t0, :, :, :]
    print(means.shape)

    fig = plt.figure(figsize=(10, 8))
    fig.suptitle(var_name + r': mean values of $f_1$, $f_2$ (t=' + np.str(time_[t0]) + ')', fontsize=20)
    for j in range(2):  # loop over variables
        plt.subplot(2, 2, j + 1)
        for k in range(z_max):
            if means[k, 0, j] < means[i, 1, j]:
                aux = means[i, 1, j]
                means[k, 1, j] = means[i, 0, j]
                means[k, 0, j] = aux
        for comp in range(ncomp):
            plt.plot(zrange_[:] * dz, means[:, comp, j], 'o-', color=colors[comp], label='comp ' + np.str(comp))
        plt.plot(zrange_[:] * dz, means_tot[:, j], 'ro-', label='total')
        plt.xlabel('height z')
        plt.ylabel(r'$<$' + var_name[j] + r'$>$')

        plt.subplot(2, 2, 2 + j + 1)
        for comp in range(ncomp):
            for i in range(z_max):
                plt.plot(zrange_[i] * dz, means[i, comp, j], 'o', color=colors[comp], label='comp '+np.str(comp))
                bar = 0.5 * np.sqrt(covars[i, comp, j, j])
                plt.plot([zrange_[i] * dz, zrange_[i] * dz], [means[i, comp, j] - bar, means[i, comp, j] + bar], color=colors[comp])
        plt.xlabel('height z')
        plt.ylabel(r'$<$' + var_name[j] + r'$>$')
    ax = plt.subplot(2, 2, 1)
    ax.legend()
    plt.savefig(
            os.path.join(fullpath_out, 'figures_EM2_bivar/') + 'means_levels_a_' + var_name + '_t' + str(np.int(time_[t0])) + '.png')
    plt.close()

    # plt.figure()
    # for comp in range(ncomp):
    #     for i in range(z_max):
    #         plt.plot(i * dz, means[i, comp, 0], 'o',color=colors[comp])
    #         bar = 0.5 * np.sqrt(covars[i, comp, 0, 0])
    #         plt.plot([i * dz,i * dz],[means[i,comp,0]-bar,means[i,comp,0]+bar],color=colors[comp])
    # plt.title('means (t=' + np.str(time_[t0]) + ')')
    # plt.xlabel('height z')
    # plt.ylabel(var_name)
    # plt.savefig(
    #     os.path.join(fullpath_out, 'figures_EM2_bivar/') + 'means_levels_a_' + var_name + '_t' + str(np.int(time_[t0])) + '.png')
    # plt.figure()
    # for comp in range(ncomp):
    #     for i in range(z_max):
    #         plt.plot(i * dz, means[i, comp, 0], 'o', color=colors[comp])
    # plt.title('means (t=' + np.str(time_[t0]) + ')')
    # plt.xlabel('height z')
    # plt.ylabel(var_name)
    # plt.savefig(
    #     os.path.join(fullpath_out, 'figures_EM2_bivar/') + 'means_levels_b_' + var_name + '_t' + str(np.int(time_[t0])) + '.png')
    return

#----------------------------------------------------------------------
def bivar_plot_covars(var_name, means_, covars_, mean_tot_, covars_tot_, z0,t0,time_):
    print('plotting covars', covars_.shape)
    print(os.path.join(fullpath_out, 'figures_EM2_bivar'))
    global time, ncomp, zrange_
    nt = time.size
    colors = ['b', 'g']
    for i1 in range(2):  # loop over variables
        for j in range(nt):
            for k in range(z_max):
                if means_[j,k, 0, i1] < means_[j,k, 1, i1]:
                    aux = means_[j,k,1,i1]
                    means_[j,k,1,i1] = means_[j,k,0,i1]
                    means_[j,k,0,i1] = aux
                    for i2 in range(2):
                        aux = covars_[j,k,0,i1,i2]
                        covars_[j,k,0,i1,i2] = covars_[j,k,1,i1,i2]
                        covars_[j,k,1,i1,i2] = aux

    # over time at given level
    means = means_[:,z0,:,:]
    covars = covars_[:,z0,:,:,:]
    covars_tot = covars_tot_[:,z0,:,:]
    print(time.shape, means.shape)
    fig = plt.figure(figsize=(12,5))
    fig.suptitle(var_name + r': covariance values of $f_1$, $f_2$ (z=' + np.str(zrange_[z0] * dz) + 'm)', fontsize=20)
    n = 1
    for i in range(2):
        for j in range(i,2):
            plt.subplot(1, 3, n)
            for comp in range(ncomp):
                plt.plot(time[:],covars[:,comp,i,j],'o-', label='comp '+np.str(comp))
                plt.xlabel('time')
                plt.ylabel(var_name[i]+var_name[j])
            plt.plot(time[:], covars_tot[:,i,j],'ro-',label='total')
            if n==1:
                plt.legend()
            n += 1
    plt.savefig(os.path.join(fullpath_out,'figures_EM2_bivar/') + 'covars_time_' + var_name + '_z' + str(np.int(zrange_[z0]*dz)) + 'm.png')
    plt.close()

    # over levels, at given time
    means = means_[t0, :, :, :]
    covars = covars_[t0, :,:, :, :]
    covars_tot = covars_tot_[t0, :, :, :]
    fig = plt.figure(figsize=(12, 5))
    fig.suptitle(var_name + r': covariance values of $f_1$, $f_2$ (t=' + np.str(time[t0]) + ')', fontsize=20)
    n = 1
    for i in range(2):
        for j in range(i, 2):
            plt.subplot(1, 3, n)
            for comp in range(ncomp):
                plt.plot(dz*zrange_, covars[:, comp, i, j], 'o-', color=colors[comp], label='comp '+np.str(comp))
                plt.xlabel('height z')
                plt.ylabel(var_name[i] + var_name[j])
            plt.plot(dz*zrange_[:], covars_tot[:, i, j], 'ro-', label='total')
            if n == 1:
                plt.legend()
            n += 1
    plt.savefig(
        os.path.join(fullpath_out, 'figures_EM2_bivar/') + 'covars_level_' + var_name + '_t' + str(np.int(time_[t0])) + '.png')
    plt.close()
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



#----------------------------------------------------------------------


if __name__ == "__main__":
    main()