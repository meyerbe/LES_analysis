import netCDF4 as nc
import argparse
import os, sys
import json as  simplejson
import pylab as plt
from matplotlib.colors import LogNorm

import pickle

import numpy as np
import itertools
from scipy import linalg
from sklearn import mixture

sys.path.append("..")
from io_read_in_files import read_in_netcdf

label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['xtick.direction']='out'
plt.rcParams['ytick.direction']='out'
plt.rcParams['legend.fontsize'] = 10

'''
means(z) = [m1(z),m2(z)]                        --> shape = nz x ncomp x nvar
covars(z) = [[c11(z),c12(z)],[c21(z),c22(z)]]   --> shape = nz x ncomp x nvar x nvar
weight(z)                                       --> shape = nz x ncomp
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
    global z_max, zrange, ncomp, nvar, var_list
    if case_name == 'DCBLSoares':
        var_list = ['w', 'u', 's']
    else:
        var_list = ['w', 'temperature', 'qt']
    var_corr_list = ['wtemperatureqt', 'wthetaliqt']

    d = files[0]
    var = var_corr_list[0]
    nc_file_name = 'EM2_trivar_' + str(d)
    fullpath_in = os.path.join(in_path, 'EM2_trivar', nc_file_name)
    means = read_in_netcdf(var, 'means', fullpath_in)
    time_ = read_in_netcdf('t', 'time', fullpath_in)
    print('time', time, time_)
    zrange = read_in_netcdf('height', 'z-profile', fullpath_in)
    print('zrange from data: ', zrange)
    z_max = means.shape[0]
    ncomp = means.shape[1]
    nvar = 3

    means_time = np.ndarray(shape=(len(files), z_max, ncomp, nvar))
    covariance_time = np.zeros(shape=(len(files), z_max, ncomp, nvar, nvar))
    weights_time = np.zeros(shape=(len(files), z_max, ncomp))
    mean_tot_time = np.ndarray(shape=(len(files), z_max, nvar))
    covar_tot_time = np.zeros(shape=(len(files), z_max, nvar, nvar))

    '''(1) read in nc-files - trivar EM2 PDF'''
    for var in var_corr_list:
        if var == 'wthetliqt':
            var_list[2] = 'thetali'
        count_t = 0
        print('')
        print('read in ' + var)
        for d in files:
            nc_file_name = 'EM2_trivar_' + str(d)
            fullpath_in = os.path.join(in_path, 'EM2_trivar', nc_file_name)
            print('fullpath_in', fullpath_in)
            means = read_in_netcdf(var, 'means', fullpath_in)
            covars = read_in_netcdf(var, 'covariances', fullpath_in)
            weights = read_in_netcdf(var, 'weights', fullpath_in)

            means_time[count_t,:,:,:] = means[:,:,:]
            covariance_time[count_t, :, :, :] = covars[:, :, :]
            weights_time[count_t,:,:] = weights[:,:]

            mean_tot_time[count_t,:,:], covar_tot_time[count_t,:,:,:] = covariance_estimate_from_multicomp_pdf(means, covars, weights)
            count_t += 1

        '''(2) sort PDF components according to weights'''
        print('Sorting A')
        for n in range(time.size):
            for k in range(z_max):
                if weights_time[n, k, 0] < weights_time[n, k, 1]:
                    aux = weights_time[n, k, 1]
                    weights_time[n, k, 1] = weights_time[n, k, 0]
                    weights_time[n, k, 0] = aux
                    for i1 in range(nvar):  # loop over variables
                        aux = means_time[n, k, 1, i1]
                        means_time[n, k, 1, i1] = means_time[n, k, 0, i1]
                        means_time[n, k, 0, i1] = aux
                        for i2 in range(nvar):
                            aux = covariance_time[n, k, 1, i1, i2]
                            covariance_time[n, k, 1, i1, i2] = covariance_time[n, k, 0, i1, i2]
                            covariance_time[n, k, 0, i1, i2] = aux
        '''(2a) plot'''
        print('Plotting')
        trivar_plot_means(var, means_time, covariance_time, weights_time, mean_tot_time, covar_tot_time, time_, 'sortA')
        trivar_plot_covars(var, means_time, covariance_time, weights_time, mean_tot_time, covar_tot_time, time_, 'sortA')
        trivar_plot_weights(var, means_time, covariance_time, weights_time, mean_tot_time, covar_tot_time, time_,'sortA')
        print('')

        '''(3) sort PDF components according to their mean'''
        print('Sorting B')
        for n in range(time.size):
            for k in range(z_max):
                if means_time[n, k, 0, 0] < means_time[n, k, 1, 0]:
                    aux = weights_time[n, k, 1]
                    weights_time[n, k, 1] = weights_time[n, k, 0]
                    weights_time[n, k, 0] = aux
                    for i1 in range(nvar):  # loop over variables
                        aux = means_time[n, k, 1, i1]
                        means_time[n, k, 1, i1] = means_time[n, k, 0, i1]
                        means_time[n, k, 0, i1] = aux
                        for i2 in range(nvar):
                            aux = covariance_time[n, k, 1, i1, i2]
                            covariance_time[n, k, 1, i1, i2] = covariance_time[n, k, 0, i1, i2]
                            covariance_time[n, k, 0, i1, i2] = aux

        '''(3a) plot'''
        print('Plotting')
        trivar_plot_means(var, means_time, covariance_time, weights_time, mean_tot_time, covar_tot_time, time_, 'sortB')
        trivar_plot_covars(var, means_time, covariance_time, weights_time, mean_tot_time, covar_tot_time, time_, 'sortB')
        trivar_plot_weights(var, means_time, covariance_time, weights_time, mean_tot_time, covar_tot_time, time_,'sortB')
        print('')

        '''(4) sort PDF components according to Var[w]'''
        print('Sorting C')
        for n in range(time.size):
            for k in range(z_max):
                if covariance_time[n, k, 0, 0, 0] < covariance_time[n, k, 1, 0, 0]:
                    aux = weights_time[n, k, 1]
                    weights_time[n, k, 1] = weights_time[n, k, 0]
                    weights_time[n, k, 0] = aux
                    for i1 in range(nvar):  # loop over variables
                        aux = means_time[n, k, 1, i1]
                        means_time[n, k, 1, i1] = means_time[n, k, 0, i1]
                        means_time[n, k, 0, i1] = aux
                        for i2 in range(nvar):
                            aux = covariance_time[n, k, 1, i1, i2]
                            covariance_time[n, k, 1, i1, i2] = covariance_time[n, k, 0, i1, i2]
                            covariance_time[n, k, 0, i1, i2] = aux
        '''(4a) plot'''
        print('Plotting')
        trivar_plot_means(var, means_time, covariance_time, weights_time, mean_tot_time, covar_tot_time, time_, 'sortC')
        trivar_plot_covars(var, means_time, covariance_time, weights_time, mean_tot_time, covar_tot_time, time_, 'sortC')
        trivar_plot_weights(var, means_time, covariance_time, weights_time, mean_tot_time, covar_tot_time, time_,
                           'sortC')
        print('')

        '''(5) sort PDF components according to E[qt] or Var[qt]'''
        print('Sorting D')
        for n in range(time.size):
            for k in range(z_max):
                if means_time[n, k, 0, 2] < means_time[n, k, 1, 2]:
                    aux = weights_time[n, k, 1]
                    weights_time[n, k, 1] = weights_time[n, k, 0]
                    weights_time[n, k, 0] = aux
                    for i1 in range(nvar):  # loop over variables
                        aux = means_time[n, k, 1, i1]
                        means_time[n, k, 1, i1] = means_time[n, k, 0, i1]
                        means_time[n, k, 0, i1] = aux
                        for i2 in range(nvar):
                            aux = covariance_time[n, k, 1, i1, i2]
                            covariance_time[n, k, 1, i1, i2] = covariance_time[n, k, 0, i1, i2]
                            covariance_time[n, k, 0, i1, i2] = aux

        '''(5a) plot'''
        print('Plotting')
        trivar_plot_means(var, means_time, covariance_time, weights_time, mean_tot_time, covar_tot_time, time_, 'sortB')
        trivar_plot_covars(var, means_time, covariance_time, weights_time, mean_tot_time, covar_tot_time, time_, 'sortB')
        trivar_plot_weights(var, means_time, covariance_time, weights_time, mean_tot_time, covar_tot_time, time_,
                            'sortD')
        print('')

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
def trivar_plot_means(var_name, means_, covars_, weights_, mean_tot_, covar_tot_, time_, sort_type):
    print('plotting means')
    global time, ncomp, zrange, nvar, var_list
    nt = time.size
    colors = ['b', 'g']

    ''' (A) over time at given level'''
    for z0 in range(zrange.shape[0]):
        means = means_[:,z0,:,:]
        covars = covars_[:,z0,:,:]
        means_tot = mean_tot_[:,z0,:]
        covar_tot = covar_tot_[:, z0, :, :]

        fig = plt.figure(figsize=(12, 10))
        fig.suptitle(var_name+ r': mean values of $f_1$, $f_2$ (z=' + np.str(zrange[z0] * dz) + 'm)', fontsize=20)
        for j in range(nvar):      # loop over variables
            plt.subplot(2,3,j+1)
            plt.plot(time[:], means[:, 0, j], 'o-')
            plt.plot(time[:], means[:, 1, j], 'o-')
            plt.plot(time[:], means_tot[:, j], 'ro-')
            plt.xlabel('time')
            plt.ylabel(r'$<$'+var_list[j]+r'$>$')

            plt.subplot(2, 3, 3+j+1)
            plt.plot(time[:],means[:,0,j],'o',color=colors[0])
            plt.plot(time[:], means[:, 1, j], 'o', color=colors[1])
            plt.plot(time[:], means_tot[:, j], 'ro')
            for comp in range(ncomp):
                for n in range(nt):
                    bar = 0.5*np.sqrt(covars[n,comp,j,j])
                    plt.plot([time[n],time[n]],[means[n,comp,j]-bar, means[n,comp,j]+bar],color=colors[comp])
                    bar = 0.5 * np.sqrt(covar_tot[n, j, j])
                    plt.plot([time[n], time[n]], [means_tot[n, j] - bar, means_tot[n, j] + bar], color='r')
            plt.xlabel('time')
            plt.ylabel(r'$<$'+var_list[j]+r'$>$')
        plt.savefig(os.path.join(fullpath_out,'EM2_trivar_figures', sort_type,
                'means_time_' + var_name + '_' + sort_type + '_z' + str(np.int(zrange[z0] * dz)) + 'm.png'))
        plt.close()

    '''over levels, at given time'''
    for t0 in range(time.size):
        means = means_[t0, :, :, :]
        covars = covars_[t0, :, :, :]
        means_tot = mean_tot_[t0, :, :]
        covar_tot = covar_tot_[t0, :, :, :]

        fig = plt.figure(figsize=(12, 10))
        fig.suptitle(var_name + r': mean values of $f_1$, $f_2$ (t=' + np.str(time_[t0]) + ')', fontsize=20)
        for j in range(nvar):  # loop over variables
            plt.subplot(2, 3, j + 1)
            for comp in range(ncomp):
                plt.plot(means[:, comp, j], zrange[:] * dz, 'o-', color=colors[comp], markersize=4, label='comp ' + np.str(comp))
            plt.plot(means_tot[:, j], zrange[:] * dz, 'ro-', markersize=4, label='total')
            plt.xlabel(r'$<$' + var_list[j] + r'$>$')
            plt.ylabel('height z')

            plt.subplot(2, 3, 3 + j + 1)
            for comp in range(ncomp):
                plt.plot(means[:, comp, j], zrange[:] * dz, 'o', markersize=4,
                         color=colors[comp], label='comp '+np.str(comp))
                bar = 0.5 * np.sqrt(covars[:, comp, j, j])
                plt.plot([means[:, comp, j] - bar, means[:, comp, j] + bar], [zrange[:] * dz, zrange[:] * dz], color=colors[comp])
                plt.xlabel(r'$<$' + var_list[j] + r'$>$')
            plt.ylabel('height z')
        ax = plt.subplot(2, 3, 1)
        ax.legend()
        plt.savefig(os.path.join(fullpath_out, 'EM2_trivar_figures', sort_type,
                'means_levels_' + var_name + '_' + sort_type + '_t' + str(np.int(time_[t0])) + '.png'))
        plt.close()

    return
#----------------------------------------------------------------------
def trivar_plot_covars(var_name, means_, covars_, weights_, mean_tot_, covars_tot_, time_, sort_type):
    print('plotting covars', covars_.shape)
    global nvar, time
    global time, ncomp, zrange
    nt = time.size
    colors = ['b', 'g']

    # ''' (A) over time at given level '''
    # for z0 in range(zrange.shape[0]):
    #     means = means_[:,z0,:,:]
    #     covars = covars_[:,z0,:,:,:]
    #     weights = weights_[:, z0, :]
    #     covars_tot = covars_tot_[:,z0,:,:]
    #     fig = plt.figure(figsize=(12,5))
    #     n = 1
    #     for i in range(nvar):
    #         for j in range(i,nvar):
    #             plt.subplot(2, 3, n)
    #             for comp in range(ncomp):
    #                 plt.plot(time[:],covars[:,comp,i,j],'o-', label='comp '+np.str(comp))
    #                 plt.xlabel('time')
    #                 plt.ylabel(var_list[i]+var_list[j])
    #             plt.plot(time[:], covars_tot[:, i, j], 'ro-', label='total')
    #             n += 1
    #     ax = plt.subplot(2,3,1)
    #     ax.legend()
    #     fig.suptitle(var_name + r': covariance values of $f_1$, $f_2$ (z=' + np.str(zrange[z0] * dz) + 'm)', fontsize=20)
    #     plt.savefig(os.path.join(fullpath_out, 'EM2_trivar_figures', sort_type,
    #             'covars_time_' + var_name + '_' + sort_type + '_z' + str(np.int(zrange[z0] * dz)) + 'm.png'))
    #     plt.close()

    ''' (B) over levels, at given time '''
    for t0 in range(time.size):
        means = means_[t0, :, :, :]
        covars = covars_[t0, :,:, :, :]
        weights = weights_[t0, :, :]
        covars_tot = covars_tot_[t0, :, :, :]
        fig = plt.figure(figsize=(12, 10))
        n = 1
        for i in range(nvar):
            for j in range(i, nvar):
                plt.subplot(2, 3, n)
                for comp in range(ncomp):
                    plt.plot(covars[:, comp, i, j], dz * zrange, 'o-', color=colors[comp], markersize=4, label='comp ' + np.str(comp))
                    plt.xlabel(var_list[i] + var_list[j])
                    plt.ylabel('height z')
                plt.plot(covars_tot[:, i, j], dz * zrange[:], 'ro-', markersize=4, label='total')
                n += 1
        ax = plt.subplot(2, 3, 1)
        ax.legend()
        fig.suptitle(var_name + r': covariance values of $f_1$, $f_2$ (t=' + np.str(time[t0]) + ')', fontsize=20)
        plt.savefig(
            os.path.join(fullpath_out, 'EM2_trivar_figures', sort_type,
                'covars_level_' + var_name + '_' + sort_type + '_t' + str(np.int(time_[t0])) + '.png'))
        plt.close()

        fig = plt.figure(figsize=(12, 10))
        n = 1
        for i in range(nvar):
            for j in range(i, nvar):
                plt.subplot(2, 3, n)
                for comp in range(ncomp):
                    plt.plot(weights[:,comp] *covars[:, comp, i, j], dz * zrange, 'o-',
                             markersize=4, color=colors[comp], label='weight*comp ' + np.str(comp))
                    plt.xlabel(var_list[i] + var_list[j])
                    plt.ylabel('height z')
                plt.plot(covars_tot[:, i, j], dz * zrange[:], 'ro-', markersize=4, label='total')
                n += 1
        ax = plt.subplot(2, 3, 1)
        ax.legend()
        fig.suptitle(var_name + r': covariance values of $f_1$, $f_2$ (t=' + np.str(time[t0]) + ')', fontsize=20)
        plt.savefig(
            os.path.join(fullpath_out, 'EM2_trivar_figures', sort_type,
                    'covars_level_alpha_' + var_name + '_' + sort_type + '_t' + str(np.int(time_[t0])) + '.png'))
        plt.close()
    return
#----------------------------------------------------------------------
def trivar_plot_weights(var_name, means_, covars_, weights_, mean_tot_, covar_tot_, time_, sort_type):
    print('plotting weights: ' + os.path.join(fullpath_out, 'figures_EM2_trivar'))
    global time, ncomp, zrange
    nt = time.size
    colors = ['b', 'g']

    ''' (A) over time at given level'''
    for z0 in range(zrange.shape[0]):
        weights = weights_[:, z0, :]

        fig = plt.figure(figsize=(8, 4))
        fig.suptitle(var_name + r': weights of $f_1$, $f_2$ (z=' + np.str(zrange[z0] * dz) + 'm)', fontsize=20)
        plt.plot(time[:], weights[:, 0], 'o-', label='comp 0')
        plt.plot(time[:], weights[:, 1], 'o-', label='comp 1')
        plt.legend()
        plt.title('weights ' + ' (z=' + np.str(zrange[z0] * dz) + 'm)')
        plt.xlabel('time')
        plt.ylabel('weights ' + var_name)

        plt.savefig(os.path.join(
            fullpath_out,'EM2_trivar_figures', sort_type,
                'weights_time_a_' + var_name + '_' + sort_type + '_z' + str(np.int(zrange[z0] * dz)) + 'm.png'))
        plt.close()

    ''' (B) over levels, at given time '''
    for t0 in range(time.size):
        weights = weights_[t0, :, :]
        means = means_[t0, :, :, :]
        covars = covars_[t0, :, :, :]

        fig = plt.figure(figsize=(10, 8))
        fig.suptitle(var_name + r': weights values of $f_1$, $f_2$ (t=' + np.str(time_[t0]) + ')', fontsize=20)
        plt.subplot(1, 3, 1)
        for comp in range(ncomp):
            plt.plot(weights[:, comp], zrange[:] * dz, 'o-', color=colors[comp], label='comp ' + np.str(comp))
        plt.xlabel('weights ' + var_name)
        plt.ylabel('height z')
        plt.subplot(1, 3, 2)
        for comp in range(ncomp):
            plt.plot(means[:, comp, 0], zrange[:] * dz, 'o-', color=colors[comp],
                     label='<' + var_name[0] + '>, comp ' + np.str(comp))
        plt.xlabel('mean ' + var_name[0])
        plt.ylabel('height z')
        plt.legend()
        plt.subplot(1, 3, 3)
        for comp in range(ncomp):
            plt.plot(covars[:, comp, 0, 0], zrange[:] * dz, 'o-', color=colors[comp],
                     label='<' + var_name[0] + var_name[0] + '>, comp ' + np.str(comp))
        plt.legend()
        plt.xlabel('covariance ' + var_name[0] + var_name[0])
        plt.ylabel('height z')

        plt.savefig(
            os.path.join(fullpath_out, 'EM2_trivar_figures', sort_type,
                'weights_levels_' + var_name + '_' + sort_type + '_t' + str(
                np.int(time_[t0])) + '.png'))
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
# # ----------------------------------------------------------------------
# def read_in_netcdf(variable_name, group_name, fullpath_in):
#     print('-------', variable_name, group_name, fullpath_in)
#     rootgrp = nc.Dataset(fullpath_in, 'r')
#     var = rootgrp.groups[group_name].variables[variable_name]
#
#     shape = var.shape
#     # print('shape:',var.shape)
#     data = np.ndarray(shape=var.shape)
#     data = var[:]
#     rootgrp.close()
#     return data
#

#----------------------------------------------------------------------
#----------------------------------------------------------------------


if __name__ == "__main__":
    main()