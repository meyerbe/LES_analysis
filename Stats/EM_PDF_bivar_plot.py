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
    var_list = ['ws','wqt','sqt']
    print(var_list)


    d = files[0]
    var = 'ws'
    nc_file_name = 'EM2_bivar_' + str(d)
    fullpath_in = os.path.join(in_path, 'EM2_bivar', nc_file_name)
    print('')
    print('fullpath_in', fullpath_in, var)
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

    '''(1) read in nc-files - bivar'''
    for var in var_list:
        count_t = 0
        print('')
        print('read in ' + var)
        for d in files:
            nc_file_name = 'EM2_bivar_' + str(d)
            fullpath_in = os.path.join(in_path, 'EM2_bivar', nc_file_name)
            means = read_in_netcdf(var, 'means', fullpath_in)
            covars = read_in_netcdf(var, 'covariances', fullpath_in)
            weights = read_in_netcdf(var, 'weights', fullpath_in)
            print('WEIGHT', weights.shape, len(files), z_max, nvar)

            means_time[count_t,:,:,:] = means[:,:,:]
            covariance_time[count_t, :, :, :] = covars[:, :, :]
            weights_time[count_t,:,:] = weights[:,:]

            mean_tot_time[count_t,:,:], covar_tot_time[count_t,:,:,:] = covariance_estimate_from_multicomp_pdf(means, covars, weights)
            count_t += 1

        '''(2) sort PDF components according to their mean'''
        print('Sorting A')
        for n in range(time.size):
            for k in range(z_max):
                for i1 in range(nvar):  # loop over variables
                    if means_time[n, k, 0, i1] < means_time[n, k, 1, i1]:
                        aux = means_time[n, k, 1, i1]
                        means_time[n, k, 1, i1] = means_time[n, k, 0, i1]
                        means_time[n, k, 0, i1] = aux
                        for i2 in range(nvar):
                            aux = covariance_time[n, k, 1, i1, i2]
                            covariance_time[n, k, 1, i1, i2] = covariance_time[n, k, 0, i1, i2]
                            covariance_time[n, k, 0, i1, i2] = aux
                        aux = weights_time[n, k, 1]
                        weights_time[n, k, 1] = weights_time[n, k, 0]
                        weights_time[n, k, 0] = aux


        '''(2a) plot'''
        print('Plotting')
        # bivar_plot_means(var, means_time, covariance_time, weights_time, mean_tot_time, covar_tot_time, time_, 'sortA')
        # bivar_plot_covars(var, means_time, covariance_time, weights_time, mean_tot_time, covar_tot_time, time_, 'sortA')
        bivar_plot_weights(var, means_time, covariance_time, weights_time, mean_tot_time, covar_tot_time, time_, 'sortA')
        print('')

        '''(3) sort PDF components according to Var[w]'''
        print('Sorting B')
        for n in range(time.size):
            for k in range(z_max):
                if weights_time[n, k, 0] > weights_time[n, k, 1]:
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
        bivar_plot_means(var, means_time, covariance_time, weights_time, mean_tot_time, covar_tot_time, time_, 'sortB')
        bivar_plot_covars(var, means_time, covariance_time, weights_time, mean_tot_time, covar_tot_time, time_, 'sortB')
        bivar_plot_weights(var, means_time, covariance_time, weights_time, mean_tot_time, covar_tot_time, time_,
                           'sortB')
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
def bivar_plot_means(var_name, means_, covars_, weights_, mean_tot_, covar_tot_, time_, sort_type):
    print('plotting means')
    print(os.path.join(fullpath_out, 'figures_EM2_bivar'))
    global time, ncomp, zrange_, nvar
    nt = time.size
    colors = ['b', 'g']

    ''' (A) over time at given level'''
    for z0 in range(zrange_.shape[0]):
        means = means_[:,z0,:,:]
        covars = covars_[:,z0,:,:]
        means_tot = mean_tot_[:,z0,:]
        covar_tot = covar_tot_[:, z0, :, :]

        fig = plt.figure(figsize=(10,8))
        fig.suptitle(var_name+ r': mean values of $f_1$, $f_2$ (z=' + np.str(zrange_[z0] * dz) + 'm)', fontsize=20)
        for j in range(nvar):      # loop over variables
            plt.subplot(2,2,j+1)
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
        plt.savefig(os.path.join(fullpath_out,'EM2_bivar_figures/') + 'means_time_a_' + var_name + '_' + sort_type + '_z' + str(np.int(zrange_[z0]*dz)) + 'm.png')
        plt.close()

    '''over levels, at given time'''
    for t0 in range(time.size):
        means = means_[t0, :, :, :]
        covars = covars_[t0, :, :, :]
        means_tot = mean_tot_[t0, :, :]
        covar_tot = covar_tot_[t0, :, :, :]

        fig = plt.figure(figsize=(10, 8))
        fig.suptitle(var_name + r': mean values of $f_1$, $f_2$ (t=' + np.str(time_[t0]) + ')', fontsize=20)
        for j in range(2):  # loop over variables
            plt.subplot(2, 2, j + 1)
            for comp in range(ncomp):
                plt.plot(means[:, comp, j], zrange_[:] * dz, 'o-', color=colors[comp], label='comp ' + np.str(comp))
            plt.plot(means_tot[:, j], zrange_[:] * dz, 'ro-', label='total')
            plt.xlabel(r'$<$' + var_name[j] + r'$>$')
            plt.ylabel('height z')

            plt.subplot(2, 2, 2 + j + 1)
            for comp in range(ncomp):
                plt.plot(means[:, comp, j], zrange_[:] * dz, 'o', markersize=4,
                         color=colors[comp], label='comp '+np.str(comp))
                bar = 0.5 * np.sqrt(covars[:, comp, j, j])
                plt.plot([means[:, comp, j] - bar, means[:, comp, j] + bar], [zrange_[:] * dz, zrange_[:] * dz], color=colors[comp])
            plt.xlabel(r'$<$' + var_name[j] + r'$>$')
            plt.ylabel('height z')
        ax = plt.subplot(2, 2, 1)
        ax.legend()
        plt.savefig(
                os.path.join(fullpath_out, 'EM2_bivar_figures/') + 'means_levels_' + var_name + '_' + sort_type + '_t' + str(np.int(time_[t0])) + '.png')
        plt.close()

    return
#----------------------------------------------------------------------
def bivar_plot_covars(var_name, means_, covars_, weights_, mean_tot_, covars_tot_, time_, sort_type):
    print('plotting covars', covars_.shape)
    print('')
    global nvar, time
    print(os.path.join(fullpath_out, 'figures_EM2_bivar'))
    global time, ncomp, zrange_
    nt = time.size
    colors = ['b', 'g']

    ''' (A) over time at given level '''
    for z0 in range(zrange_.shape[0]):
        means = means_[:,z0,:,:]
        covars = covars_[:,z0,:,:,:]
        weights = weights_[:, z0, :]
        covars_tot = covars_tot_[:,z0,:,:]
        fig = plt.figure(figsize=(12,5))
        n = 1
        for i in range(nvar):
            for j in range(i,nvar):
                plt.subplot(1, 3, n)
                for comp in range(ncomp):
                    plt.plot(time[:],covars[:,comp,i,j],'o-', label='comp '+np.str(comp))
                    plt.xlabel('time')
                    plt.ylabel(var_name[i]+var_name[j])
                plt.plot(time[:], covars_tot[:, i, j], 'ro-', label='total')
                n += 1
        ax = plt.subplot(1,3,1)
        ax.legend()
        fig.suptitle(var_name + r': covariance values of $f_1$, $f_2$ (z=' + np.str(zrange_[z0] * dz) + 'm)', fontsize=20)
        plt.savefig(os.path.join(fullpath_out,'EM2_bivar_figures/') + 'covars_time_' + var_name + '_' + sort_type + '_z' + str(np.int(zrange_[z0]*dz)) + 'm.png')
        plt.close()

    ''' (B) over levels, at given time '''
    for t0 in range(time.size):
        means = means_[t0, :, :, :]
        covars = covars_[t0, :,:, :, :]
        weights = weights_[t0, :, :]
        covars_tot = covars_tot_[t0, :, :, :]
        fig = plt.figure(figsize=(12, 10))
        n = 1
        for i in range(2):
            for j in range(i, 2):
                plt.subplot(2, 3, n)
                for comp in range(ncomp):
                    plt.plot(covars[:, comp, i, j], dz*zrange_, 'o-', color=colors[comp], label='comp '+np.str(comp))
                    plt.xlabel(var_name[i] + var_name[j])
                    plt.ylabel('height z')
                plt.plot(covars_tot[:, i, j], dz*zrange_[:], 'ro-', label='total')
                plt.subplot(2, 3, n + 3)
                for comp in range(ncomp):
                    plt.plot(weights[:,comp]*covars[:, comp, i, j], dz * zrange_, 'o-', color=colors[comp], label='weight*comp ' + np.str(comp))
                    plt.xlabel(var_name[i] + var_name[j])
                    plt.ylabel('height z')
                plt.plot(covars_tot[:, i, j], dz * zrange_[:], 'ro-', label='total')
                n += 1
        ax = plt.subplot(2, 3, 1)
        ax.legend()
        ax = plt.subplot(2, 3, 4)
        ax.legend()
        fig.suptitle(var_name + r': covariance values of $f_1$, $f_2$ (t=' + np.str(time[t0]) + ')', fontsize=20)
        plt.savefig(
            os.path.join(fullpath_out, 'EM2_bivar_figures/') + 'covars_level_' + var_name + '_' + sort_type + '_t' + str(np.int(time_[t0])) + '.png')
        plt.close()
    return
#----------------------------------------------------------------------
def bivar_plot_weights(var_name, means_, covars_, weights_, mean_tot_, covar_tot_, time_, sort_type):
    print('plotting weights')
    print(os.path.join(fullpath_out, 'figures_EM2_bivar'))
    global time, ncomp, zrange_
    nt = time.size
    colors = ['b', 'g']

    ''' (A) over time at given level'''
    for z0 in range(zrange_.shape[0]):
        weights = weights_[:,z0,:]

        fig = plt.figure(figsize=(10, 5))
        fig.suptitle(var_name + r': weights of $f_1$, $f_2$ (z=' + np.str(zrange_[z0] * dz) + 'm)', fontsize=20)
        plt.plot(time[:], weights[:, 0], 'o-', label='comp '+np.str(comp))
        plt.plot(time[:], weights[:, 1], 'o-', label='comp '+np.str(comp))
        plt.legend()
        plt.title('weights '+' (z=' + np.str(zrange_[z0] * dz) + 'm)')
        plt.xlabel('time')
        plt.ylabel('weights ' + var_name)

        plt.savefig(os.path.join(fullpath_out, 'EM2_bivar_figures/') + 'weights_time_a_' + var_name + '_' + sort_type + '_z' + str(
            np.int(zrange_[z0] * dz)) + 'm.png')
        plt.close()

    ''' (B) over levels, at given time '''
    for t0 in range(time.size):
        weights = weights_[t0, :, :]
        means = means_[t0, :, :, :]
        covars = covars_[t0, :, :, :]
        # means_tot = mean_tot_[t0, :, :]
        # covar_tot = covar_tot_[t0, :, :, :]

        fig = plt.figure(figsize=(10, 8))
        fig.suptitle(var_name + r': weights values of $f_1$, $f_2$ (t=' + np.str(time_[t0]) + ')', fontsize=20)
        plt.subplot(1, 3, 1)
        for comp in range(ncomp):
            plt.plot(weights[:, comp], zrange_[:] * dz, 'o-', color=colors[comp], label='comp ' + np.str(comp))
    #         plt.plot(means_tot[:, j], zrange_[:] * dz, 'ro-', label='total')
        plt.xlabel('weights ' + var_name)
        plt.ylabel('height z')
        plt.subplot(1, 3, 2)
        for comp in range(ncomp):
            plt.plot(means[:, comp, 0], zrange_[:] * dz, 'o-', color=colors[comp], label='<w>, comp ' + np.str(comp))
        plt.xlabel('means ' + var_name)
        plt.ylabel('height z')
        plt.legend()
        plt.subplot(1, 3, 3)
        for comp in range(ncomp):
            plt.plot(covars[:, comp, 0, 0], zrange_[:] * dz, 'o-', color=colors[comp], label='<ww>, comp ' + np.str(comp))
        plt.xlabel('covariances ' + var_name)
        plt.ylabel('height z')

        plt.savefig(
            os.path.join(fullpath_out, 'EM2_bivar_figures/') + 'weights_levels_' + var_name + '_' + sort_type + '_t' + str(
                np.int(time_[t0])) + '.png')
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