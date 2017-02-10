import netCDF4 as nc
import argparse
import os, sys
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
plt.rcParams['lines.linewidth'] = 2
# plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y']) +
#                            cycler('linestyle', ['-', '--', ':', '-.'])))
# plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y'])))
# plt.rc('axes', prop_cycle=(cycler('color', colmap)))

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
    global path, fullpath_in
    path = args.path
    print('path:', path)
    fullpath_in = os.path.join(in_path, 'EM2_trivar_alltimes', 'EM2_trivar_alltimes.nc')
    # ______________________
    files = os.listdir(os.path.join(args.path,'fields'))
    N = len(files)
    # print('Found the following directories', files, N)
    # ______________________
    '''
        zrange:     z-values for which the PDF is fitted
        var_list:   list of variables that are included in (multi-variate) PDF
        amp:        list of normalisation factors for data
    '''
    print('')
    global time
    # time = np.zeros((1))
    # for d in files:
    #     time = np.sort(np.append(time, np.int(d[0:-3])))
    time = read_in_netcdf('t', 'time', fullpath_in)
    print('time', time)
    # ______________________
    print('')
    global zrange, var_list
    # zrange = np.asarray([5, 10, 20, 30, 50, 70])
    zrange = read_in_netcdf('height', 'z-profile', fullpath_in)
    print('zrange from data: ', zrange * dz)
    # ______________________

    if case_name == 'DCBLSoares':
        var_list = ['w', 'u', 's']
    else:
        var_list = ['w', 'temperature', 'qt']
    var_corr_list = ['wtemperatureqt', 'wthetaliqt']



    # Test File
    print('')
    global z_max, ncomp, nvar, amp
    means = read_in_netcdf(var_corr_list[0], 'means', fullpath_in)
    z_max = means.shape[0]
    ncomp = means.shape[1]
    nvar = means.shape[2]
    amp = np.ones(nvar)
    print('ncomp, nvar, amp: ', ncomp, nvar, amp)

    # data = np.ndarray(shape=((nx * ny), nz, nvar))
    # data = np.ndarray(shape=((nx * ny), nvar))
    # data_all = np.ndarray(shape=(0, nvar))
    '''(1) read in 3d variable fields: print data'''
    for var in var_corr_list:
        print('var', var)
        if var == 'wthetaliqt':
            var_list[1] = 'thetali'

        plot_data(files)


    '''(2) read in trivar EM2 PDF parameters'''
    nc_file_name = 'EM2_trivar_alltimes.nc'
    for var in var_corr_list:
        if var == 'wthetaliqt':
            var_list[1] = 'thetali'
        elif var == 'wtemperatureqt':
            var_list[1] = 'temperature'
        # count_t = 0
        means = read_in_netcdf(var, 'means', fullpath_in)
        covars = read_in_netcdf(var, 'covariances', fullpath_in)
        weights = read_in_netcdf(var, 'weights', fullpath_in)
        mean_tot, covars_tot = covariance_estimate_from_multicomp_pdf(means,covars,weights)
        print('')


        '''(3) sort PDF components according to their weight'''
        print('Sorting A: weights')
        for k in range(z_max):
            if weights[k, 0] < weights[k, 1]:
                aux = weights[k, 1]
                weights[k, 1] = weights[k, 0]
                weights[k, 0] = aux
                for i1 in range(nvar):  # loop over variables
                    aux = means[k, 1, i1]
                    means[k, 1, i1] = means[k, 0, i1]
                    means[k, 0, i1] = aux
                    for i2 in range(nvar):
                        aux = covars[k, 1, i1, i2]
                        covars[k, 1, i1, i2] = covars[k, 0, i1, i2]
                        covars[k, 0, i1, i2] = aux
        '''(3a) plot'''
        print('Plotting')
        trivar_plot_means(var, weights, means, covars, mean_tot, covars_tot, time, 'sortA')
        trivar_plot_covars(var, weights, means, covars, mean_tot, covars_tot, time, 'sortA')
        trivar_plot_weights(var, weights, means, covars, mean_tot, covars_tot, time, 'sortA')
        print('')


        '''(4) sort PDF components according to <qt>'''
        print('Sorting B: E[qt]')
        for k in range(z_max):
            if means[k, 0, 2] < means[k, 1, 2]:
                aux = weights[k, 1]
                weights[k, 1] = weights[k, 0]
                weights[k, 0] = aux
                for i1 in range(nvar):  # loop over variables
                    aux = means[k, 1, i1]
                    means[k, 1, i1] = means[k, 0, i1]
                    means[k, 0, i1] = aux
                    for i2 in range(nvar):
                        aux = covars[k, 1, i1, i2]
                        covars[k, 1, i1, i2] = covars[k, 0, i1, i2]
                        covars[k, 0, i1, i2] = aux
        '''(4a) plot'''
        print('Plotting')
        trivar_plot_means(var, weights, means, covars, mean_tot, covars_tot, time, 'sortB')
        trivar_plot_covars(var, weights, means, covars, mean_tot, covars_tot, time, 'sortB')
        trivar_plot_weights(var, weights, means, covars, mean_tot, covars_tot, time, 'sortB')
        print('')


        '''(5) sort PDF components according to <w>'''
        print('Sorting C: E[w]')
        for k in range(z_max):
            if means[k, 0, 0] < means[k, 1, 0]:
                aux = weights[k, 1]
                weights[k, 1] = weights[k, 0]
                weights[k, 0] = aux
                for i1 in range(nvar):  # loop over variables
                    aux = means[k, 1, i1]
                    means[k, 1, i1] = means[k, 0, i1]
                    means[k, 0, i1] = aux
                    for i2 in range(nvar):
                        aux = covars[k, 1, i1, i2]
                        covars[k, 1, i1, i2] = covars[k, 0, i1, i2]
                        covars[k, 0, i1, i2] = aux
        '''(5a) plot'''
        print('Plotting')
        trivar_plot_means(var, weights, means, covars, mean_tot, covars_tot, time, 'sortC')
        trivar_plot_covars(var, weights, means, covars, mean_tot, covars_tot, time, 'sortC')
        trivar_plot_weights(var, weights, means, covars, mean_tot, covars_tot, time, 'sortC')
        print('')

    return

#----------------------------------------------------------------------
def covariance_estimate_from_multicomp_pdf(means, covars, weights):
    print('')
    print('Computing total Mean & Covariance')
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
# ----------------------------------------------------------------------
def plot_data(files):
    global path, zrange, var_list, amp
    fig = plt.figure(figsize=(15,10))
    colmap = cm.get_cmap('viridis')
    print('')
    print('plotting data: '+str(var_list))
    for d in files:
        fullpath_in = os.path.join(path, 'fields', d)
        data1_ = read_in_netcdf_fields(var_list[0], fullpath_in).reshape((nx * ny, nz))
        data2_ = read_in_netcdf_fields(var_list[1], fullpath_in).reshape((nx * ny, nz))
        data3_ = read_in_netcdf_fields(var_list[2], fullpath_in).reshape((nx * ny, nz))

        for i in range(len(zrange)):
            iz = np.int(zrange[i])
            # print('i = ' + np.str(iz) + ': ' + str(zrange))
            # plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y']) +
            #                            cycler('linestyle', ['-', '--', ':', '-.'])))
            # plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y'])))
            # plt.rc('axes', prop_cycle=(cycler('color', colmap)))
            plt.subplot(1, 3, 1)
            plt.scatter(data1_[:, iz], data2_[:, iz], s=2, alpha=0.3, color = colmap(0.5*i/z_max), label='z='+str(iz*dz)+'m')
            axis_label(var_list[0], var_list[1], amp[0], amp[1])
            plt.subplot(1, 3, 2)
            plt.scatter(data1_[:, iz], data3_[:, iz], s=2, alpha=0.3, color = colmap(0.5*i/z_max), label='z='+str(iz*dz)+'m')
            axis_label(var_list[0], var_list[2], amp[0], amp[2])
            plt.subplot(1, 3, 3)
            plt.scatter(data2_[:, iz], data3_[:, iz], s=2, alpha=0.3, color = colmap(0.5*i/z_max), label='z='+str(iz*dz)+'m')
            axis_label(var_list[1], var_list[2], amp[1], amp[2])
        if d == files[0]:
            if var_list[1] == 'thetali':
                plt.legend(loc=1)
            else:
                plt.legend(loc=4)

    fig.suptitle('LES Data: ' + var_list[0] + ' ' + var_list[1] + ' ' + var_list[2] , fontsize=20)
    plt.savefig(os.path.join(
        path, 'EM2_trivar_alltimes_figures',
        'Data_' + var_list[0] + var_list[1] + var_list[2] + '.png')
    )
    return
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
def trivar_plot_means(var_name, weights, means, covars, means_tot, covar_tot, time_, sort_type):
    print('plotting means')
    global time, ncomp, zrange, nvar, var_list, path
    nt = time.size
    colmap = cm.get_cmap('viridis')
    colors = ['b', 'g']

    '''over levels'''
    fig = plt.figure(figsize=(12, 10))
    fig.suptitle(var_name + r': mean values of $f_1$, $f_2$', fontsize=20)
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
    plt.savefig(os.path.join(path, 'EM2_trivar_alltimes_figures', sort_type,
                             'means_levels_' + var_name + '_' + sort_type + '.png'))
    plt.close()

    return
#----------------------------------------------------------------------
def trivar_plot_covars(var_name, weights, means, covars, mean_tot, covars_tot, time_, sort_type):
    global nvar, time
    global time, ncomp, zrange, path
    print('plotting covars:', var_name, var_list)
    nt = time.size
    colmap = cm.get_cmap('viridis')
    colors = [colmap(0), colmap(10)]
    colors = ['b', 'g']

    ''' over levels '''
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
    fig.suptitle(var_name + r': covariance values of $f_1$, $f_2$', fontsize=20)
    plt.savefig(
        os.path.join(path, 'EM2_trivar_alltimes_figures', sort_type,
                     'covars_level_' + var_name + '_' + sort_type + '.png'))
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
    fig.suptitle(var_name + r': covariance values of $f_1$, $f_2$', fontsize=20)
    plt.savefig(
        os.path.join(path, 'EM2_trivar_alltimes_figures', sort_type,
                     'covars_level_alpha_' + var_name + '_' + sort_type + '.png'))
    plt.close()
    return
#----------------------------------------------------------------------
def trivar_plot_weights(var_name, weights, means, covars, mean_tot, covar_tot, time_, sort_type):
    print('plotting weights')
    global time, ncomp, zrange, path
    nt = time.size
    colmap = cm.get_cmap('viridis')
    colors = [colmap(0), colmap(.5)]
    colors = ['b', 'g']

    ''' (B) over levels, at given time '''
    fig = plt.figure(figsize=(10, 8))
    fig.suptitle(var_name + r': weights values of $f_1$, $f_2$', fontsize=20)
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
        os.path.join(path, 'EM2_trivar_alltimes_figures', sort_type,
                     'weights_levels_' + var_name + '_' + sort_type + '.png'))
    plt.close()
    return
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



#----------------------------------------------------------------------
#----------------------------------------------------------------------


if __name__ == "__main__":
    main()