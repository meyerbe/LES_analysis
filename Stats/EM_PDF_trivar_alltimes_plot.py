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


# #----------------------------------------------------------------------
# def plot_PDF_samples(data, means, covariances, var_name, time, z):
#     import matplotlib.mlab as mlab
#     import matplotlib.cm as cm
#
#     amp = 100
#     if var_name[0]== 'qt' or var_name[1]== 'qt' or var_name[2] == 'qt':
#         pass
#     #     # print('plot PDF samples qt: amp = ' + np.str(amp))
#     #     data[:, 0] = data[:, 0] * amp
#     #     data[:, 1] = data[:, 1]
#     # elif var_name2 == 'qt':
#     #     # print('plot PDF samples qt: amp = ' + np.str(amp))
#     #     data[:, 0] = data[:, 0]
#     #     data[:, 1] = data[:, 1] * amp
#     #
#     det1 = np.linalg.det(covariances[0, :, :])
#     det2 = np.linalg.det(covariances[1, :, :])
#     det_ = min(det1, det2)
#     fact_ = 1. / np.sqrt((2 * np.pi) ** 2 * det_)
#     # fact = 1.1 * clf.weights_[0] * 1. / np.sqrt((2 * np.pi) ** 2 * det1) + clf.weights_[1] * 1. / np.sqrt(
#     #     (2 * np.pi) ** 2 * det2)
#
#     # # Plotting
#     n_sample = 300
#     x1_max = np.amax(data[:, 0])
#     x1_min = np.amin(data[:, 0])
#     x2_max = np.amax(data[:, 1])
#     x2_min = np.amin(data[:, 1])
#     x3_max = np.amax(data[:, 2])
#     x3_min = np.amin(data[:, 2])
#     x = np.linspace(x1_min, x1_max, n_sample)
#     y = np.linspace(x2_min, x2_max, n_sample)
#     X, Y = np.meshgrid(x, y)
#     XX = np.array([X.ravel(), Y.ravel()]).T
#     Z = clf.score_samples(XX).reshape(X.shape)
#     # mx1 = clf.means_[0, 0]
#     # my1 = clf.means_[0, 1]
#     # sx1 = np.sqrt(clf.covariances_[0, 0, 0])
#     # sy1 = np.sqrt(clf.covariances_[0, 1, 1])
#     # sxy1 = clf.covariances_[0, 1, 0]
#     # mx2 = clf.means_[1, 0]
#     # my2 = clf.means_[1, 1]
#     # sx2 = np.sqrt(clf.covariances_[1, 0, 0])
#     # sy2 = np.sqrt(clf.covariances_[1, 1, 1])
#     # sxy2 = clf.covariances_[1, 1, 0]
#     # Z1 = mlab.bivariate_normal(X, Y, sigmax=sx1, sigmay=sy1, mux=mx1, muy=my1, sigmaxy=sxy1)
#     # Z2 = mlab.bivariate_normal(X, Y, sigmax=sx2, sigmay=sy2, mux=mx2, muy=my2, sigmaxy=sxy2)
#     #
#     fig = plt.figure(figsize=(10, 12))
#     # levels_tot = np.linspace(0, fact, 10)
#     # if fact <= 2:
#     #     levels_cont = np.arange(0, fact, 0.2)
#     #     levels_contf = np.arange(0, fact, 0.2)
#     # elif fact <= 10:
#     #     levels_cont = np.arange(0,fact,0.5)
#     #     levels_contf = np.arange(0,fact,0.5)
#     # elif fact <= 20:
#     #     levels_cont = np.arange(0,fact,2)
#     #     levels_contf = np.arange(0,fact,2)
#     # elif fact <= 50:
#     #     levels_cont = np.arange(0,fact,5)
#     #     levels_contf = np.arange(0,fact,5)
#     # else:
#     #     levels_cont = np.arange(0,fact,20)
#     #     levels_contf = np.arange(0,fact,20)
#     # levels_comp = np.linspace(0, fact_, 7)
#     #
#     # plt.subplot(3, 2, 1)
#     # plt.scatter(data[:, 0], data[:, 1], s=2, alpha=0.3)
#     # ax1 = plt.contour(X, Y, Z, levels=np.linspace(10, 20, 2))
#     # plt.colorbar(ax1,shrink=0.8)
#     # plt.xlim([x1_min,x1_max])
#     # plt.ylim([x2_min, x2_max])
#     # plt.title(var_name1 + var_name2 + ' (data), t=' + str(time) + ', z=' + str(z))
#     # plt.xlim([x1_min,x1_max])
#     # plt.ylim([x2_min, x2_max])
#     # # plt.title(var_name1 + var_name2 + ' (data), t=' + str(time) + ', z=' + str(z))
#     # plt.title(var_name1 + var_name2 + ' (data)')
#     # axis_label(var_name1, var_name2, amp)
#     #
#     # plt.subplot(3, 2, 3)
#     # ax1 = plt.hist2d(data[:, 0], data[:, 1], bins=30, normed=True)
#     # plt.colorbar(shrink=0.8)
#     # ax2 = plt.contour(X, Y, np.exp(Z), levels=levels_cont, linewidths=1, colors='w')
#     # axis_label(var_name1,var_name2, amp)
#     # plt.title('data histogram')
#     # plt.subplot(3, 2, 4)
#     # ax1 = plt.contourf(X, Y, np.exp(Z),levels=levels_contf)
#     # plt.colorbar(ax1, shrink=0.8)
#     # plt.title('EM PDF')
#     # axis_label(var_name1, var_name2, amp)
#     # plt.subplot(3, 2, 5)
#     # plt.scatter(data[:, 0], data[:, 1], s=2, alpha=0.05)
#     # ax1 = plt.contour(X, Y, Z1, levels=levels_comp, linewidths=1.5)
#     # ax2 = plt.contour(X, Y, Z2, levels=levels_comp, linewidths=1.5)
#     # plt.plot([clf.means_[0, 0]], [clf.means_[0, 1]], 'wo', markersize=6)
#     # plt.plot([clf.means_[1, 0]], [clf.means_[1, 1]], 'wo', markersize=6)
#     # plt.colorbar(ax1, shrink=0.8)
#     # plt.title('f1, f2')
#     # plt.xlim([x1_min, x1_max])
#     # plt.ylim([x2_min, x2_max])
#     # axis_label(var_name1, var_name2, amp)
#     # plt.subplot(3, 2, 6)
#     # plt.scatter(data[:, 0], data[:, 1], s=2, alpha=0.05)
#     # ax1 = plt.contour(X, Y, np.exp(Z), levels=levels_cont, linewidths=1.5)
#     # plt.plot([clf.means_[0, 0]], [clf.means_[0, 1]], 'wo', markersize=6)
#     # plt.plot([clf.means_[1, 0]], [clf.means_[1, 1]], 'wo', markersize=6)
#     # plt.colorbar(ax1, shrink=0.8)
#     # plt.xlim([x1_min, x1_max])
#     # plt.ylim([x2_min, x2_max])
#     # axis_label(var_name1, var_name2, amp)
#     # plt.title('f = f1 + f2')
#
#     # fig.suptitle('EM2 PDF: ' + var_name1 + var_name2 + ' (t=' + str(time) + ', z=' + str(z) +'m)', fontsize=20)
#     # plt.savefig(os.path.join(
#     #     fullpath_out,'EM2_trivar_figures','EM2_PDF_trivariate_' + var_name1 + '_' + var_name2 + '_' + str(time) + '_z' + str(
#     #         np.int(z)) + 'm.png')
#     # )
#     #
#     plt.close()
#     return
# #----------------------------------------------------------------------
# def covariance_estimate_from_multicomp_pdf(means, covars, weights):
#     '''
#     Input: clf - Expectation Maximization (EM) Gaussian Mixture Model (GMM) with weights, and parameters of all PDF components
#     Output: parameters of total PDF
#     '''
#     z_max, ncomp, nvar = means.shape
#     mean_tot = np.zeros((z_max, nvar))
#     covar_tot = np.zeros((z_max, nvar,nvar))
#     for i in range(ncomp):
#         for z in range(z_max):
#             mean_tot[z,:] += weights[z,i]*means[z,i,:]
#             covar_tot[z,:] += weights[z,i]*covars[z,i,:,:]
#
#     return mean_tot, covar_tot
#
#
# #----------------------------------------------------------------------
# def bivar_plot_means(var_name, means_, covars_, mean_tot_, covar_tot_, z0,t0, time_):
#     print('plotting means')
#     print(os.path.join(fullpath_out, 'figures_EM2_bivar'))
#     global time, ncomp, zrange_
#     nt = time.size
#     colors = ['b', 'g']
#
#     # for k in range(z_max):
#     #     for n in range(nt):
#     #         for j in range(2):
#     #             if means_[n, k, 0, j] < means_[n, k, 1, j]:
#     #                 aux = means_[n, k, 1, j]
#     #                 means_[n, k, 1, j] = means_[n, k, 0, j]
#     #                 means_[n, k, 0, j] = aux
#
#     # over time at given level
#     means = means_[:,z0,:,:]
#     covars = covars_[:,z0,:,:]
#     means_tot = mean_tot_[:,z0,:]
#     covar_tot = covar_tot_[:, z0, :, :]
#
#     fig = plt.figure(figsize=(10,8))
#     fig.suptitle(var_name+ r': mean values of $f_1$, $f_2$ (z=' + np.str(zrange_[z0] * dz) + 'm)', fontsize=20)
#     for j in range(2):      # loop over variables
#         plt.subplot(2,2,j+1)
#         plt.plot(time[:], means[:, 0, j], 'o-')
#         plt.plot(time[:], means[:, 1, j], 'o-')
#         plt.plot(time[:], means_tot[:, j], 'ro-')
#         # plt.title('means '+var_name[j]+' (z=' + np.str(zrange_[z0] * dz) + 'm)')
#         plt.xlabel('time')
#         plt.ylabel(r'$<$'+var_name[j]+r'$>$')
#
#         plt.subplot(2, 2, 2+j+1)
#         plt.plot(time[:],means[:,0,j],'o',color=colors[0])
#         plt.plot(time[:], means[:, 1, j], 'o', color=colors[1])
#         plt.plot(time[:], means_tot[:, j], 'ro')
#         for comp in range(ncomp):
#             for n in range(nt):
#                 bar = 0.5*np.sqrt(covars[n,comp,j,j])
#                 plt.plot([time[n],time[n]],[means[n,comp,j]-bar, means[n,comp,j]+bar],color=colors[comp])
#                 bar = 0.5 * np.sqrt(covar_tot[n, j, j])
#                 plt.plot([time[n], time[n]], [means_tot[n, j] - bar, means_tot[n, j] + bar], color='r')
#         # plt.title('means ' + var_name[j] + ' (z=' + np.str(zrange_[z0] * dz) + 'm)')
#         plt.xlabel('time')
#         plt.ylabel(r'$<$' + var_name[j] + r'$>$')
#     plt.savefig(os.path.join(fullpath_out,'EM2_bivar_figures/') + 'means_time_a_' + var_name + '_z' + str(np.int(zrange_[z0]*dz)) + 'm.png')
#     plt.close()
#
#     # over levels, at given time
#     means = means_[t0, :, :, :]
#     covars = covars_[t0, :, :, :]
#     means_tot = mean_tot_[t0, :, :]
#     covar_tot = covar_tot_[t0, :, :, :]
#
#     fig = plt.figure(figsize=(10, 8))
#     fig.suptitle(var_name + r': mean values of $f_1$, $f_2$ (t=' + np.str(time_[t0]) + ')', fontsize=20)
#     for j in range(2):  # loop over variables
#         plt.subplot(2, 2, j + 1)
#         for comp in range(ncomp):
#             plt.plot(means[:, comp, j], zrange_[:] * dz, 'o-', color=colors[comp], label='comp ' + np.str(comp))
#         plt.plot(means_tot[:, j], zrange_[:] * dz, 'ro-', label='total')
#         plt.xlabel(r'$<$' + var_name[j] + r'$>$')
#         plt.ylabel('height z')
#
#         plt.subplot(2, 2, 2 + j + 1)
#         for comp in range(ncomp):
#             plt.plot(means[:, comp, j], zrange_[:] * dz, 'o', markersize=4,
#                      color=colors[comp], label='comp '+np.str(comp))
#             bar = 0.5 * np.sqrt(covars[:, comp, j, j])
#             plt.plot([means[:, comp, j] - bar, means[:, comp, j] + bar], [zrange_[:] * dz, zrange_[:] * dz], color=colors[comp])
#         plt.xlabel(r'$<$' + var_name[j] + r'$>$')
#         plt.ylabel('height z')
#     ax = plt.subplot(2, 2, 1)
#     ax.legend()
#     plt.savefig(
#             os.path.join(fullpath_out, 'EM2_bivar_figures/') + 'means_levels_' + var_name + '_t' + str(np.int(time_[t0])) + '.png')
#     plt.close()
#
#     # plt.figure()
#     # for comp in range(ncomp):
#     #     for i in range(z_max):
#     #         plt.plot(i * dz, means[i, comp, 0], 'o',color=colors[comp])
#     #         bar = 0.5 * np.sqrt(covars[i, comp, 0, 0])
#     #         plt.plot([i * dz,i * dz],[means[i,comp,0]-bar,means[i,comp,0]+bar],color=colors[comp])
#     # plt.title('means (t=' + np.str(time_[t0]) + ')')
#     # plt.xlabel('height z')
#     # plt.ylabel(var_name)
#     # plt.savefig(
#     #     os.path.join(fullpath_out, 'EM2_bivar_figures/') + 'means_levels_a_' + var_name + '_t' + str(np.int(time_[t0])) + '.png')
#     # plt.figure()
#     # for comp in range(ncomp):
#     #     for i in range(z_max):
#     #         plt.plot(i * dz, means[i, comp, 0], 'o', color=colors[comp])
#     # plt.title('means (t=' + np.str(time_[t0]) + ')')
#     # plt.xlabel('height z')
#     # plt.ylabel(var_name)
#     # plt.savefig(
#     #     os.path.join(fullpath_out, 'EM2_bivar_figures/') + 'means_levels_b_' + var_name + '_t' + str(np.int(time_[t0])) + '.png')
#     return
#
# #----------------------------------------------------------------------
# def bivar_plot_covars(var_name, means_, covars_, mean_tot_, covars_tot_, z0,t0,time_):
#     print('plotting covars', covars_.shape)
#     print(os.path.join(fullpath_out, 'figures_EM2_bivar'))
#     global time, ncomp, zrange_
#     nt = time.size
#     colors = ['b', 'g']
#     # for i1 in range(2):  # loop over variables
#     #     for n in range(nt):
#     #         for k in range(z_max):
#     #             if means_[n,k, 0, i1] < means_[n,k, 1, i1]:
#     #                 aux = means_[n,k,1,i1]
#     #                 means_[n,k,1,i1] = means_[n,k,0,i1]
#     #                 means_[n,k,0,i1] = aux
#     #                 for i2 in range(2):
#     #                     aux = covars_[n,k,0,i1,i2]
#     #                     covars_[n,k,0,i1,i2] = covars_[n,k,1,i1,i2]
#     #                     covars_[n,k,1,i1,i2] = aux
#
#     # over time at given level
#     means = means_[:,z0,:,:]
#     covars = covars_[:,z0,:,:,:]
#     covars_tot = covars_tot_[:,z0,:,:]
#     print(time.shape, means.shape)
#     fig = plt.figure(figsize=(12,5))
#     fig.suptitle(var_name + r': covariance values of $f_1$, $f_2$ (z=' + np.str(zrange_[z0] * dz) + 'm)', fontsize=20)
#     n = 1
#     for i in range(2):
#         for j in range(i,2):
#             plt.subplot(1, 3, n)
#             for comp in range(ncomp):
#                 plt.plot(time[:],covars[:,comp,i,j],'o-', label='comp '+np.str(comp))
#                 plt.xlabel('time')
#                 plt.ylabel(var_name[i]+var_name[j])
#             plt.plot(time[:], covars_tot[:,i,j],'ro-',label='total')
#             n += 1
#     ax = plt.subplot(1,3,1)
#     ax.legend()
#     plt.savefig(os.path.join(fullpath_out,'EM2_bivar_figures/') + 'covars_time_' + var_name + '_z' + str(np.int(zrange_[z0]*dz)) + 'm.png')
#     plt.close()
#
#     # over levels, at given time
#     means = means_[t0, :, :, :]
#     covars = covars_[t0, :,:, :, :]
#     covars_tot = covars_tot_[t0, :, :, :]
#     fig = plt.figure(figsize=(12, 5))
#     fig.suptitle(var_name + r': covariance values of $f_1$, $f_2$ (t=' + np.str(time[t0]) + ')', fontsize=20)
#     n = 1
#     for i in range(2):
#         for j in range(i, 2):
#             plt.subplot(1, 3, n)
#             for comp in range(ncomp):
#                 plt.plot(covars[:, comp, i, j], dz*zrange_, 'o-', color=colors[comp], label='comp '+np.str(comp))
#                 plt.xlabel(var_name[i] + var_name[j])
#                 plt.ylabel('height z')
#             plt.plot(covars_tot[:, i, j], dz*zrange_[:], 'ro-', label='total')
#             n += 1
#     ax = plt.subplot(1, 3, 1)
#     ax.legend()
#     plt.savefig(
#         os.path.join(fullpath_out, 'EM2_bivar_figures/') + 'covars_level_' + var_name + '_t' + str(np.int(time_[t0])) + '.png')
#     plt.close()
#     return


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