import netCDF4 as nc
import argparse
import os
import json as  simplejson
import sys
import pylab as plt
from matplotlib import colors, ticker, cm

import pickle

import numpy as np
import itertools
from scipy import linalg
from sklearn import mixture

sys.path.append("..")
from io_read_in_files import read_in_netcdf

from thermodynamics import sat_adj
from thermodynamics import thetali
# from thermodynamics_1 import sat_adj


label_size = 10
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['xtick.direction']='out'
plt.rcParams['ytick.direction']='out'
'''
TO DO:
- only in levels that are in cloud layer (condition on ql!= 0)

(1) compute th_l from entropy (??? even in output?)
(2) compute bivariate Gaussian PDF for (th_l, q_t)
    (a) compute mean and covariance
(3) test quality of this PDF fitting --> how?
    (a) fit Kernel-Density-Estimated PDF & compare to Gaussian PDF (relative entropy minimisation)
    ???? is optimal number of kernels chosen ???
    !!!!! output / compute optimal number of kernels !!!!
    !!!!! (b) compare two PDFs (measure difference) !!!!!
'''


'''
Cloud Closer Schemes:

'''

def main():
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("path")
    parser.add_argument("casename")
    args = parser.parse_args()
    print('_______________________')
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
    # zrange = np.arange(0,36,2)
    zrange = np.arange(6, 20, 8)
    print('zrange', zrange*dz)
    print('_______________________')
    if case_name == 'DCBLSoares':
        var_list = ['u','w','s']
    else:
        var_list = ['w','s','qt']
        # var_list = ['w', 's']


    global ncomp
    global nvar
    ncomp = 1
    nvar = 2
    data_all = np.ndarray(shape=(0, nvar))
    nc_file_name_out = 'CC_' + str(d[0:-3]) + '_alltime.nc'
    create_statistics_file(os.path.join(fullpath_out, 'CloudClosure'), nc_file_name_out, ncomp, nvar, len(zrange))

    for i in range(len(zrange)):
        iz = zrange[i]
        for d in files:
            '''(1) compute liquid potential temperature from temperature and moisture'''
            p0 = 1e5
            # T, ql, qi = sat_adj(p0, 6500, 1e-3)

            nc_file_name = str(d)
            fullpath_in = os.path.join(in_path, 'fields', nc_file_name)
            print('fullpath_in', fullpath_in)
            T = read_in_netcdf('temperature', 'fields', fullpath_in)
            qt = read_in_netcdf('qt', 'fields', fullpath_in)
            ql = read_in_netcdf('ql', 'fields', fullpath_in)
            qi = np.zeros(shape=T.shape)
            theta_l = thetali(p0,T,qt,ql,qi)

            data = np.ndarray(shape=((nx * ny), nvar))
            means_ = np.ndarray(shape=(len(zrange), ncomp, nvar))
            covariance_ = np.zeros(shape=(len(zrange), ncomp, nvar, nvar))

            data1_ = theta_l.reshape((nx * ny), nz)
            data2_ = qt.reshape((nx * ny), nz)
            data[:, 0] = data1_[:, iz]
            data[:, 1] = data2_[:, iz]
            data_all = np.append(data_all, data, axis=0)

        '''(2) Compute bivariate Gaussian PDF (theta_l, qt) '''
        # means, covariance, weights = Gaussian_mixture_bivariate(data, var1, var2, np.int(d[0:-3]), iz*dz)
        clf = Gaussian_bivariate(data, 'T', 'qt', np.int(d[0:-3]), iz * dz)
        means_[i, :, :] = clf.means_[:, :]
        covariance_[i,:,:,:] = clf.covariances_[:,:,:]

        '''(3) Compute Kernel-Estimate PDF '''
        kde, kde_aux = Kernel_density_estimate(data, 'T', 'qt', np.int(d[0:-3]), iz * dz)

        relative_entropy(data, clf, kde)

        '''(4) Save Gaussian Mixture PDFs '''
    dump_variable(os.path.join(fullpath_out, 'CloudClosure', nc_file_name_out), 'means', means_, 'qtT', ncomp, nvar, len(zrange))
    dump_variable(os.path.join(fullpath_out, 'CloudClosure', nc_file_name_out), 'covariances', covariance_, 'qtT', ncomp, nvar, len(zrange))

    return


#----------------------------------------------------------------------
def relative_entropy(data, clf, kde):
    print('Relative Entropy')
    # rel_ent_clf = D(p_clf || p_kde)
    # rel_ent_kde = D(p_kde || p_clf )

    rel_ent_clf = 0
    rel_ent_kdf = 0
    n_sample = 50
    x_ = np.linspace(np.amin(data[:, 0]), np.amax(data[:, 0]), n_sample)
    y_ = np.linspace(np.amin(data[:, 1]), np.amax(data[:, 1]), n_sample)
    X, Y = np.meshgrid(x_, y_)
    XX = np.array([X.ravel(), Y.ravel()]).T
    Z_clf = np.exp(clf.score_samples(XX)).reshape(X.shape)
    Z_kde = np.exp(kde.score_samples(XX)).reshape(X.shape)
    for i in range(n_sample):
        for j in range(n_sample):
            rel_ent_clf += Z_clf[i,j] * np.log(Z_clf[i,j] / Z_kde[i,j])
            rel_ent_kdf += Z_kde[i, j] * np.log(Z_kde[i, j] / Z_clf[i, j])

    print('rel entr D(clf || kdf): ', rel_ent_clf)
    print('rel entr D(kdf || clf): ', rel_ent_kdf)
    # rel_ent_clf_np = np.sum(Z_clf * np.log(Z_clf / Z_kde))
    # rel_ent_kdf_np = np.sum(Z_kde * np.log(Z_kde / Z_clf))

    return
#----------------------------------------------------------------------
def Gaussian_bivariate(data, var_name1, var_name2, time, z):
    global ncomp
    clf = mixture.GaussianMixture(n_components=ncomp,covariance_type='full')
    # clf = sklearn.mixture.GaussianMixture(n_components=2, covariance_type='full')
    clf.fit(data)
    print('')

    if var_name1 == 'qt' or var_name2 == 'qt':
        plot_PDF_samples_qt(data, var_name1, var_name2, clf, time, z)
        print('!!!! qt: factor 100')
    else:
        plot_PDF_samples(data, var_name1, var_name2, clf, time, z)
    # plot_PDF_samples(data, var_name1, var_name2, clf, time, z)

    # return clf.means_, clf.covariances_, clf.weights_
    return clf
#----------------------------------------------------------------------
def Gaussian_univariate(data, var_name, time, iz):
    for i in range(1):
        print('Gaussian mixture: '+ var_name + ', height: '+np.str(iz*dz))
        clf = mixture.GaussianMixture(n_components=1,covariance_type='full')

        clf.fit(data[:,i].reshape(nx*ny,1))

        # print(var_name + ': means=' + np.str(clf.means_))
        # print(var_name + ': covar=' + np.str(clf.covariances_))
        # print(var_name + ': ', clf.means_.shape, clf.covariances_.shape)

        plot_PDF_samples(data, var_name, clf, time, iz)

    return clf.means_, clf.covariances_, clf.weights_
#----------------------------------------------------------------------
def Kernel_density_estimate(data, var_name1, var_name2, time, z):
    from sklearn.neighbors.kde import KernelDensity
    ''' Kerne Density Estimation:
    from sklearn.neighbors import KernelDensity

    Parameters:
    - bandwidth: The bandwidth here acts as a smoothing parameter, controlling the tradeoff between bias and variance
    in the result. A large bandwidth leads to a very smooth (i.e. high-bias) density distribution.
    A small bandwidth leads to an unsmooth (i.e. high-variance) density distribution.
    'metric': 'euclidean' (distance metric to use. Note that not all metrics are valid with all algorithms.)
    'atol': 0 (The desired absolute tolerance of the result.)
    'leaf_size': 40
    'kernel': 'gaussian'
    'rtol': 0 (The desired relative tolerance of the result. )
    'breadth_first': True
    'metric_params': None
    'algorithm': 'auto'
    '''
    amp = 100
    data_aux = np.ndarray(shape=((nx * ny), nvar))
    data_aux[:, 0] = data[:, 0]
    data_aux[:, 1] = data[:, 1] * amp

    # construct a kernel density estimate of the distribution
    print(" - computing KDE in spherical coordinates")
    # kde = KernelDensity(bandwidth=0.04, metric='haversine',
    #                     kernel='gaussian', algorithm='ball_tree')
    # kde.fit(Xtrain[ytrain == i])

    # Plotting
    n_sample = 100
    x_ = np.linspace(np.amin(data[:, 0]), np.amax(data[:, 0]), n_sample)
    y_ = np.linspace(np.amin(data[:, 1]), np.amax(data[:, 1]), n_sample)
    X, Y = np.meshgrid(x_, y_)
    XX = np.array([X.ravel(), Y.ravel()]).T

    x_aux = np.linspace(np.amin(data_aux[:, 0]), np.amax(data_aux[:, 0]), n_sample)
    y_aux = np.linspace(np.amin(data_aux[:, 1]), np.amax(data_aux[:, 1]), n_sample)
    X_aux, Y_aux = np.meshgrid(x_aux, y_aux)
    XX_aux = np.array([X_aux.ravel(), Y_aux.ravel()]).T



    fig = plt.figure(figsize=(12, 16))
    plt.subplot(3, 2, 1)
    bw = 5e-2
    kde = KernelDensity(kernel='gaussian', bandwidth=bw).fit(data_aux)
    # kde.score_samples(data)
    # Z = np.exp(kde.score_samples(XX_aux)).reshape(X.shape)
    Z_log = kde.score_samples(XX_aux).reshape(X.shape)
    plt.scatter(data_aux[:, 0], data_aux[:, 1], s=5, alpha=0.2)
    ax1 = plt.contour(X_aux, Y_aux, Z_log)
    plt.colorbar(ax1, shrink=0.8)
    labeling(var_name1, var_name2, amp)
    plt.title('bw = '+str(bw))

    plt.subplot(3, 2, 2)
    bw = 3e-2
    kde = KernelDensity(kernel='gaussian', bandwidth=bw).fit(data_aux)
    # kde.score_samples(data)
    # Z = np.exp(kde.score_samples(XX_aux)).reshape(X.shape)
    Z_log = kde.score_samples(XX_aux).reshape(X.shape)
    plt.scatter(data_aux[:, 0], data_aux[:, 1], s=5, alpha=0.2)
    ax1 = plt.contour(X_aux, Y_aux, Z_log)
    plt.colorbar(ax1, shrink=0.8)
    labeling(var_name1, var_name2, amp)
    plt.title('bw = ' + str(bw))

    plt.subplot(3, 2, 3)
    bw = 1e-2
    kde = KernelDensity(kernel='gaussian', bandwidth=bw).fit(data_aux)
    # kde.score_samples(data)
    # Z = np.exp(kde.score_samples(XX_aux)).reshape(X.shape)
    Z_log = kde.score_samples(XX_aux).reshape(X.shape)
    plt.scatter(data_aux[:, 0], data_aux[:, 1], s=5, alpha=0.2)
    ax1 = plt.contour(X_aux, Y_aux, Z_log)
    plt.colorbar(ax1, shrink=0.8)
    labeling(var_name1, var_name2, amp)
    plt.title('bw = ' + str(bw))

    plt.subplot(3, 2, 4)
    bw = 8e-3
    kde = KernelDensity(kernel='gaussian', bandwidth=bw).fit(data_aux)
    # kde.score_samples(data)
    # Z = np.exp(kde.score_samples(XX_aux)).reshape(X.shape)
    Z_log = kde.score_samples(XX_aux).reshape(X.shape)
    plt.scatter(data_aux[:, 0], data_aux[:, 1], s=5, alpha=0.2)
    ax1 = plt.contour(X_aux, Y_aux, Z_log)
    plt.colorbar(ax1, shrink=0.8)
    labeling(var_name1, var_name2, amp)
    plt.title('bw = ' + str(bw))

    plt.subplot(3, 2, 5)
    bw = 5e-3
    kde = KernelDensity(kernel='gaussian', bandwidth=bw).fit(data_aux)
    # kde.score_samples(data)
    # Z = np.exp(kde.score_samples(XX_aux)).reshape(X.shape)
    Z_log = kde.score_samples(XX_aux).reshape(X.shape)
    plt.scatter(data_aux[:, 0], data_aux[:, 1], s=5, alpha=0.2)
    ax1 = plt.contour(X_aux, Y_aux, Z_log)
    plt.colorbar(ax1, shrink=0.8)
    labeling(var_name1, var_name2, amp)
    plt.title('bw = ' + str(bw))

    plt.subplot(3, 2, 6)
    bw = 2e-3
    kde = KernelDensity(kernel='gaussian', bandwidth=bw).fit(data_aux)
    # kde.score_samples(data)
    # Z = np.exp(kde.score_samples(XX_aux)).reshape(X.shape)
    Z_log = kde.score_samples(XX_aux).reshape(X.shape)
    plt.scatter(data_aux[:, 0], data_aux[:, 1], s=5, alpha=0.2)
    ax1 = plt.contour(X_aux, Y_aux, Z_log)
    plt.colorbar(ax1, shrink=0.8)
    labeling(var_name1, var_name2, amp)
    plt.title('bw = ' + str(bw))

    fig.suptitle('Cloud Closure: Kernel Density Estimate (gaussian)', fontsize=20)
    plt.savefig(os.path.join(fullpath_out,'CloudClosure_figures','CC_' + var_name1 + '_' + var_name2 + '_z' + str(np.int(z)) + 'm_KDE_alltime.png'))
    plt.close()

    print('KDE shapes: ', kde.score_samples(XX).shape, X.shape)
    print(kde.get_params())

    return kde, kde

#----------------------------------------------------------------------
def labeling(var_name1, var_name2, amp):
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    if var_name1 == 'qt':
        plt.xlabel(var_name1 + ' ( * ' + np.str(amp) + ')')
        plt.ylabel(var_name2)
    else:
        plt.xlabel(var_name1)
        plt.ylabel(var_name2 + ' ( * ' + np.str(amp) + ')')

    return
#----------------------------------------------------------------------
def covariance_estimate_from_multicomp_pdf(clf):
    '''
    Input:
        clf: Expectation Maximization (EM) Gaussian Mixture Model (GMM) with weights, and parameters of all PDF components
    Output:
        parameters of total PDF
    '''

    ncomp, nvar = clf.means_.shape
    mean_tot = np.zeros(nvar)
    covar_tot = np.zeros((nvar,nvar))
    for i in range(ncomp):
        mean_tot[:] += clf.weights_[i]*clf.means_[i,:]
        covar_tot[:,:] += clf.weights_[i]*clf.covariances_[i,:,:]

    return mean_tot, covar_tot


#----------------------------------------------------------------------
# Exner Function
def exner_c(p0):
    p_tilde = 1.0e5
    kappa = 0.285956175299
    return np.pow((p0/p_tilde),kappa)
# Dry potential temperature
def theta_c(p, T):
    return T / exner_c(p)
# Liquid ice potential temperature consistent with Triopoli and Cotton (1981)
def thetali_c(p, T, qt, ql, qi, L):
    cpd = 1004.0
    return theta_c(p, T) * np.exp(-L*(ql/(1.0 - qt) + qi/(1.0 -qt))/(T*cpd))


#----------------------------------------------------------------------
# def plot_PDF_samples_qt(data, data_aux, var_name1, var_name2, clf, clf_aux, time, z):
def plot_PDF_samples_qt(data, var_name1, var_name2, clf, time, z):
    global ncomp
    amp = 1e2
    print('')
    print('plot PDF samples qt, factor='+np.str(amp))

    import matplotlib.mlab as mlab
    import matplotlib.cm as cm

    data_aux = np.ndarray(shape=((nx * ny), nvar))
    data_aux[:, 0] = data[:, 0]
    data_aux[:, 1] = data[:, 1] * amp
    clf_aux = mixture.GaussianMixture(n_components=ncomp, covariance_type='full')
    clf_aux.fit(data_aux)

    det_ = np.linalg.det(clf.covariances_[0, :, :])
    fact_ = 1. / np.sqrt((2 * np.pi) ** 2 * det_)
    det_aux_ = np.linalg.det(clf_aux.covariances_[0, :, :])
    fact_aux_ = 1. / np.sqrt((2 * np.pi) ** 2 * det_aux_)

    # Plotting
    n_sample = 100
    x1_max = np.amax(data[:, 0])
    x1_min = np.amin(data[:, 0])
    x2_max = np.amax(data[:, 1])
    x2_min = np.amin(data[:, 1])
    x = np.linspace(x1_min, x1_max, n_sample)
    y = np.linspace(x2_min, x2_max, n_sample)
    X, Y = np.meshgrid(x, y)
    XX = np.array([X.ravel(), Y.ravel()]).T
    Z = clf.score_samples(XX).reshape(X.shape)
    x_aux = np.linspace(np.amin(data_aux[:,0]), np.amax(data_aux[:,0]), n_sample)
    y_aux = np.linspace(np.amin(data_aux[:, 1]), np.amax(data_aux[:, 1]), n_sample)
    X_aux, Y_aux = np.meshgrid(x_aux, y_aux)
    XX_aux = np.array([X_aux.ravel(), Y_aux.ravel()]).T
    Z_aux = clf_aux.score_samples(XX_aux).reshape(X_aux.shape)

    # mx1 = clf.means_[0, 0]
    # my1 = clf.means_[0, 1]
    # sx1 = np.sqrt(clf.covariances_[0, 0, 0])
    # sy1 = np.sqrt(clf.covariances_[0, 1, 1])
    # sxy1 = clf.covariances_[0, 1, 0]
    # Z1 = mlab.bivariate_normal(X, Y, sigmax=sx1, sigmay=sy1, mux=mx1, muy=my1, sigmaxy=sxy1)
    # mx1_aux = clf_aux.means_[0, 0]
    # my1_aux = clf_aux.means_[0, 1]
    # sx1_aux = np.sqrt(clf_aux.covariances_[0, 0, 0])
    # sy1_aux = np.sqrt(clf_aux.covariances_[0, 1, 1])
    # sxy1_aux = clf_aux.covariances_[0, 1, 0]
    # Z1_aux = mlab.bivariate_normal(X_aux, Y_aux, sigmax=sx1_aux, sigmay=sy1_aux, mux=mx1_aux, muy=my1_aux, sigmaxy=sxy1_aux)

    fig = plt.figure(figsize=(12, 16))
    plt.subplot(3, 2, 1)
    plt.scatter(data[:, 0], data[:, 1], s=5, alpha=0.2)
    ax1 = plt.contour(X, Y, Z, colors='w', levels=np.linspace(10, 20, 11))
    plt.colorbar(ax1,shrink=0.8)
    plt.xlim([x1_min,x1_max])
    plt.ylim([x2_min, x2_max])
    plt.title(var_name1 + var_name2 + ' (data), t=' + str(time) + ', z=' + str(z))
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    plt.subplot(3, 2, 2)
    plt.scatter(data[:, 0], data[:, 1], s=5, alpha=0.2)
    ax1 = plt.contour(X_aux, Y_aux, Z_aux, colors='w', levels=np.linspace(10, 20, 11))
    plt.colorbar(ax1, shrink=0.8)
    plt.xlim([x1_min, x1_max])
    plt.ylim([x2_min, x2_max])
    plt.title(var_name1 + var_name2 + ' (data), t=' + str(time) + ', z=' + str(z))
    labeling(var_name1, var_name2, amp)

    plt.subplot(3, 2, 3)
    plt.scatter(data[:, 0], data[:, 1], s=5, alpha=0.2)
    levels = np.linspace(0, fact_, 10)
    ax1 = plt.contour(X, Y, np.exp(Z), levels=levels, linewidths=2)
    plt.plot([clf.means_[0, 0]], [clf.means_[0, 1]], 'o', markersize=8)
    plt.colorbar(ax1, shrink=0.8)
    plt.title('EM PDF')
    plt.xlim([x1_min, x1_max])
    plt.ylim([x2_min, x2_max])
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    plt.subplot(3, 2, 4)
    plt.scatter(data_aux[:, 0], data_aux[:, 1], s=5, alpha=0.2)
    # ax1 = plt.contour(X_aux, Y_aux, np.exp(Z_aux), levels=levels, linewidths=2)
    ax1 = plt.contour(X_aux, Y_aux, np.exp(Z_aux), linewidths=2)
    plt.plot([clf_aux.means_[0, 0]], [clf_aux.means_[0, 1]], 'o', markersize=8)
    plt.colorbar(ax1, shrink=0.8)
    plt.title('EM PDF')
    labeling(var_name1, var_name2, amp)

    plt.subplot(3, 2, 5)
    ax1 = plt.contourf(X, Y, np.exp(Z))
    # plt.scatter(X, Y, s=2, alpha=0.5)
    plt.colorbar(ax1, shrink=0.8)
    plt.title('EM PDF')
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    plt.subplot(3, 2, 6)
    ax1 = plt.contourf(X_aux, Y_aux, np.exp(Z_aux))
    plt.colorbar(ax1, shrink=0.8)
    plt.title('EM PDF')
    labeling(var_name1, var_name2, amp)

    fig.suptitle('Cloud Closure: Univariate Gaussian PDF fit', fontsize=20)
    plt.savefig(
        os.path.join(
            fullpath_out,'CloudClosure_figures/CC_' + var_name1 + '_' + var_name2 + '_z'
                         + str(np.int(z)) + 'm_bivariate_alltime.png')
    )

    plt.close()
    return
#----------------------------------------------------------------------
def plot_PDF_samples(data, var_name1, var_name2, clf, time, z):
    import matplotlib.mlab as mlab
    import matplotlib.cm as cm

    det_ = np.linalg.det(clf.covariances_[0, :, :])
    fact_ = 1. / np.sqrt((2 * np.pi) ** 2 * det_)

    # Plotting
    n_sample = 300
    x1_max = np.amax(data[:, 0])
    x1_min = np.amin(data[:, 0])
    x2_max = np.amax(data[:, 1])
    x2_min = np.amin(data[:, 1])
    x = np.linspace(x1_min, x1_max, n_sample)
    y = np.linspace(x2_min, x2_max, n_sample)
    X, Y = np.meshgrid(x, y)
    XX = np.array([X.ravel(), Y.ravel()]).T
    Z = clf.score_samples(XX).reshape(X.shape)
    mx1 = clf.means_[0, 0]
    my1 = clf.means_[0, 1]
    sx1 = np.sqrt(clf.covariances_[0, 0, 0])
    sy1 = np.sqrt(clf.covariances_[0, 1, 1])
    sxy1 = clf.covariances_[0, 1, 0]
    Z1 = mlab.bivariate_normal(X, Y, sigmax=sx1, sigmay=sy1, mux=mx1, muy=my1, sigmaxy=sxy1)

    plt.figure(figsize=(12, 12))
    levels_tot = np.linspace(0, fact_, 10)
    if fact_ <= 2:
        levels_cont = np.arange(0, fact_, 0.2)
        levels_contf = np.arange(0, fact_, 0.2)
    elif fact_ <= 10:
        levels_cont = np.arange(0,fact_,0.5)
        levels_contf = np.arange(0,fact_,0.5)
    elif fact_ <= 20:
        levels_cont = np.arange(0,fact_,2)
        levels_contf = np.arange(0,fact_,2)
    elif fact_ <= 50:
        levels_cont = np.arange(0,fact_,5)
        levels_contf = np.arange(0,fact_,5)
    else:
        levels_cont = np.arange(0,fact_,20)
        levels_contf = np.arange(0,fact_,20)
    levels_comp = np.linspace(0, fact_, 7)

    plt.subplot(3, 2, 1)
    plt.scatter(data[:, 0], data[:, 1], s=2, alpha=0.05)
    ax1 = plt.contour(X, Y, Z, levels=np.linspace(10, 20, 2))
    plt.colorbar(ax1,shrink=0.8)
    plt.xlim([x1_min,x1_max])
    plt.ylim([x2_min, x2_max])
    plt.title(var_name1 + var_name2 + ' (data), t=' + str(time) + ', z=' + str(z))
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    plt.colorbar(ax1,shrink=0.8)
    plt.xlim([x1_min,x1_max])
    plt.ylim([x2_min, x2_max])
    plt.title(var_name1 + var_name2 + ' (data), t=' + str(time) + ', z=' + str(z))
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)

    plt.subplot(3, 2, 3)
    ax1 = plt.hist2d(data[:, 0], data[:, 1], bins=30, normed=True)
    plt.colorbar(shrink=0.8)
    ax2 = plt.contour(X, Y, np.exp(Z), levels=levels_cont, linewidths=1, colors='w')
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    plt.title('data histogram')
    plt.subplot(3, 2, 4)
    ax1 = plt.contourf(X, Y, np.exp(Z),levels=levels_contf)
    plt.colorbar(ax1, shrink=0.8)
    plt.title('EM PDF')
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    plt.subplot(3, 2, 5)
    plt.scatter(data[:, 0], data[:, 1], s=2, alpha=0.05)
    ax1 = plt.contour(X, Y, Z1, levels=levels_comp, linewidths=1.5)
    plt.plot([clf.means_[0, 0]], [clf.means_[0, 1]], 'wo', markersize=6)
    plt.colorbar(ax1, shrink=0.8)
    plt.title('f1, f2')
    plt.xlim([x1_min, x1_max])
    plt.ylim([x2_min, x2_max])
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    plt.subplot(3, 2, 6)
    plt.scatter(data[:, 0], data[:, 1], s=2, alpha=0.05)
    ax1 = plt.contour(X, Y, np.exp(Z), levels=levels_cont, linewidths=1.5)
    # ax1 = plt.contour(X, Y, Z1+Z2, linewidths=1.5)
    plt.plot([clf.means_[0, 0]], [clf.means_[0, 1]], 'wo', markersize=6)
    plt.colorbar(ax1, shrink=0.8)
    plt.xlim([x1_min, x1_max])
    plt.ylim([x2_min, x2_max])
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    plt.title('f = f1 + f2')

    plt.savefig(fullpath_out+'CloudClosure_figures/CC_bivariate_' + var_name1 + '_' + var_name2 + '_z' + str(
        np.int(z)) + 'm_alltime.png')

    plt.close()
    return



#----------------------------------------------------------------------
def create_statistics_file(path,file_name, ncomp, nvar, nz_):
    # ncomp: number of Gaussian components in EM
    # nvar: number of variables of multi-variate Gaussian components
    global time, zrange
    print('create file:', path, file_name)
    rootgrp = nc.Dataset(os.path.join(path,file_name), 'w', format='NETCDF4')
    dimgrp = rootgrp.createGroup('dims')
    means_grp = rootgrp.createGroup('means')
    means_grp.createDimension('nz', nz_)
    means_grp.createDimension('ncomp', ncomp)
    means_grp.createDimension('nvar', nvar)
    cov_grp = rootgrp.createGroup('covariances')
    cov_grp.createDimension('nz', nz_)
    cov_grp.createDimension('ncomp', ncomp)
    cov_grp.createDimension('nvar', nvar)
    weights_grp = rootgrp.createGroup('weights')
    weights_grp.createDimension('nz', nz_)
    weights_grp.createDimension('EM2', 2)
    ts_grp = rootgrp.createGroup('time')
    ts_grp.createDimension('nt',len(time)-1)
    var = ts_grp.createVariable('t','f8',('nt'))
    for i in range(len(time)-1):
        var[i] = time[i+1]
    z_grp = rootgrp.createGroup('z-profile')
    z_grp.createDimension('nz', len(zrange))
    var = z_grp.createVariable('height', 'f8', ('nz'))
    for i in range(len(zrange)):
        var[i] = zrange[i]
    rootgrp.close()
    # print('create file end')
    return

def dump_variable(path, group_name, data_, var_name, ncomp, nvar, nz_):
    print('-------- dump variable --------', var_name, group_name, path)
    # print('dump variable', path, group_name, var_name, data_.shape, ncomp, nvar)
    if group_name == 'means':
        add_means(path, var_name, ncomp, nvar)
        data = np.empty((nz_,ncomp,nvar), dtype=np.double, order='c')
        for i in range(nz_):
            for j in range(ncomp):
                for k in range(nvar):
                    data[i,j,k] = data_[i,j,k]
        write_mean(path, group_name, data, var_name)

    elif group_name == 'covariances':
        add_covariance(path, var_name, ncomp, nvar)
        data = np.empty((nz_, ncomp, nvar, nvar), dtype=np.double, order='c')
        for i in range(nz_):
            for j in range(ncomp):
                for k1 in range(nvar):
                    for k2 in range(nvar):
                        data[i, j, k1, k2] = data_[i, j, k1, k2]
        write_covar(path, group_name, data, var_name)

    elif group_name == 'weights':
        add_weights(path, var_name, ncomp, nvar)
        data = np.empty((nz_, ncomp), dtype=np.double, order='c')
        for i in range(nz_):
            for j in range(ncomp):
                data[i, j] = data_[i, j]
        write_weights(path, group_name, data, var_name)

    # write_field(path, group_name, data, var_name)
    # print('--------')
    return

def add_means(path, var_name, ncomp, nvar):
    print('add means: ', var_name, path)
    # rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
    rootgrp = nc.Dataset(path, 'r+')
    group = rootgrp.groups['means']
    var = group.createVariable(var_name, 'f8', ('nz', 'ncomp', 'nvar'))
    rootgrp.close()
    return

def add_covariance(path, var_name, ncomp, nvar):
    # print('add covariance: ', var_name, path)
    # rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
    rootgrp = nc.Dataset(path, 'r+')
    group = rootgrp.groups['covariances']
    var = group.createVariable(var_name, 'f8', ('nz', 'ncomp', 'nvar', 'nvar'))
    rootgrp.close()
    return

# def add_weights(path, var_name, ncomp, nvar):
#     # print('add weights: ', var_name, path)
#     # rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
#     rootgrp = nc.Dataset(path, 'r+')
#     group = rootgrp.groups['weights']
#     var = group.createVariable(var_name, 'f8', ('nz', 'EM2'))
#     rootgrp.close()
#     return


def write_mean(path, group_name, data, var_name):
    print('write mean:', path, var_name, data.shape)
    rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
    fieldgrp = rootgrp.groups[group_name]
    var = fieldgrp.variables[var_name]
    var[:, :, :] = data[:,:,:]
    rootgrp.close()
    return

def write_covar(path, group_name, data, var_name):
    # print('write covar:', path, var_name, data.shape)
    rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
    fieldgrp = rootgrp.groups[group_name]
    var = fieldgrp.variables[var_name]
    var[:, :, :,:] = data[:, :, :,:]
    rootgrp.close()
    return

# def write_weights(path, group_name, data, var_name):
#     # print('write weights:', path, var_name, data.shape)
#     rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
#     fieldgrp = rootgrp.groups[group_name]
#     var = fieldgrp.variables[var_name]
#     print(var.shape, data.shape)
#     var[:, :] = data[:,:]
#     rootgrp.close()
#     return

def write_field(path, group_name, data, var_name):
    # print('')
    # print('write field:', path, var_name, data.shape, group_name)
    rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
    fieldgrp = rootgrp.groups[group_name]
    # print('fieldgrp', fieldgrp)
    var = fieldgrp.variables[var_name]
    # var = data
       # var[:] = np.array(data)
    # var[:, :, :] = data
    rootgrp.close()
    return

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
    profiles_group = f[group_name]
    #Get access to the variable dataset
    variable_dataset = profiles_group[variable_name]
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