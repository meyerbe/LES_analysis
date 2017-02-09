import netCDF4 as nc
import argparse
import os
import sys
import json as  simplejson
import pylab as plt
from matplotlib.colors import LogNorm
from matplotlib import colors, ticker, cm

import pickle

import numpy as np
import itertools
from scipy import linalg
from sklearn import mixture

sys.path.append("..")
from io_read_in_files import read_in_netcdf_fields

label_size = 10
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['xtick.direction']='out'
plt.rcParams['ytick.direction']='out'


'''
Bivariate Gaussian Mixture Model = Superposition of multiple Gaussian Distributions
(1) Fit A GMM to data points to get means, covariance matrices and relative weights
    --> use the expectation-maximization algorithm
(2) check out to which mode each data point belongs (class probability)

http://scikit-learn.org/stable/modules/generated/sklearn.mixture.GaussianMixture.html#sklearn.mixture.GaussianMixture
The GaussianMixture object implements the expectation-maximization (EM) algorithm for fitting mixture-of-Gaussian models.
- GaussianMixture.fit method: learns a Gaussian Mixture Model from train data.
- GaussianMixture.predict: Given test data, it can assign to each sample the Gaussian it mostly probably belong to.
- different options to constrain the covariance of the difference classes estimated: spherical, diagonal, tied or full covariance.

 First one assumes random components (randomly centered on data points, learned from k-means, or even just normally distributed
  around the origin) and computes for each point a probability of being generated by each component of the model.
  Then, one tweaks the parameters to maximize the likelihood of the data given those assignments.
  Repeating this process is guaranteed to always converge to a local optimum.

BIC: Bayesian Information Criterion

sklearn.mixture.GaussianMixture(n_compontens=1,covariance_type='full',
tol=0.001,reg_covar=1e-06,max_iter=100,n_init=1,
init_params='kmeans',weights_init=None,mans_init=None,precisions_init=None)


Methods:
- aic
- bic
- fit(data):  train model with data
- get_params
- predict(X[, y]): Predict the labels for the data samples in X using trained model.
- predict_proba
- sample
- score
- score_samples(X): Generate random samples from the fitted Gaussian distribution.
- set_params

--------------------
Input:
3D fields at given time-step

Output:
means(z) = [m1(z),m2(z)] --> shape = nz x ncomp x nvar
covars(z) = [[c11(z),c12(z)],[c21(z),c22(z)]] --> shape = nz x ncomp x nvar x nvar

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
    # ______________________
    '''
    zrange:     z-values for which the PDF is fitted
    var_list:   list of variables that are included in (multi-variate) PDF
    '''
    global zrange
    zrange = np.arange(0,72,2)
    print('zrange', zrange)
    print('_______________________')
    if case_name == 'DCBLSoares':
        var_list = ['u','w','s']
    else:
        var_list = ['w','s','qt']
        # var_list = ['w', 's']



    '''
    Bi-variate PDF for (s,qt,w)
    '''
    global ncomp
    global nvar
    ncomp = 2
    nvar = 2
    data = np.ndarray(shape=((nx * ny), nvar))
    means_ = np.ndarray(shape=(len(zrange), ncomp, nvar))
    covariance_ = np.zeros(shape=(len(zrange), ncomp, nvar, nvar))
    weights_ = np.zeros(shape=(len(zrange), 2))
    mean_tot = np.ndarray(shape=(len(zrange), nvar))
    covariance_tot = np.zeros(shape=(len(zrange), nvar, nvar))

    count_t = 0
    for d in files:
        nc_file_name = 'EM2_bivar_' + str(d)
        create_statistics_file(os.path.join(fullpath_out,'EM2_bivar'), nc_file_name, ncomp, nvar, len(zrange))

        fullpath_in = os.path.join(args.path, 'fields', d)
        print(fullpath_in)
        for n1 in range(len(var_list)):
            var1 = var_list[n1]
            data1_ = read_in_netcdf_fields(var1,fullpath_in).reshape((nx*ny,nz))
            for n2 in range(n1,len(var_list)):
                var2 = var_list[n2]
                print(var1, var2)
                if var1 == var2:
                    continue
                data2_ = read_in_netcdf_fields(var2, fullpath_in).reshape((nx*ny,nz))
                for i in range(len(zrange)):
                    iz = zrange[i]
                    data[:, 0] = data1_[:, iz]
                    data[:, 1] = data2_[:, iz]

                    # means, covariance, weights = Gaussian_mixture_bivariate(data, var1, var2, np.int(d[0:-3]), iz*dz)
                    clf = Gaussian_mixture_bivariate(data, var1, var2, np.int(d[0:-3]), iz * dz)
                    means_[i, :, :] = clf.means_[:, :]
                    covariance_[i,:,:,:] = clf.covariances_[:,:,:]
                    weights_[i,:] = clf.weights_[:]

                    mean_tot[i,:], covariance_tot[i,:,:] = covariance_estimate_from_multicomp_pdf(clf)

                dump_variable(os.path.join(fullpath_out, 'EM2_bivar', nc_file_name), 'means', means_, var1+var2, ncomp, nvar, len(zrange))
                dump_variable(os.path.join(fullpath_out, 'EM2_bivar', nc_file_name), 'covariances', covariance_, var1+var2, ncomp, nvar, len(zrange))
                dump_variable(os.path.join(fullpath_out, 'EM2_bivar', nc_file_name), 'weights', weights_, var1+var2, ncomp, nvar, len(zrange))

                # dump_variable(os.path.join(fullpath_out, nc_file_name), 'means', mean_tot, var1 + var2+'tot', ncomp, nvar,len(zrange))
                # dump_variable(os.path.join(fullpath_out, nc_file_name), 'covariances', covariance_tot, var1 + var2+'tot', ncomp, nvar,len(zrange))

        if var1 == 'w' and var2 == 's':
            means_time_ws[count_t,:,:,:] = means_[:,:,:]
            covariance_time_ws[count_t, :, :, :] = covariance_[:, :, :]
        count_t += 1
    # z0 = 2
    # plot_means(means_time_ws, 'w', 's', z0)
    return


#----------------------------------------------------------------------
def Gaussian_mixture_bivariate(data, var_name1, var_name2, time, z):
    print('Gaussian Mixture bivariate: ' + var_name1 + var_name2)
    clf = mixture.GaussianMixture(n_components=2,covariance_type='full')
    # clf = sklearn.mixture.GaussianMixture(n_components=2, covariance_type='full')
    clf.fit(data)
    print('')
    # mean_tot, covariance_tot = covariance_estimate_from_multicomp_pdf(clf)

    amp_qt = 1e2
    amp_w = 1e0
    data_aux = np.ndarray(shape=((nx * ny), nvar))

    if var_name1 == 'qt':
        data_aux[:, 0] = data[:, 0] * amp_qt
    else:
        data_aux[:, 0] = data[:, 0]
    if var_name2 == 'qt':
        data_aux[:, 1] = data[:, 1] * amp_qt
    else:
        data_aux[:, 1] = data[:, 1]

    clf_aux = mixture.GaussianMixture(n_components=2, covariance_type='full')
    clf_aux.fit(data_aux)

    if var_name1 == 'qt' or var_name2 == 'qt':
        plot_PDF_samples(data_aux, var_name1, var_name2, clf_aux, amp_qt, time, z)
        # plot_PDF_samples_log(data_aux, var_name1, var_name2, clf_aux, amp_qt, time, z)
        # return clf_aux.means_, clf_aux.covariances_
        return clf_aux
    else:
        plot_PDF_samples(data, var_name1, var_name2, clf, amp_qt, time, z)
        plot_PDF_samples_log(data, var_name1, var_name2, clf, amp_qt, time, z)
        return clf

    # plot_PDF_samples(data, var_name1, var_name2, clf, time, z)
    # return clf.means_, clf.covariances_, clf.weights_
    # return clf

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
def plot_PDF_samples(data, var_name1, var_name2, clf, amp, time, z):
    import matplotlib.mlab as mlab
    import matplotlib.cm as cm

    # if var_name1 == 'qt':
    #     # print('plot PDF samples qt: amp = ' + np.str(amp))
    #     data[:, 0] = data[:, 0] * amp
    #     data[:, 1] = data[:, 1]
    # elif var_name2 == 'qt':
    #     # print('plot PDF samples qt: amp = ' + np.str(amp))
    #     data[:, 0] = data[:, 0]
    #     data[:, 1] = data[:, 1] * amp

    det1 = np.linalg.det(clf.covariances_[0, :, :])
    det2 = np.linalg.det(clf.covariances_[1, :, :])
    det_ = min(det1, det2)
    fact_ = 1. / np.sqrt((2 * np.pi) ** 2 * det_)
    fact = 1.1 * clf.weights_[0] * 1. / np.sqrt((2 * np.pi) ** 2 * det1) + clf.weights_[1] * 1. / np.sqrt(
        (2 * np.pi) ** 2 * det2)

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
    mx2 = clf.means_[1, 0]
    my2 = clf.means_[1, 1]
    sx2 = np.sqrt(clf.covariances_[1, 0, 0])
    sy2 = np.sqrt(clf.covariances_[1, 1, 1])
    sxy2 = clf.covariances_[1, 1, 0]
    Z1 = mlab.bivariate_normal(X, Y, sigmax=sx1, sigmay=sy1, mux=mx1, muy=my1, sigmaxy=sxy1)
    Z2 = mlab.bivariate_normal(X, Y, sigmax=sx2, sigmay=sy2, mux=mx2, muy=my2, sigmaxy=sxy2)

    fig = plt.figure(figsize=(10, 12))
    levels_tot = np.linspace(0, fact, 10)
    if fact <= 2:
        levels_cont = np.arange(0, fact, 0.2)
        levels_contf = np.arange(0, fact, 0.2)
    elif fact <= 10:
        levels_cont = np.arange(0,fact,0.5)
        levels_contf = np.arange(0,fact,0.5)
    elif fact <= 20:
        levels_cont = np.arange(0,fact,2)
        levels_contf = np.arange(0,fact,2)
    elif fact <= 50:
        levels_cont = np.arange(0,fact,5)
        levels_contf = np.arange(0,fact,5)
    else:
        levels_cont = np.arange(0,fact,20)
        levels_contf = np.arange(0,fact,20)
    levels_comp = np.linspace(0, fact_, 7)

    plt.subplot(3, 2, 1)
    plt.scatter(data[:, 0], data[:, 1], s=2, alpha=0.3)
    ax1 = plt.contour(X, Y, Z, levels=np.linspace(10, 20, 2))
    plt.colorbar(ax1,shrink=0.8)
    plt.xlim([x1_min,x1_max])
    plt.ylim([x2_min, x2_max])
    plt.title(var_name1 + var_name2 + ' (data), t=' + str(time) + ', z=' + str(z))
    plt.xlim([x1_min,x1_max])
    plt.ylim([x2_min, x2_max])
    # plt.title(var_name1 + var_name2 + ' (data), t=' + str(time) + ', z=' + str(z))
    plt.title(var_name1 + var_name2 + ' (data)')
    axis_label(var_name1, var_name2, amp)

    plt.subplot(3, 2, 3)
    ax1 = plt.hist2d(data[:, 0], data[:, 1], bins=30, normed=True)
    plt.colorbar(shrink=0.8)
    ax2 = plt.contour(X, Y, np.exp(Z), levels=levels_cont, linewidths=1, colors='w')
    axis_label(var_name1,var_name2, amp)
    plt.title('data histogram')
    plt.subplot(3, 2, 4)
    ax1 = plt.contourf(X, Y, np.exp(Z), levels=levels_contf)
    plt.colorbar(ax1, shrink=0.8)
    plt.title('EM PDF')
    axis_label(var_name1, var_name2, amp)
    plt.subplot(3, 2, 5)
    plt.scatter(data[:, 0], data[:, 1], s=2, alpha=0.05)
    ax1 = plt.contour(X, Y, Z1, levels=levels_comp, linewidths=1.5)
    ax2 = plt.contour(X, Y, Z2, levels=levels_comp, linewidths=1.5)
    plt.plot([clf.means_[0, 0]], [clf.means_[0, 1]], 'wo', markersize=6)
    plt.plot([clf.means_[1, 0]], [clf.means_[1, 1]], 'wo', markersize=6)
    plt.colorbar(ax1, shrink=0.8)
    plt.title('f1, f2')
    plt.xlim([x1_min, x1_max])
    plt.ylim([x2_min, x2_max])
    axis_label(var_name1, var_name2, amp)
    plt.subplot(3, 2, 6)
    plt.scatter(data[:, 0], data[:, 1], s=2, alpha=0.05)
    ax1 = plt.contour(X, Y, np.exp(Z), levels=levels_cont, linewidths=1.5)
    plt.plot([clf.means_[0, 0]], [clf.means_[0, 1]], 'wo', markersize=6)
    plt.plot([clf.means_[1, 0]], [clf.means_[1, 1]], 'wo', markersize=6)
    plt.colorbar(ax1, shrink=0.8)
    plt.xlim([x1_min, x1_max])
    plt.ylim([x2_min, x2_max])
    axis_label(var_name1, var_name2, amp)
    plt.title('f = f1 + f2')

    fig.suptitle('EM2 PDF: ' + var_name1 + var_name2 + ' (t=' + str(time) + ', z=' + str(z) +'m)', fontsize=20)
    plt.savefig(os.path.join(
        fullpath_out,'EM2_bivar_figures','EM2_PDF_bivariate_' + var_name1 + '_' + var_name2 + '_' + str(time) + '_z' + str(
            np.int(z)) + 'm.png')
    )

    plt.close()
    return
#----------------------------------------------------------------------
def plot_PDF_samples_log(data, var_name1, var_name2, clf, amp, time, z):
    import matplotlib.mlab as mlab
    import matplotlib.cm as cm

    det1 = np.linalg.det(clf.covariances_[0, :, :])
    det2 = np.linalg.det(clf.covariances_[1, :, :])
    det_ = min(det1, det2)
    fact_ = 1. / np.sqrt((2 * np.pi) ** 2 * det_)
    fact = 1.1 * clf.weights_[0] * 1. / np.sqrt((2 * np.pi) ** 2 * det1) + clf.weights_[1] * 1. / np.sqrt(
        (2 * np.pi) ** 2 * det2)

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
    mx2 = clf.means_[1, 0]
    my2 = clf.means_[1, 1]
    sx2 = np.sqrt(clf.covariances_[1, 0, 0])
    sy2 = np.sqrt(clf.covariances_[1, 1, 1])
    sxy2 = clf.covariances_[1, 1, 0]
    Z1 = mlab.bivariate_normal(X, Y, sigmax=sx1, sigmay=sy1, mux=mx1, muy=my1, sigmaxy=sxy1)
    Z2 = mlab.bivariate_normal(X, Y, sigmax=sx2, sigmay=sy2, mux=mx2, muy=my2, sigmaxy=sxy2)

    fig = plt.figure(figsize=(10, 12))
    levels_tot = np.linspace(0, fact, 10)
    if fact <= 2:
        levels_cont = np.arange(0, fact, 0.2)
        levels_contf = np.arange(0, fact, 0.2)
    elif fact <= 10:
        levels_cont = np.arange(0,fact,0.5)
        levels_contf = np.arange(0,fact,0.5)
    elif fact <= 20:
        levels_cont = np.arange(0,fact,2)
        levels_contf = np.arange(0,fact,2)
    elif fact <= 50:
        levels_cont = np.arange(0,fact,5)
        levels_contf = np.arange(0,fact,5)
    else:
        levels_cont = np.arange(0,fact,20)
        levels_contf = np.arange(0,fact,20)
    levels_comp = np.linspace(0, fact_, 7)

    plt.subplot(3, 2, 1)
    plt.scatter(data[:, 0], data[:, 1], s=2, alpha=0.3)
    ax1 = plt.contour(X, Y, Z, levels=np.linspace(10, 20, 2))
    plt.colorbar(ax1,shrink=0.8)
    plt.xlim([x1_min,x1_max])
    plt.ylim([x2_min, x2_max])
    plt.title(var_name1 + var_name2 + ' (data), t=' + str(time) + ', z=' + str(z))
    plt.xlim([x1_min,x1_max])
    plt.ylim([x2_min, x2_max])
    plt.title(var_name1 + var_name2 + ' (data)')
    axis_label(var_name1, var_name2, amp)

    plt.subplot(3, 2, 3)
    ax1 = plt.hist2d(data[:, 0], data[:, 1], bins=30, normed=True, norm=colors.LogNorm())
    plt.colorbar(shrink=0.8)
    ax2 = plt.contour(X, Y, np.exp(Z), norm=colors.LogNorm(), linewidths=2, colors='w')
    axis_label(var_name1,var_name2, amp)
    plt.title('data histogram')
    plt.subplot(3, 2, 4)
    ax1 = plt.contourf(X, Y, np.exp(Z), norm=colors.LogNorm())
    plt.colorbar(ax1, shrink=0.8)
    plt.title('EM PDF')
    axis_label(var_name1, var_name2, amp)
    plt.subplot(3, 2, 5)
    plt.scatter(data[:, 0], data[:, 1], s=2, alpha=0.05)
    ax1 = plt.contour(X, Y, Z1, norm = colors.LogNorm(), linewidths=1.5)
    ax2 = plt.contour(X, Y, Z2, norm = colors.LogNorm(), linewidths=1.5)
    plt.plot([clf.means_[0, 0]], [clf.means_[0, 1]], 'wo', markersize=6)
    plt.plot([clf.means_[1, 0]], [clf.means_[1, 1]], 'wo', markersize=6)
    plt.colorbar(ax1, shrink=0.8)
    plt.title('f1, f2')
    plt.xlim([x1_min, x1_max])
    plt.ylim([x2_min, x2_max])
    axis_label(var_name1, var_name2, amp)
    plt.subplot(3, 2, 6)
    plt.scatter(data[:, 0], data[:, 1], s=2, alpha=0.05)
    ax1 = plt.contour(X, Y, np.exp(Z), norm = colors.LogNorm(), linewidths=1.5)
    plt.plot([clf.means_[0, 0]], [clf.means_[0, 1]], 'wo', markersize=6)
    plt.plot([clf.means_[1, 0]], [clf.means_[1, 1]], 'wo', markersize=6)
    plt.colorbar(ax1, shrink=0.8)
    plt.xlim([x1_min, x1_max])
    plt.ylim([x2_min, x2_max])
    axis_label(var_name1, var_name2, amp)
    plt.title('f = f1 + f2')

    fig.suptitle('EM2 PDF: ' + var_name1 + var_name2 + ' (t=' + str(time) + ', z=' + str(z) +'m)', fontsize=20)
    plt.savefig(os.path.join(
        fullpath_out,'EM2_bivar_figures','EM2_PDF_bivariate_' + var_name1 + '_' + var_name2 + '_' + str(time) + '_z' + str(
            np.int(z)) + 'm_log.png')
    )

    plt.close()
    return
#----------------------------------------------------------------------
def plot_PDF_samples_qt(data, var_name1, var_name2, clf, time, z):
    print('')
    amp = 100
    print('plot PDF samples qt: amp = '+np.str(amp))

    import matplotlib.mlab as mlab
    import matplotlib.cm as cm

    data_aux = np.ndarray(shape=((nx * ny), nvar))
    if var_name1 == 'qt':
        data_aux[:, 0] = data[:, 0] * amp
        data_aux[:, 1] = data[:, 1]
    elif var_name2 == 'qt':
        data_aux[:, 0] = data[:, 0]
        data_aux[:, 1] = data[:, 1] * amp

    clf_aux = mixture.GaussianMixture(n_components=2, covariance_type='full')
    clf_aux.fit(data_aux)

    det1 = np.linalg.det(clf.covariances_[0, :, :])
    det2 = np.linalg.det(clf.covariances_[1, :, :])
    # print(det1, det2)
    det_ = min(det1, det2)
    fact_ = 1. / np.sqrt((2 * np.pi) ** 2 * det_)
    fact = clf.weights_[0] * 1. / np.sqrt((2 * np.pi) ** 2 * det1) + clf.weights_[1] * 1. / np.sqrt(
        (2 * np.pi) ** 2 * det2)
    det1_aux = np.linalg.det(clf_aux.covariances_[0, :, :])
    det2_aux = np.linalg.det(clf_aux.covariances_[1, :, :])
    det_aux_ = min(det1_aux, det2_aux)
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
    print(np.amin(Z_aux), np.amax(Z_aux))

    mx1 = clf.means_[0, 0]
    my1 = clf.means_[0, 1]
    sx1 = np.sqrt(clf.covariances_[0, 0, 0])
    sy1 = np.sqrt(clf.covariances_[0, 1, 1])
    sxy1 = clf.covariances_[0, 1, 0]
    mx2 = clf.means_[1, 0]
    my2 = clf.means_[1, 1]
    sx2 = np.sqrt(clf.covariances_[1, 0, 0])
    sy2 = np.sqrt(clf.covariances_[1, 1, 1])
    sxy2 = clf.covariances_[1, 1, 0]
    Z1 = mlab.bivariate_normal(X, Y, sigmax=sx1, sigmay=sy1, mux=mx1, muy=my1, sigmaxy=sxy1)
    Z2 = mlab.bivariate_normal(X, Y, sigmax=sx2, sigmay=sy2, mux=mx2, muy=my2, sigmaxy=sxy2)
    mx1_aux = clf_aux.means_[0, 0]
    my1_aux = clf_aux.means_[0, 1]
    sx1_aux = np.sqrt(clf_aux.covariances_[0, 0, 0])
    sy1_aux = np.sqrt(clf_aux.covariances_[0, 1, 1])
    sxy1_aux = clf_aux.covariances_[0, 1, 0]
    mx2_aux = clf_aux.means_[1, 0]
    my2_aux = clf_aux.means_[1, 1]
    sx2_aux = np.sqrt(clf_aux.covariances_[1, 0, 0])
    sy2_aux = np.sqrt(clf_aux.covariances_[1, 1, 1])
    sxy2_aux = clf_aux.covariances_[1, 1, 0]
    Z1_aux = mlab.bivariate_normal(X_aux, Y_aux, sigmax=sx1_aux, sigmay=sy1_aux, mux=mx1_aux, muy=my1_aux, sigmaxy=sxy1_aux)
    Z2_aux = mlab.bivariate_normal(X_aux, Y_aux, sigmax=sx2_aux, sigmay=sy2_aux, mux=mx2_aux, muy=my2_aux, sigmaxy=sxy2_aux)

    fig = plt.figure(figsize=(12, 16))

    plt.subplot(4, 2, 1)
    plt.scatter(data[:, 0], data[:, 1], s=5, alpha=0.5)
    ax1 = plt.contour(X, Y, Z, levels=np.linspace(10, 20, 11))
    plt.colorbar(ax1,shrink=0.8)
    plt.xlim([x1_min,x1_max])
    plt.ylim([x2_min, x2_max])
    plt.title(var_name1 + var_name2 + ' (data), t=' + str(time) + ', z=' + str(z))
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)

    plt.subplot(4, 2, 3)
    ax1 = plt.contourf(X, Y, np.exp(Z))
    # plt.scatter(X, Y, s=2, alpha=0.5)
    plt.colorbar(ax1, shrink=0.8)
    plt.title('EM PDF')
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    plt.subplot(4, 2, 4)
    ax1 = plt.contourf(X_aux, Y_aux, np.exp(Z_aux))
    plt.colorbar(ax1, shrink=0.8)
    plt.title('EM PDF')
    if var_name1 == 'qt':
        plt.xlabel(var_name1 + '  (* ' + np.str(amp) + ')')
        plt.ylabel(var_name2)
    else:
        plt.xlabel(var_name1)
        plt.ylabel(var_name2 + '  (* ' + np.str(amp) + ')')

    plt.subplot(4, 2, 5)
    plt.scatter(data[:, 0], data[:, 1], s=5, alpha=0.3)
    levels = np.linspace(0, fact, 10)
    ax1 = plt.contour(X, Y, np.exp(Z), levels=levels, linewidths=2)
    plt.plot([clf.means_[0, 0]], [clf.means_[0, 1]], 'o', markersize=8)
    plt.plot([clf.means_[1, 0]], [clf.means_[1, 1]], 'o', markersize=8)
    plt.colorbar(ax1, shrink=0.8)
    plt.title('EM PDF')
    # plt.title(np.str(clf.means_))
    plt.xlim([x1_min, x1_max])
    plt.ylim([x2_min, x2_max])
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    plt.subplot(4, 2, 6)
    plt.scatter(data_aux[:, 0], data_aux[:, 1], s=5, alpha=0.3)
    # ax1 = plt.contour(X_aux, Y_aux, np.exp(Z_aux), levels=levels, linewidths=2)
    ax1 = plt.contour(X_aux, Y_aux, np.exp(Z_aux), linewidths=2)
    plt.plot([clf_aux.means_[0, 0]], [clf_aux.means_[0, 1]], 'o', markersize=8)
    plt.plot([clf_aux.means_[1, 0]], [clf_aux.means_[1, 1]], 'o', markersize=8)
    plt.colorbar(ax1, shrink=0.8)
    # plt.xlim([x1_min, x1_max])
    # plt.ylim([x2_min, x2_max])
    plt.title('EM PDF')
    if var_name1 == 'qt':
        plt.xlabel(var_name1 + '  (* ' + np.str(amp) + ')')
        plt.ylabel(var_name2)
    else:
        plt.xlabel(var_name1)
        plt.ylabel(var_name2 + '  (* ' + np.str(amp) + ')')

    plt.subplot(4, 2, 7)
    levels = np.linspace(0, fact_, 7)
    plt.scatter(data[:, 0], data[:, 1], s=2, alpha=0.2)
    ax1 = plt.contour(X, Y, Z1, levels=levels, linewidths=2)
    ax2 = plt.contour(X, Y, Z2, levels=levels, linewidths=2)
    plt.plot([clf.means_[0, 0]], [clf.means_[0, 1]], 'o', markersize=8)
    plt.plot([clf.means_[1, 0]], [clf.means_[1, 1]], 'o', markersize=8)
    # ax1 = plt.contour(X, Y, clf.weights_[0] * Z1 + clf.weights_[0] * Z2, levels=levels,linewidths=2)
    plt.colorbar(ax1, shrink=0.8)
    plt.title('f=f1+f2')
    plt.xlim([x1_min, x1_max])
    plt.ylim([x2_min, x2_max])
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    plt.subplot(4, 2, 8)
    levels = np.linspace(0, fact_aux_, 7)
    plt.scatter(data_aux[:, 0], data_aux[:, 1], s=2, alpha=0.2)
    ax1 = plt.contour(X_aux, Y_aux, Z1_aux, levels=levels, linewidths=2)
    ax2 = plt.contour(X_aux, Y_aux, Z2_aux, levels=levels, linewidths=2)
    plt.plot([clf_aux.means_[0, 0]], [clf_aux.means_[0, 1]], 'o', markersize=8)
    plt.plot([clf_aux.means_[1, 0]], [clf_aux.means_[1, 1]], 'o', markersize=8)
    # ax1 = plt.contour(X_aux, Y_aux, clf_aux.weights_[0] * Z1_aux + clf_aux.weights_[0] * Z2_aux, levels=levels,
    #                   linewidths=2)
    plt.colorbar(ax1, shrink=0.8)
    plt.title('f=f1+f2')
    if var_name1 == 'qt':
        plt.xlabel(var_name1 + '  (* ' + np.str(amp) + ')')
        plt.ylabel(var_name2)
    else:
        plt.xlabel(var_name1)
        plt.ylabel(var_name2 + '  (* ' + np.str(amp) + ')')

    fig.suptitle('EM2 PDF: '+var_name1+var_name2, fontsize=20)
    plt.savefig(os.path.join(
        fullpath_out,'EM2_bivar_figures','EM2_PDF_bivariate_' + var_name1 + '_' + var_name2 + '_' + str(time) + '_z' + str(
        np.int(z)) + 'm.png')
    )

    plt.close()
    return
#----------------------------------------------------------------------
def axis_label(var_name1, var_name2, amp):
    if var_name1 == 'qt':
        plt.xlabel(var_name1 + '  (* ' + np.str(amp) + ')')
        plt.ylabel(var_name2)
    elif var_name2 == 'qt':
        plt.xlabel(var_name1)
        plt.ylabel(var_name2 + '  (* ' + np.str(amp) + ')')
    else:
        plt.xlabel(var_name1)
        plt.ylabel(var_name2)
    return
#----------------------------------------------------------------------





# ____________________
def dump_pickle(data,out_path,file_name):
    data_ = (1.4,42)
    # output = open(os.path.join(out_path,'data.pkl'), 'w')
    output = open(os.path.join(out_path, file_name), 'w')
    pickle.dump(data, output)
    output.close()
    return
def test_pickle(in_path,file_name):
    print('')
    print('------- test pickle ------')
    fullpath_in = os.path.join(in_path,file_name)
    f = open(fullpath_in)
    data = pickle.load(f)
    print(data)
    print('')
    var = data['w']
    print(var)
    print
    ''
    means_ = var['means']
    print(means_)
    print('-------------------------')
    print('')
    return

# ____________________
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
    print('-------- dump variable --------', var_name, group_name)
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
    # print('add means: ', var_name, path)
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

def add_weights(path, var_name, ncomp, nvar):
    # print('add weights: ', var_name, path)
    # rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
    rootgrp = nc.Dataset(path, 'r+')
    group = rootgrp.groups['weights']
    var = group.createVariable(var_name, 'f8', ('nz', 'EM2'))
    rootgrp.close()
    return


def write_mean(path, group_name, data, var_name):
    # print('write mean:', path, var_name, data.shape)
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

def write_weights(path, group_name, data, var_name):
    # print('write weights:', path, var_name, data.shape)
    rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
    fieldgrp = rootgrp.groups[group_name]
    var = fieldgrp.variables[var_name]
    print(var.shape, data.shape)
    var[:, :] = data[:,:]
    rootgrp.close()
    return

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