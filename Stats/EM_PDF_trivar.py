import netCDF4 as nc
import argparse
import os
import sys
import json as  simplejson
import pylab as plt
from matplotlib.colors import LogNorm

import pickle

import numpy as np
import itertools
from scipy import linalg
from sklearn import mixture

sys.path.append("..")
from io_read_in_files import read_in_netcdf_fields


label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size


'''
Gaussian Mixture Model = Superposition of multiple Gaussian Distributions
(1) Fit A GMM to data points to get means, covariance matrices and relaive weights
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
    # ______________________
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
    zrange = np.arange(10, 11, 5)
    print('zrange', zrange)
    print('_______________________')
    if case_name == 'DCBLSoares':
        var_list = ['w','u','s']
    else:
        var_list = ['w','s','qt']
    # var_list = ['s']


    '''
    Tri - variate PDF for (s, qt, w)
    '''
    global nvar, ncomp
    ncomp = 2
    nvar = 3
    data = np.ndarray(shape=((nx * ny), nvar))
    means_ = np.ndarray(shape=(len(zrange), ncomp, nvar))
    covariance_ = np.zeros(shape=(len(zrange), ncomp, nvar, nvar))
    weights_ = np.zeros(shape=(len(zrange), 2))
    mean_tot = np.ndarray(shape=(len(zrange), nvar))
    covariance_tot = np.zeros(shape=(len(zrange), nvar, nvar))
    # ---
    for d in files:
        print('- - - - time: ' + str(d) + '- - - -')
        nc_file_name = 'EM2_trivar_' + str(d)
        create_statistics_file(os.path.join(fullpath_out, 'EM2_trivar'), nc_file_name, ncomp, nvar, len(zrange))

        fullpath_in = os.path.join(args.path, 'fields', d)
        print(fullpath_in)
        for n1 in range(len(var_list)):
            var1 = var_list[n1]
            # data1_ = read_in_netcdf_fields(var1, fullpath_in).reshape((nx * ny, nz))
            for n2 in range(n1, len(var_list)):
                var2 = var_list[n2]
                # data2_ = read_in_netcdf_fields(var2, fullpath_in).reshape((nx * ny, nz))
                for n3 in range(n2, len(var_list)):
                    var3 = var_list[n3]

                    var1 = 'w'
                    var2 = 'temperature'
                    var3 = 'qt'
                    data1_ = read_in_netcdf_fields(var1, fullpath_in).reshape((nx * ny, nz))
                    data2_ = read_in_netcdf_fields(var2, fullpath_in).reshape((nx * ny, nz))
                    data3_ = read_in_netcdf_fields(var3, fullpath_in).reshape((nx * ny, nz))
                    print(np.min(data2_)), np.max(data2_)
                    print('----', var1, var2, var3)
                    if var1 == var2 or var2 == var3 or var1 == var3:
                        continue
                    for i in range(len(zrange)):
                        iz = zrange[i]
                        data[:, 0] = data1_[:, iz]
                        data[:, 1] = data2_[:, iz]
                        data[:, 2] = data3_[:, iz]

                        clf = Gaussian_mixture_trivariate(data, var1, var2, var3, np.int(d[0:-3]), iz*dz)
                        means_[i, :, :] = clf.means_[:, :]
                        covariance_[i, :, :, :] = clf.covariances_[:, :, :]
                        weights_[i, :] = clf.weights_[:]

                        print('')
                        print('checking: ')
                        if var1 == 'w':
                            print('w max (z=', np.str(iz), ': ', np.amax(data1_[:, iz]), np.amax(data[:, 0]))
                        if var3 == 'qt':
                            print('qt max (z=', np.str(iz), ': ', np.amax(data3_[:, iz]), np.amax(data[:, 2]))

                        mean_tot[i, :], covariance_tot[i, :, :] = covariance_estimate_from_multicomp_pdf(clf)

        dump_variable(os.path.join(fullpath_out, 'EM2_trivar', nc_file_name), 'means', means_, var1+var2+var3,
                                      ncomp, nvar, len(zrange))
        dump_variable(os.path.join(fullpath_out, 'EM2_trivar', nc_file_name), 'covariances', covariance_,
                                      var1 + var2 + var3, ncomp, nvar, len(zrange))
        dump_variable(os.path.join(fullpath_out, 'EM2_trivar', nc_file_name), 'weights', weights_, var1+var2+var3,
                                      ncomp, nvar, len(zrange))
        dump_variable(os.path.join(fullpath_out, 'EM2_trivar', nc_file_name), 'mean_tot', mean_tot, var1+var2+var3,
                                    ncomp, nvar,len(zrange))
        dump_variable(os.path.join(fullpath_out, 'EM2_trivar', nc_file_name), 'covariances_tot', covariance_tot, var1+var2+var3,
                                    ncomp, nvar,len(zrange))

    return

#----------------------------------------------------------------------
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

#----------------------------------------------------------------------
# compute matrices
def compute_Zmatrix(data,var_name1,var_name2,var_name3,clf,amp_qt,amp_w,time,z):
    # data normal (without amplification works well for var1, var2)
    # data aux: ZZ seems to be out of range (for var1, var2)

    n_sample = 300
    x_ = np.linspace(np.amin(data[:, 0]), np.amax(data[:, 0]), n_sample)
    y_ = np.linspace(np.amin(data[:, 1]), np.amax(data[:, 1]), n_sample)
    z_ = np.linspace(np.amin(data[:, 2]), np.amax(data[:, 2]), n_sample)
    X, Y, Z = np.meshgrid(x_, y_, z_)
    XX = np.array([X.ravel(), Y.ravel(), Z.ravel()]).T

    print('Computing Z matrix')
    if var_name1 == 'w' and var_name3 == 'qt':
        print('max w: ', np.amax(data[:, 0]), np.amin(data[:, 0]), np.amax(x_), np.amin(x_), np.amax(X))
        print('max s: ', np.amax(data[:, 1]), np.amin(data[:, 1]), np.amax(y_), np.amin(y_), np.amax(Y))
        print('max qt: ', np.amax(data[:, 2]), np.amin(data[:, 2]), np.amax(z_), np.amin(z_), np.amax(Z))

    ZZ = clf.score_samples(XX).reshape(X.shape)
    print ''
    print('ZZ: ', ZZ.shape, XX.shape, X.shape, x_.shape)


    plt.figure()

    # plt.contourf(x_, y_, np.exp(ZZ[:, :, 0]))
    # plt.colorbar()
    # plt.contour(x_, y_, np.exp(ZZ[:, :, 100]), colors='k', linewidth=1)
    # plt.scatter(data[:, 0], data[:, 1], s=2, alpha=0.8)
    # plt.xlabel(var_name1 + ' *' + np.str(amp_w))
    # plt.ylabel(var_name2)
    # plt.savefig(os.path.join(
    #     fullpath_out, 'EM2_trivar_figures',
    #     'test_EM2_PDF_trivariate_' + var_name1 + var_name2 + '_' + str(time) + '_z' + str(
    #         np.int(z)) + 'm.png')
    # )
    #
    # plt.figure()
    # plt.contourf(x_, z_, np.exp(ZZ[:, 0, :]))
    # plt.colorbar()
    # plt.contour(x_, z_, np.exp(ZZ[:, 100, :]), colors='k', linewidth=2)
    # plt.scatter(data[:, 0], data[:, 2], s=2, alpha=0.8)
    # plt.xlabel(var_name1 + ' *'+np.str(amp_w))
    # plt.ylabel(var_name3 + ' *'+np.str(amp_qt))
    # plt.savefig(os.path.join(
    #     fullpath_out, 'EM2_trivar_figures',
    #     'test_EM2_PDF_trivariate_' + var_name1 + var_name3 + '_' + str(time) + '_z' + str(
    #         np.int(z)) + 'm.png')
    # )

    a = (np.sum(np.exp(ZZ),axis=2))
    print('----a', np.str(np.amax(np.exp(ZZ[:,:,100]))), np.str(np.amax(a)), ZZ.shape, a.shape)
    plt.figure()
    plt.contourf(x_, y_, a)
    plt.colorbar()
    # plt.contour(x_, z_, np.exp(ZZ[:, 100, :]), colors='k', linewidth=2)
    plt.scatter(data[:, 0], data[:, 1], s=2, alpha=0.8)
    plt.xlabel(var_name1 + ' *' + np.str(amp_w))
    plt.ylabel(var_name2)
    plt.savefig(os.path.join(
        fullpath_out, 'EM2_trivar_figures',
        'test_EM2_PDF_trivariate_' + var_name1 + var_name2 + '_' + str(time) + '_sum_z' + str(
            np.int(z)) + 'm.png')
    )

    a = (np.sum(np.exp(ZZ), axis=1))
    print('----a', np.str(np.amax(np.exp(ZZ[:, 100, :]))), np.str(np.amax(a)), ZZ.shape, a.shape)
    plt.figure()
    plt.contourf(x_, z_, a)
    plt.colorbar()
    # plt.contour(x_, z_, np.exp(ZZ[:, 100, :]), colors='k', linewidth=2)
    plt.scatter(data[:, 0], data[:, 2], s=2, alpha=0.8)
    plt.xlabel(var_name1 + ' *' + np.str(amp_w))
    plt.ylabel(var_name3 + ' *' + np.str(amp_qt))
    plt.savefig(os.path.join(
        fullpath_out, 'EM2_trivar_figures',
        'test_EM2_PDF_trivariate_' + var_name1 + var_name3 + '_' + str(time) + '_sum_z' + str(
            np.int(z)) + 'm.png')
    )

    a = (np.sum(np.exp(ZZ), axis=0))
    print('----a', np.str(np.amax(np.exp(ZZ[100, :, :]))), np.str(np.amax(a)), ZZ.shape, a.shape)
    plt.figure()
    plt.contourf(y_, z_, a)
    plt.colorbar()
    # plt.contour(x_, z_, np.exp(ZZ[:, 100, :]), colors='k', linewidth=2)
    plt.scatter(data[:, 1], data[:, 2], s=2, alpha=0.8)
    plt.xlabel(var_name2)
    plt.ylabel(var_name3+ ' *' + np.str(amp_qt))
    plt.savefig(os.path.join(
        fullpath_out, 'EM2_trivar_figures',
        'test_EM2_PDF_trivariate_' + var_name2 + var_name3 + '_' + str(time) + '_sum_z' + str(
            np.int(z)) + 'm.png')
    )
    # plt.show()

    return ZZ

#----------------------------------------------------------------------
def Gaussian_mixture_trivariate(data, var_name1, var_name2, var_name3, time, z):
    # !!! when assigning data_aux = data --> data is changed as soon as data_aux is changed! (no effect on other modules)
    # even when (1) data_ori = data, (2) data_aux = data, (3) data_aux *= 1e2 --> data, data_ori and data_aux are changed
    # Note: not given if element-wise (or column-wise) assignment: data_aux[:,1] = data[:,1]
    print('')
    print('Gaussian Mixture: ')

    amp_qt = 1e2
    amp_w = 1e0

    clf = mixture.GaussianMixture(n_components=2, covariance_type='full')
    clf.fit(data)

    data_aux = np.ndarray(shape=((nx * ny), nvar))
    data_aux = np.array(data,copy=True)
    if var_name1 == 'w':
        print('normalising w')
        data_aux[:,0] = data[:,0] * amp_w
    if var_name2 == 's':
        s_mean = np.average(data[:,1])
        print('normalising s', s_mean)
        data_aux[:,1] = data[:,1] - s_mean
    if var_name3 == 'qt':
        print('normalising qt')
        data_aux[:, 2] = data[:, 2] * amp_qt

    clf_aux = mixture.GaussianMixture(n_components=2, covariance_type='full')
    clf_aux.fit(data_aux)
    plot_PDF_samples(data_aux, var_name1, var_name2, var_name3, clf_aux, amp_qt, amp_w, time, z)
    return clf_aux

#----------------------------------------------------------------------
def plot_PDF_samples(data, var_name1, var_name2, var_name3, clf, amp_qt, amp_w, time, z):
    import matplotlib.mlab as mlab
    import matplotlib.cm as cm
    print('Plot: '+var_name1+var_name2+var_name3)

    print('in plotting max w: ', np.amax(data[:, 0]), np.amin(data[:,0]))
    print('in plotting max temp: ', np.amax(data[:, 1]), np.amin(data[:, 1]))
    print('in plotting max qt: ', np.amax(data[:, 2]), np.amin(data[:,2]))
    print('')

    # Plotting
    n_sample = 100
    x_ = np.linspace(np.amin(data[:,0]), np.amax(data[:,0]), n_sample)
    y_ = np.linspace(np.amin(data[:,1]), np.amax(data[:,1]), n_sample)
    z_ = np.linspace(np.amin(data[:,2]), np.amax(data[:,2]), n_sample)
    # x_ = np.linspace(0,4, n_sample)
    # y_ = np.linspace(10, 40, n_sample)
    # z_ = np.linspace(100, 400, n_sample)

    # X, Y, Z = np.meshgrid(x_, y_, z_)
    # XX = np.array([X.ravel(), Y.ravel(), Z.ravel()]).T
    # # ZZ = clf.score_samples(XX).reshape(X.shape)
    # ZZ = clf.score_samples(XX)
    # print(zrange.shape, nvar, ncomp, n_sample, data.shape)
    # print('plotting: ZZ', ZZ.shape, XX.shape)
    # ZZ = ZZ.reshape(X.shape)
    # print(ZZ.shape)

    XX_ = np.ndarray(shape=(n_sample**nvar,nvar))
    delta_i = n_sample*n_sample
    delta_j = n_sample
    for i in range(n_sample):
        for j in range(n_sample):
            for k in range(n_sample):
                shift = i*delta_i + j*delta_j + k
                XX_[shift, 0] = x_[i]
                XX_[shift, 1] = y_[j]
                XX_[shift, 2] = z_[k]
                # print('shift, i, j, k: ', shift, i, j, k)

    print('')
    # print x_
    # print y_
    # print z_
    # print('')
    # print(XX_.shape)
    # print(XX_)
    ZZ_ = clf.score_samples(XX_)
    # print('')
    # print('ZZ_: ', ZZ_.shape)

    ZZ = np.ndarray(shape=(n_sample, n_sample, n_sample))
    for k in range(n_sample**nvar):
        k_shift = np.mod(k,delta_j)
        j_shift = (np.mod(k,delta_i) - k_shift) / delta_j
        i_shift = (k - k_shift - n_sample*j_shift) / delta_i
        # print('k, ishift, jshift, kshift: ', k, i_shift, j_shift, k_shift)
        ZZ[i_shift, j_shift, k_shift] = ZZ_[k]

    # print('')
    # print('ZZ: ', ZZ.shape)
    # print('')
    # print(ZZ_)
    # print('')
    # print(ZZ)
    # print(np.exp(ZZ))

    fig = plt.figure(figsize=(18, 15))

    ZZ12 = np.sum(ZZ,axis=2)
    plt.subplot(3, 5, 1)
    plt.scatter(data[:, 0], data[:, 1], s=2, alpha=0.3)
    axis_label(var_name1, var_name2, amp_qt, amp_w)
    plt.subplot(3, 5, 2)
    ax1 = plt.hist2d(data[:, 0], data[:, 1], bins=30, normed=True)
    plt.colorbar(shrink=0.8)
    axis_label(var_name1, var_name2, amp_qt, amp_w)
    plt.title('data histogram')
    plt.subplot(3, 5, 3)
    # xxx
    ax1 = plt.contourf(x_, y_, np.sum(np.exp(ZZ),axis=2).T)
    # xxx
    plt.colorbar(ax1, shrink=0.8)
    axis_label(var_name1, var_name2, amp_qt, amp_w)
    plt.subplot(3, 5, 4)
    ax1 = plt.hist2d(data[:, 0], data[:, 1], bins=30, normed=True)
    # xxx
    ax2 = plt.contour(x_, y_, np.sum(np.exp(ZZ), axis=2).T, colors='w', linedwidht=2)
    # xxx
    # plt.colorbar(shrink=0.8)
    axis_label(var_name1, var_name2, amp_qt, amp_w)


    plt.subplot(3, 5, 6)
    plt.scatter(data[:, 0], data[:, 2], s=2, alpha=0.3)
    axis_label(var_name1, var_name3, amp_qt, amp_w)
    plt.subplot(3, 5, 7)
    ax1 = plt.hist2d(data[:, 0], data[:, 2], bins=30, normed=True)
    plt.colorbar(shrink=0.8)
    axis_label(var_name1, var_name3, amp_qt, amp_w)
    plt.title('data histogram')
    plt.subplot(3, 5, 8)
    # xxx
    ax1 = plt.contourf(x_, z_, np.sum(np.exp(ZZ), axis=1).T)
    # xxx
    plt.colorbar(ax1, shrink=0.8)
    axis_label(var_name1, var_name3, amp_qt, amp_w)
    plt.subplot(3, 5, 9)
    ax1 = plt.hist2d(data[:,0], data[:, 2], bins=30, normed=True)
    # xxx
    ax2 = plt.contour(x_, z_ , np.sum(np.exp(ZZ), axis=1).T, colors='w', linewidth=2)
    # xxx
    # plt.colorbar(ax1, shrink=0.8)
    axis_label(var_name1, var_name3, amp_qt, amp_w)

    plt.subplot(3, 5, 11)
    plt.scatter(data[:, 1], data[:, 2], s=2, alpha=0.3)
    axis_label(var_name2, var_name3, amp_qt, amp_w)
    plt.subplot(3, 5, 12)
    ax1 = plt.hist2d(data[:, 1], data[:, 2], bins=30, normed=True)
    plt.colorbar(shrink=0.8)
    plt.title('data histogram')
    plt.subplot(3, 5, 13)
    # xxx
    ax1 = plt.contourf(y_, z_, np.sum(np.exp(ZZ), axis=0).T)
    # xxx
    axis_label(var_name2, var_name3, amp_qt, amp_w)
    plt.colorbar(ax1, shrink=0.8)
    plt.subplot(3, 5, 14)
    ax1 = plt.hist2d(data[:,1], data[:,2], bins=30, normed=True)
    # xxx
    ax2 = plt.contour(y_, z_, np.sum(np.exp(ZZ), axis=0).T, colors='w', linewidth=2)
    # xxx
    axis_label(var_name2, var_name3, amp_qt, amp_w)
    # plt.colorbar(ax1, shrink=0.8)


    fig.suptitle('EM2 PDF: ' + var_name1 + var_name2 + var_name3 + ' (t=' + str(time) + ', z=' + str(z) +'m)', fontsize=20)
    plt.savefig(os.path.join(
        fullpath_out,'EM2_trivar_figures','EM2_PDF_trivariate_' + var_name1 + var_name2 + var_name3 + '_' + str(time) + '_z' + str(
            np.int(z)) + 'm.png')
    )

    plt.close()



    # # YY = XX.T
    # # print(XX[:,0].reshape(X.shape))
    # # Z = mlab.bivariate_normal(X, Y, sigmax=sx[0], sigmay=sy[0], mux=mu1[0], muy=mu1[1],sigmaxy=sxy[0])
    #
    # # levels_tot = np.linspace(0, fact, 10)
    # # if fact <= 2:
    # #     levels_cont = np.arange(0, fact, 0.2)
    # #     levels_contf = np.arange(0, fact, 0.2)
    # # elif fact <= 10:
    # #     levels_cont = np.arange(0,fact,0.5)
    # #     levels_contf = np.arange(0,fact,0.5)
    # # elif fact <= 20:
    # #     levels_cont = np.arange(0,fact,2)
    # #     levels_contf = np.arange(0,fact,2)
    # #     levels_contf = np.arange(0,fact,2)
    # # elif fact <= 50:
    # #     levels_cont = np.arange(0,fact,5)
    # #     levels_contf = np.arange(0,fact,5)
    # # else:
    # #     levels_cont = np.arange(0,fact,20)
    # #     levels_contf = np.arange(0,fact,20)
    # # levels_comp = np.linspace(0, fact_, 7)
    #
    # plt.subplot(3, 5, 1)
    # plt.scatter(data[:, 0], data[:, 1], s=2, alpha=0.3)
    # plt.title(var_name1 + var_name2 + ' (data)')
    # axis_label(var_name1, var_name2, amp_qt, amp_w)

    #
    #
    # plt.subplot(3, 5, 6)
    # plt.scatter(data[:, 0], data[:, 2], s=2, alpha=0.3)
    # plt.title(var_name1 + var_name3 + ' (data)')
    # axis_label(var_name1, var_name3, amp_qt, amp_w)
    # plt.subplot(3, 5, 7)
    # ax1 = plt.hist2d(data[:, 0], data[:, 2], bins=30, normed=True)
    # plt.colorbar(shrink=0.8)
    # # ax2 = plt.contour(X, Y, np.exp(Z), levels=levels_cont, linewidths=1, colors='w')
    # axis_label(var_name1, var_name3, amp_qt, amp_w)
    # plt.title('data histogram')
    # plt.subplot(3, 5, 9)
    # ax1 = plt.contourf(x_, z_, np.sum(np.exp(ZZ).T, axis=1))
    # plt.subplot(3, 5, 10)
    # # ax1 = plt.contourf(X, Y, np.exp(Z), levels=levels_contf)
    # # ax1 = plt.contourf(X[:, 0, :], Y[:, 0, :], np.exp(ZZ[:, 0, :]))
    # ax1 = plt.contourf(x_, z_, np.sum(np.exp(ZZ),axis=1))
    # plt.colorbar(ax1, shrink=0.8)
    # plt.title('EM PDF')
    # axis_label(var_name1, var_name3, amp_qt, amp_w)
    #
    # plt.subplot(3, 5, 11)
    # plt.scatter(data[:, 1], data[:, 2], s=2, alpha=0.3)
    # plt.title(var_name2 + var_name3 + ' (data)')
    # axis_label(var_name2, var_name3, amp_qt, amp_w)
    # plt.subplot(3, 5, 12)
    # ax1 = plt.hist2d(data[:, 1], data[:, 2], bins=30, normed=True)
    # plt.colorbar(shrink=0.8)
    # # ax2 = plt.contour(X, Y, np.exp(Z), levels=levels_cont, linewidths=1, colors='w')
    # axis_label(var_name2,var_name3, amp_qt, amp_w)
    # plt.title('data histogram')
    # plt.subplot(3, 5, 14)
    # ax1 = plt.contourf(y_, z_, np.sum(np.exp(ZZ).T, axis=0))
    # axis_label(var_name2, var_name3, amp_qt, amp_w)
    # plt.subplot(3, 5, 15)
    # # ax1 = plt.contourf(X, Y, np.exp(Z), levels=levels_contf)
    # # ax1 = plt.contourf(X[0, :, :], Y[0, :, :], np.exp(ZZ[0, :, :]))
    # ax1 = plt.contourf(y_, z_, np.sum(np.exp(ZZ),axis=0))
    # plt.colorbar(ax1, shrink=0.8)
    # plt.title('EM PDF')

    return



# ----------------------------------------------------------------------
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
# ____________________

def create_statistics_file(path,file_name, ncomp, nvar, nz_):
    # ncomp: number of Gaussian components in EM
    # nvar: number of variables of multi-variate Gaussian components
    print('create file:', path, file_name)
    rootgrp = nc.Dataset(os.path.join(path,file_name), 'w', format='NETCDF4')
    dimgrp = rootgrp.createGroup('dims')
    means_grp = rootgrp.createGroup('means')
    cov_grp = rootgrp.createGroup('covariances')
    weights_grp = rootgrp.createGroup('weights')
    mean_tot_grp = rootgrp.createGroup('mean_tot')
    cov_tot_grp = rootgrp.createGroup('covariances_tot')
    means_grp.createDimension('nz', nz_)
    means_grp.createDimension('ncomp', ncomp)
    means_grp.createDimension('nvar', nvar)
    cov_grp.createDimension('nz', nz_)
    cov_grp.createDimension('ncomp', ncomp)
    cov_grp.createDimension('nvar', nvar)
    weights_grp.createDimension('nz', nz_)
    weights_grp.createDimension('EM2', 2)
    mean_tot_grp.createDimension('nz', nz_)
    mean_tot_grp.createDimension('nvar', nvar)
    cov_tot_grp.createDimension('nz', nz_)
    cov_tot_grp.createDimension('nvar', nvar)
    ts_grp = rootgrp.createGroup('time')
    ts_grp.createDimension('nt', len(time) - 1)
    var = ts_grp.createVariable('t', 'f8', ('nt'))
    for i in range(len(time) - 1):
        var[i] = time[i + 1]
    z_grp = rootgrp.createGroup('z-profile')
    z_grp.createDimension('nz', len(zrange))
    var = z_grp.createVariable('height', 'f8', ('nz'))
    rootgrp.close()
    print('create file end')
    return

def dump_variable(path, group_name, data_, var_name, ncomp, nvar, nz_):
    print('--------')
    print('dump variable', path, group_name, var_name, data_.shape, ncomp, nvar)
    if group_name == 'means':
        add_means(path, var_name)
        data = np.empty((nz_,ncomp,nvar), dtype=np.double, order='c')
        for i in range(nz_):
            for j in range(ncomp):
                for k in range(nvar):
                    data[i,j,k] = data_[i,j,k]
        write_mean(path, group_name, data, var_name)

    elif group_name == 'covariances':
        add_covariance(path, var_name)
        data = np.empty((nz_, ncomp, nvar, nvar), dtype=np.double, order='c')
        for i in range(nz_):
            for j in range(ncomp):
                for k1 in range(nvar):
                    for k2 in range(nvar):
                        data[i, j, k1, k2] = data_[i, j, k1, k2]
        write_covar(path, group_name, data, var_name)

    elif group_name == 'weights':
        print('dump weights')
        add_weights(path, var_name)
        data = np.empty((nz_, ncomp), dtype=np.double, order='c')
        for i in range(nz_):
            for j in range(ncomp):
                data[i, j] = data_[i, j]
        write_weights(path, group_name, data, var_name)

    elif group_name == 'mean_tot':
        print('dump mean_tot')
        add_mean_tot(path, var_name)
        data = np.empty((nz_,nvar), dtype=np.double, order='c')
        for k in range(nz_):
            for i in range(nvar):
                data[k, i] = data_[k, i]
        write_mean_tot(path, group_name, data, var_name)

    elif group_name == 'covariances_tot':
        print('dump covars_tot')
        add_covariance_tot(path, var_name)
        data = np.empty((nz_,nvar,nvar), dtype=np.double, order='c')
        for k in range(nz_):
            for i1 in range(nvar):
                for i2 in range(nvar):
                    data[k, i1,i2] = data_[k, i1,i2]
        write_covar_tot(path, group_name, data, var_name)

    print('--------')
    return


def add_means(path, var_name):
    print('add means: ', var_name, path)
    # rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
    rootgrp = nc.Dataset(path, 'r+')
    group = rootgrp.groups['means']
    var = group.createVariable(var_name, 'f8', ('nz', 'ncomp', 'nvar'))
    rootgrp.close()
    return

def add_covariance(path, var_name):
    print('add covariance: ', var_name, path)
    # rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
    rootgrp = nc.Dataset(path, 'r+')
    group = rootgrp.groups['covariances']
    var = group.createVariable(var_name, 'f8', ('nz', 'ncomp', 'nvar', 'nvar'))
    rootgrp.close()
    return

def add_weights(path, var_name):
    print('add weights: ', var_name, path)
    # rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
    rootgrp = nc.Dataset(path, 'r+')
    group = rootgrp.groups['weights']
    var = group.createVariable(var_name, 'f8', ('nz', 'EM2'))
    rootgrp.close()
    return

def add_mean_tot(path, var_name):
    print('add mean tot: ', var_name, path)
    # rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
    rootgrp = nc.Dataset(path, 'r+')
    group = rootgrp.groups['mean_tot']
    var = group.createVariable(var_name, 'f8', ('nz', 'nvar'))
    rootgrp.close()
    return

def add_covariance_tot(path, var_name):
    print('add covariance: ', var_name, path)
    # rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
    rootgrp = nc.Dataset(path, 'r+')
    group = rootgrp.groups['covariances_tot']
    var = group.createVariable(var_name, 'f8', ('nz', 'nvar', 'nvar'))
    rootgrp.close()
    return


def write_mean(path, group_name, data, var_name):
    print('write mean:', path, var_name, data.shape)
    rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
    fieldgrp = rootgrp.groups[group_name]
    var = fieldgrp.variables[var_name]
    var[:, :, :] = data[:,:,:]
    rootgrp.close()
    return

def write_covar(path, group_name, data, var_name):
    print('write covar:', path, var_name, data.shape)
    rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
    fieldgrp = rootgrp.groups[group_name]
    var = fieldgrp.variables[var_name]
    var[:, :, :,:] = data[:, :, :,:]
    rootgrp.close()
    return

def write_weights(path, group_name, data, var_name):
    print('write weights:', path, var_name, data.shape)
    rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
    fieldgrp = rootgrp.groups[group_name]
    var = fieldgrp.variables[var_name]
    print(var.shape, data.shape)
    var[:, :] = data[:,:]
    rootgrp.close()
    return

def write_mean_tot(path, group_name, data, var_name):
    print('write mean tot:', path, var_name, data.shape)
    rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
    fieldgrp = rootgrp.groups[group_name]
    var = fieldgrp.variables[var_name]
    var[:, :] = data[:,:]
    rootgrp.close()
    return

def write_covar_tot(path, group_name, data, var_name):
    print('write covar tot:', path, var_name, data.shape)
    rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
    fieldgrp = rootgrp.groups[group_name]
    var = fieldgrp.variables[var_name]
    var[:, :,:] = data[:, :,:]
    rootgrp.close()
    return

def write_field(path, group_name, data, var_name):
    print('')
    print('write field:', path, var_name, data.shape, group_name)
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
def read_in_netcdf_fields(variable_name, fullpath_in):
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups['fields'].variables[variable_name]

    shape = var.shape
    # print('shape:',var.shape)
    data = np.ndarray(shape = var.shape)
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