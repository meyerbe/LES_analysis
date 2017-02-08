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
from sklearn import mixture, preprocessing

sys.path.append("..")
from io_read_in_files import read_in_netcdf

from thermodynamics import sat_adj, sat_adj_theta
from thermodynamics import theta_li


label_size = 10
plt.rcParams['xtick.direction']='out'
plt.rcParams['ytick.direction']='out'
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['lines.linewidth'] = 2

'''
Using sklearn.preprocessing.StandardScaler to normalise input data onto normally distributed data (mean zero and standard deviation between 0 and 1)


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
    zrange = np.arange(15, 30, 5)
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

    for d in files:
        p0 = 1e5
        # T, ql, qi = sat_adj(p0, 6500, 1e-3)
        nc_file_name = str(d)
        fullpath_in = os.path.join(in_path, 'fields', nc_file_name)
        print('fullpath_in', fullpath_in)
        nc_file_name = 'CC_s_' + str(d)
        # create_statistics_file(os.path.join(fullpath_out, 'CloudClosure'), nc_file_name, ncomp, nvar, len(zrange))



        '''(0) Read in Data & compute liquid potential temperature from temperature and moisture'''
        s = read_in_netcdf('s', 'fields', fullpath_in).reshape((nx * ny), nz)
        T = read_in_netcdf('temperature', 'fields', fullpath_in).reshape((nx * ny), nz)
        qt = read_in_netcdf('qt', 'fields', fullpath_in).reshape((nx * ny), nz)
        ql = read_in_netcdf('ql', 'fields', fullpath_in).reshape((nx * ny), nz)
        qi = np.zeros(shape=T.shape)
        theta_l = theta_li(p0,T,qt,ql,qi).reshape((nx * ny), nz)

        data_th = np.ndarray(shape=((nx * ny), nvar))
        data_s = np.ndarray(shape=((nx * ny), nvar))

        for i in range(len(zrange)):
            iz = zrange[i]
            data_th[:, 0] = theta_l[:, iz]
            data_th[:, 1] = qt[:, iz]
            data_s[:, 0] = s[:, iz]
            data_s[:, 1] = qt[:, iz]
            print('!!! zi ql: ', np.amax(ql[:,iz]), iz*dz, d)

            '''(1) Normalise Data (zero mean and normalised standard deviation)'''
            data_th_norm = preprocessing.StandardScaler().fit_transform(data_th)
            data_s_norm = preprocessing.StandardScaler().fit_transform(data_s)
            # X_test = X_scaler.transform(X_test)
            plot_data_comp(data_th, data_s, data_th_norm, data_s_norm, ql[:,iz], np.int(d[0:-3]), iz * dz)

            '''(2) Compute bivariate Gaussian PDF (theta_l, qt) '''
            # means, covariance, weights = Gaussian_mixture_bivariate(data, var1, var2, np.int(d[0:-3]), iz*dz)
            clf_th = Gaussian_bivariate(data_th, 'theta_l', 'qt', np.int(d[0:-3]), iz * dz)
            clf_th_norm = Gaussian_bivariate(data_th_norm, 'theta_l', 'qt', np.int(d[0:-3]), iz * dz)
            plot_PDF_samples_qt(data_th, data_th_norm, 'theta_l', 'qt', clf_th, clf_th_norm, np.int(d[0:-3]), iz * dz)
            # plot_PDF_samples(data_th_norm, data_th_norm, 'theta_l', 'qt', clf_th_norm, clf_th_norm, np.int(d[0:-3]), iz * dz)
            clf_s = Gaussian_bivariate(data_s, 's', 'qt', np.int(d[0:-3]), iz * dz)
            clf_s_norm = Gaussian_bivariate(data_s_norm, 's', 'qt', np.int(d[0:-3]), iz * dz)
            plot_PDF_samples_qt(data_s, data_s_norm, 's', 'qt', clf_s, clf_s_norm, np.int(d[0:-3]), iz * dz)
            # plot_PDF_samples(data_s_norm, data_s_norm, 's', 'qt', clf_s_norm, clf_s_norm, np.int(d[0:-3]), iz * dz)
            # plot_both_PDF(data_th_norm, data_s_norm, clf_th_norm, clf_s_norm, np.int(d[0:-3]), iz * dz)

            # means_[i, :, :] = clf.means_[:, :]
            # covariance_[i,:,:,:] = clf.covariances_[:,:,:]

            '''(3) Compute Kernel-Estimate PDF '''
            # kde, kde_aux = Kernel_density_estimate(data, 'T', 'qt', np.int(d[0:-3]), iz * dz)

            '''(4) Compute Relative Entropy '''
            # relative_entropy(data, clf, kde)

            '''(5) Compute Liquid Water '''
            nn = np.int(1e2)
            S, y = clf_s.sample(n_samples=nn)
            S_norm, y = clf_s_norm.sample(n_samples=nn)
            Th, y = clf_th.sample(n_samples=nn)
            Th_norm, y = clf_th_norm.sample(n_samples=nn)
            print('....', S_norm.shape, y.shape)
            print(S[:,0])
            print(S[:,1])
            print('')

            ql_computed = np.zeros(nn)
            for i in range(nn):
                if S[i,1] <= 0.016:
                    print(S[i,0], S[i,1])
                    ql_computed[i] = compute_ql(p0, S[i,0], Th[i, 0], S[i,1])


            # print(ql_computed)

            if np.amax(ql_computed) > 0:
                plt.figure()
                lim = 1e-5
                if np.amax(ql) > 0.0:
                    plt.hist(ql_computed, bins=50, range=(lim, np.amax(ql)))
                else:
                    plt.hist(ql_computed, normed=True)
                plt.suptitle('ql computed ('+str(iz*dz)+', '+str(np.int(d[0:-3])))
                plt.show()
            else:
                print('unsat')

        '''(4) Save Gaussian Mixture PDFs '''
        # dump_variable(os.path.join(fullpath_out, 'CloudClosure', nc_file_name), 'means', means_, 'qtT', ncomp, nvar, len(zrange))
        # dump_variable(os.path.join(fullpath_out, 'CloudClosure', nc_file_name), 'covariances', covariance_, 'qtT', ncomp, nvar, len(zrange))

    return


#----------------------------------------------------------------------
def compute_ql(p0, s, theta_l, qt):
    print('Compute ql: '+str(s)+', '+str(qt))

    ql = np.array(qt, copy=True)
    # T, ql = sat_adj_theta(p0, theta_l, qt)

    T, ql, qi = sat_adj(p0, s, qt)
    # sat_adj(p0, s, qt)


    return ql
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
    clf.fit(data)
    print('')

    return clf
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

    fig = plt.figure(figsize=(12, 18))
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

    plt.subplot(3, 2, 3)
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

    plt.subplot(3, 2, 4)
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

    plt.subplot(3, 2, 5)
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

    plt.subplot(3, 2, 6)
    bw = 1e-3
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
    plt.savefig(os.path.join(fullpath_out,'CloudClosure_figures_norm','CC_' + var_name1 + '_' + var_name2 + '_' + str(
        time) + '_z' + str(np.int(z)) + 'm_KDE.png'))
    plt.close()

    print('KDE shapes: ', kde.score_samples(XX).shape, X.shape)
    print(kde.get_params())

    return kde, kde
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
def plot_both_PDF(data1, data2, clf1, clf2, time, z):
    n_sample = 300
    x1 = np.linspace(np.amin(data1[:, 0]), np.amax(data1[:, 0]), n_sample)
    y1 = np.linspace(np.amin(data1[:, 1]), np.amax(data1[:, 1]), n_sample)
    x2 = np.linspace(np.amin(data2[:, 0]), np.amax(data2[:, 0]), n_sample)
    y2 = np.linspace(np.amin(data2[:, 1]), np.amax(data2[:, 1]), n_sample)

    XX1_ = np.ndarray(shape=(n_sample ** nvar, nvar))
    XX2_ = np.ndarray(shape=(n_sample ** nvar, nvar))
    for i in range(n_sample):
        for j in range(n_sample):
            shift = i * n_sample + j
            XX1_[shift, 0] = x1[i]
            XX1_[shift, 1] = y1[j]
            XX2_[shift, 0] = x2[i]
            XX2_[shift, 1] = y2[j]

    Z1_ = clf1.score_samples(XX1_)
    Z2_ = clf2.score_samples(XX2_)
    Z1 = np.ndarray(shape=(n_sample, n_sample))
    Z2 = np.ndarray(shape=(n_sample, n_sample))
    for k in range(n_sample**nvar):
        j_shift = np.mod(k, n_sample)
        i_shift = (k - j_shift) / n_sample
        Z1[i_shift, j_shift] = Z1_[k]
        Z2[i_shift, j_shift] = Z2_[k]

    fig = plt.figure(figsize=(12,18))
    plt.subplot(3, 2, 1)
    plt.scatter(data1[:,0], data1[:,1], alpha=0.2, s=5)
    ax1 = plt.contour(x1, y1, Z1)
    plt.colorbar(ax1)
    plt.title('data1')
    labeling('th_l', 'qt', 0)
    plt.subplot(3, 2, 2)
    plt.scatter(data2[:, 0], data2[:, 1], alpha=0.2, s=5)
    ax1 = plt.contour(x2, y2, Z2)
    plt.colorbar(ax1)
    plt.title('data2')
    labeling('s', 'qt', 0)
    plt.subplot(3,2,3)
    ax1 = plt.contourf(x1,y1,Z1)
    plt.colorbar(ax1)
    plt.title('Z1')
    plt.subplot(3, 2, 4)
    ax1 = plt.contourf(x2, y2, Z2)
    plt.title('Z2')
    plt.colorbar(ax1)
    plt.subplot(3, 2, 5)
    plt.contourf(x1, y1, Z1)
    ax1 = plt.contour(x2, y2,  Z2)
    plt.colorbar(ax1)
    plt.subplot(3, 2, 6)
    ax1 = plt.contour(x1, y1, Z1)
    plt.contourf(x2, y2, Z2)
    plt.colorbar(ax1)
    fig.suptitle('')
    plt.savefig(
        os.path.join(
            fullpath_out, 'CloudClosure_figures_norm/CC_comp_' + str(time) + '_z' + str(np.int(z)) + 'm.png')
    )
    plt.close()
    return
#--------------------------
def plot_data_comp(data_th, data_s, data_th_norm, data_s_norm, ql, time, z):
    plt.figure(figsize=(12,18))
    plt.subplot(4, 2, 1)
    lim = 1e-5
    if np.amax(ql) > 0.0:
        plt.hist(ql[:], bins=50, range=(lim, np.amax(ql)))
    else:
        plt.hist(ql, normed=True)
    plt.title(r'$q_l$')
    plt.subplot(4, 2, 2)
    plt.scatter(data_th[:,1], ql[:], alpha=0.2, s=5)
    plt.xlim([np.amin(data_th[:,1]),np.amax(data_th[:,1])])
    plt.ylim([-0.0001, np.amax(ql)])
    plt.plot([np.amin(data_th[:,1]),np.amax(data_th[:,1])],[lim, lim],'k',linewidth=1)
    lim = 0.016
    plt.plot([lim, lim], [np.amin(ql), np.amax(ql)], 'r', linewidth=2)
    labeling('qt', 'ql', 1)
    plt.title(r'$q_l$')

    plt.subplot(4, 2, 3)
    plt.scatter(data_th[:,0], data_th[:,1], alpha=0.2, s=5)
    plt.title(r'$\theta$')
    labeling('th_l', 'qt', 1)
    plt.subplot(4, 2, 4)
    plt.scatter(data_th_norm[:, 0], data_th_norm[:, 1], alpha=0.2, s=5)
    labeling('th_l', 'qt', 1)
    plt.title(r'$\theta$ normalised')
    plt.subplot(4, 2, 5)
    plt.scatter(data_s[:,0],data_s[:,1], alpha=0.2, s=5)
    plt.title(r'$s$')
    labeling('s', 'qt', 1)
    plt.subplot(4, 2, 6)
    plt.scatter(data_s_norm[:, 0], data_s_norm[:, 1], alpha=0.2, s=5)
    plt.title(r'$s$ normalised')
    labeling('s', 'qt', 1)
    plt.subplot(4, 2, 7)
    ax1 = plt.contourf(data_th[:,0].reshape(nx,ny))
    ax2 = plt.contour(data_th_norm[:, 0].reshape(nx, ny), colors = 'k', linewidths=1)
    plt.colorbar(ax1)
    plt.colorbar(ax2)
    plt.title(r'$\theta$')
    plt.subplot(4, 2, 8)
    ax1 = plt.contourf(data_s[:, 0].reshape(nx, ny))
    ax2 = plt.contour(data_s_norm[:, 0].reshape(nx, ny), colors='k', linewidths=1)
    plt.colorbar(ax1)
    plt.colorbar(ax2)
    plt.title(r'$s$')
    plt.savefig(
        os.path.join(
            fullpath_out, 'CloudClosure_figures_norm/data_' + str(time) + '_z' + str(np.int(z)) + 'm.png')
    )
    plt.close()
    return
#--------------------------
def plot_PDF_samples_qt(data, data_norm, var_name1, var_name2, clf, clf_norm, time, z):
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
    x = np.linspace(np.amin(data[:, 0]), np.amax(data[:, 0]), n_sample)
    y = np.linspace(np.amin(data[:, 1]), np.amax(data[:, 1]), n_sample)
    x_norm = np.linspace(np.amin(data_norm[:, 0]), np.amax(data_norm[:, 0]), n_sample)
    y_norm = np.linspace(np.amin(data_norm[:, 1]), np.amax(data_norm[:, 1]), n_sample)
    # x = np.linspace(x1_min, x1_max, n_sample)
    # y = np.linspace(x2_min, x2_max, n_sample)
    X, Y = np.meshgrid(x, y)
    XX = np.array([X.ravel(), Y.ravel()]).T
    Z = clf.score_samples(XX).reshape(X.shape)

    x_aux = np.linspace(np.amin(data_aux[:,0]), np.amax(data_aux[:,0]), n_sample)
    y_aux = np.linspace(np.amin(data_aux[:, 1]), np.amax(data_aux[:, 1]), n_sample)
    X_aux, Y_aux = np.meshgrid(x_aux, y_aux)
    XX_aux = np.array([X_aux.ravel(), Y_aux.ravel()]).T
    Z_aux = clf_aux.score_samples(XX_aux).reshape(X_aux.shape)

    X, Y = np.meshgrid(x_norm, y_norm)
    XX = np.array([X.ravel(), Y.ravel()]).T
    Z_norm = clf_norm.score_samples(XX).reshape(X.shape)

    fig = plt.figure(figsize=(12, 18))
    plt.subplot(3, 2, 1)
    plt.scatter(data_aux[:, 0], data_aux[:, 1], s=5, alpha=0.2)
    ax1 = plt.contour(X, Y, Z, colors='w', levels=np.linspace(10, 20, 11))
    plt.colorbar(ax1,shrink=0.8)
    plt.xlim([np.amin(data_aux[:, 0]), np.amax(data_aux[:, 0])])
    plt.ylim([np.amin(data_aux[:, 1]), np.amax(data_aux[:, 1])])
    plt.title(var_name1 + var_name2 + ' (data)')
    labeling(var_name1, var_name2, amp)
    plt.subplot(3, 2, 2)
    plt.scatter(data_norm[:, 0], data_norm[:, 1], s=5, alpha=0.2)
    ax1 = plt.contour(X_aux, Y_aux, Z_aux, colors='w', levels=np.linspace(10, 20, 11))
    plt.colorbar(ax1, shrink=0.8)
    plt.xlim([np.amin(data_norm[:,0]), np.amax(data_norm[:,0])])
    plt.ylim([np.amin(data_norm[:, 1]), np.amax(data_norm[:, 1])])
    plt.title(var_name1 + var_name2 + ' (data normalised)')
    labeling(var_name1, var_name2, 1)

    plt.subplot(3, 2, 3)
    plt.scatter(data_aux[:, 0], data_aux[:, 1], s=5, alpha=0.2)
    levels = np.linspace(0, fact_, 10)
    ax1 = plt.contour(X_aux, Y_aux, np.exp(Z_aux), linewidths=2)
    plt.plot([clf_aux.means_[0, 0]], [clf_aux.means_[0, 1]], 'o', markersize=8)
    plt.colorbar(ax1, shrink=0.8)
    plt.title('EM PDF')
    labeling(var_name1, var_name2, amp)
    plt.subplot(3, 2, 4)
    plt.scatter(data_norm[:, 0], data_norm[:, 1], s=5, alpha=0.2)
    ax1 = plt.contour(x_norm, y_norm, np.exp(Z_norm), linewidths=2)
    plt.plot([clf_norm.means_[0, 0]], [clf_norm.means_[0, 1]], 'o', markersize=8)
    plt.colorbar(ax1, shrink=0.8)
    plt.title('EM PDF')
    labeling(var_name1, var_name2, 1)

    plt.subplot(3, 2, 5)
    ax1 = plt.contourf(x_aux, y_aux, np.exp(Z_aux))
    # plt.scatter(X, Y, s=2, alpha=0.5)
    plt.colorbar(ax1, shrink=0.8)
    plt.title('EM PDF')
    labeling(var_name1, var_name2, amp)
    plt.subplot(3, 2, 6)
    ax1 = plt.contourf(x_norm, y_norm, np.exp(Z_norm))
    plt.colorbar(ax1, shrink=0.8)
    plt.title('EM PDF')
    labeling(var_name1, var_name2, 1)

    fig.suptitle('Cloud Closure: Bivariate Gaussian PDF fit (t=' + str(time) + ', z=' + str(z)+')', fontsize=20)
    plt.savefig(
        os.path.join(
            fullpath_out,'CloudClosure_figures_norm/CC_' + var_name1 + '_' + var_name2 + '_' + str(time) + '_z'
                         + str(np.int(z)) + 'm_bivariate.png')
    )

    plt.close()
    return
#--------------------------
def plot_PDF_samples(data, data_norm, var_name1, var_name2, clf, clf_norm, time, z):
    print('pdf samples', var_name1, var_name2)
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

    fig = plt.figure(figsize=(12, 18))
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
    plt.title(var_name1 + var_name2 + ' (data)')
    # plt.subplot(3, 2, 2)
    # plt.scatter(data_norm[:, 0], data_norm[:, 1], s=2, alpha=0.05)
    # ax1 = plt.contour(X, Y, Z, levels=np.linspace(10, 20, 2))
    # plt.colorbar(ax1, shrink=0.8)
    # plt.xlim([x1_min, x1_max])
    # plt.ylim([x2_min, x2_max])
    # plt.title(var_name1 + var_name2 + ' (data normalised)')

    plt.subplot(3, 2, 3)
    ax1 = plt.hist2d(data[:, 0], data[:, 1], bins=30, normed=True)
    plt.colorbar(shrink=0.8)
    ax2 = plt.contour(X, Y, np.exp(Z), levels=levels_cont, linewidths=1, colors='w')
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    plt.title('data histogram')
    plt.subplot(3, 2, 4)
    # ax1 = plt.contourf(X, Y, np.exp(Z),levels=levels_contf)
    ax1 = plt.contourf(X, Y, np.exp(Z))
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
    # ax1 = plt.contour(X, Y, np.exp(Z), levels=levels_cont, linewidths=1.5)
    ax1 = plt.contour(X, Y, np.exp(Z))
    # ax1 = plt.contour(X, Y, Z1+Z2, linewidths=1.5)
    plt.plot([clf.means_[0, 0]], [clf.means_[0, 1]], 'wo', markersize=6)
    plt.colorbar(ax1, shrink=0.8)
    plt.xlim([x1_min, x1_max])
    plt.ylim([x2_min, x2_max])
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    plt.title('f = f1 + f2')

    fig.suptitle('Cloud Closure: Bivariate Gaussian PDF fit ('
                 +  var_name1 + var_name2 + '; t=' + str(time) + ', z=' + str(z)+')'
                 , fontsize=20)

    plt.savefig(os.path.join(
        fullpath_out, 'CloudClosure_figures_norm/CC_bivariate_'
        + var_name1 + '_' + var_name2 + '_' + str(time) + '_z' + str(np.int(z)) + 'm.png'
    ))

    plt.close()
    return
#--------------------------
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