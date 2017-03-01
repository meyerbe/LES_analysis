import netCDF4 as nc
import argparse
import os, sys
import pylab as plt
import numpy as np
from sklearn import mixture, preprocessing
import itertools

from matplotlib.colors import LogNorm

from read_in_files import read_in_nml
from read_in_files import read_in_netcdf_fields

label_size = 12
# plt.rcParams['xtick.direction']='out'
# plt.rcParams['ytick.direction']='out'
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['axes.labelsize'] = label_size
plt.rcParams['legend.fontsize'] = label_size
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['lines.linewidth'] = 2

"""
=======================================
NOTES:
=======================================
Computes an n-component Gaussian mixture model, being multivariate in the sense of X=(x(z1), x(z2), ..., x(zn)), where x
is the same variable (e.g. w)
"""

#----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("path")
    parser.add_argument("casename")
    args = parser.parse_args()
    case_name = args.casename
    dx, nx, dt = read_in_nml(args.path, case_name)

    global fullpath_out
    fullpath_out = args.path
    print('')
    print('fullpath_out:', fullpath_out)
    print('_______________________')
    files = os.listdir(os.path.join(args.path, 'fields'))
    N = len(files)
    print('Found the following directories', files, N)

    global time
    time = np.zeros((1))
    for d in files:
        time = np.sort(np.append(time, np.int(d[0:-3])))

    print('_______________________')
    '''
    Parameters:
        krange:                 indices of vertical levels for which the PDF is fitted
        zrange = krange*dz      z-values for which the PDF is fitted
        var_list:               list of variables that are included in (multi-variate) PDF

        nvar:                   number of z-values included for PDF fit
            univar PDF: nvar = len(krange)

    Modules:
        Gaussian_mixture(data, ncomp, var_name, t):
                compute Gaussian mixture model (GMM) for a given number of components 'ncomp' at given time t
                and plot data with PDF components and contours
                    - data: trainingsdata
                    - ncomp: number of components of GMM
                    - t: time
    '''
    global krange, zrange
    # zrange = np.arange(5, 21, 5)
    krange = np.arange(5, 31, 5)
    zrange = krange * dx[2]
    print('zrange: ', krange, zrange)
    if case_name == 'DCBLSoares':
        var_list = ['u', 'w', 's']
    elif case_name == 'Bomex':
        var_list = ['w', 's', 'qt']
    else:
        var_list = ['w', 's', 'qt', 'thetali', 'temperature']
        # var_list = ['w']


    '''UNIVAR: only for one variable'''
    global nvar
    nvar = len(krange)
    zmax = len(krange)
    # files_ = files[0:3]
    files_ = files
    print(files_)
    for var in var_list:
        data_all = np.ndarray(shape=(0, nvar))
        # data_all = 999999.9*np.ones(shape=(0, nvar))
        for d in files_:
            t = np.int(d[0:-3])
            fullpath_in = os.path.join(args.path, 'fields', d)
            print(fullpath_in)
            # nc_file_name = 'EM_nonlocal_' + str(d)
            # create_statistics_file(fullpath_out, nc_file_name, ncomp, nvar, len(zrange))
            try:
                data_ = read_in_netcdf_fields(var, fullpath_in).reshape((nx[0] * nx[1], nx[2]))
                print(var + ' in list')
            except:
                print(var + ' NOT in list')
                continue
            data = np.ndarray(shape=(nx[0] * nx[1], nvar))
            # data = -9999.9 * np.ones(shape=(nx[0] * nx[1], nvar))
            for k in range(nvar):
                data[:, k] = data_[:, krange[k]]
            data_all = np.append(data_all, data, axis=0)
            # print('.....', t, data.shape, data_all.shape, nx[0]*nx[1], nvar)

        #     if np.isnan(data).any():
        #         print('huiuiui')
        #         sys.exit()
        #
        #     for i in range(nx[0]*nx[1]):
        #         for k in range(nvar):
        #             if data[i,k] < -1000.0 or data[i,k] > 100000:
        #                 print('ohoh')
        #                 sys.exit()
        # for l in range(len(files_)*nx[0]*nx[1]):
        #     for k in range(nvar):
        #         if data_all[l,k] > 100000:
        #             print('OHOH', l, l/(nx[0]*nx[1]), k)
        #             sys.exit()
        #     if np.isnan(data_all).any():
        #         print('HUIUII')
        #         sys.exit()


        # (a) Gaussian Mixture Model: ncomp = 2 (univar)
        print('Gaussian Mixture Model: ncomp = 2, var = '+var)
        ncomp = 2
        # means, covariance, weights = Gaussian_mixture_univar(data, ncomp, var, np.int(d[0:-3]), zrange[0] * dx[2])
        clf = Gaussian_mixture(data_all, ncomp, var)
        # print('shapes: ', clf.means_.shape, clf.covariances_.shape, clf.weights_.shape, len(zrange))
        mean_tot, covar_tot = covariance_estimate_from_multicomp_pdf(clf)
        corr, corr_tot = compute_correlation(clf.covariances_, covar_tot, ncomp)
        plot_covar_matrix(covar_tot, corr_tot, zrange, zrange, var, ncomp)

        print('')

        # (b) Gaussian Mixture Model: ncomp = 3 (univar)
        print('Gaussian Mixture Model: ncomp = 3')
        ncomp = 3
        # means, covariance, weights = Gaussian_mixture_univar(data, ncomp, var, np.int(d[0:-3]), zrange[0] * dx[2])
        clf = Gaussian_mixture(data_all, ncomp, var)

        # (c) Gaussian Mixture Model with flexible #components
        print('Gaussian Mixture Model: BIC')
        # nvar = len(var_list)*len(zrange)
        # zmax = 2
        # Gaussian_mixture_ncompselection(data_all[:,0:nvar],var, t)

        # (d) Bayesian Gaussian Mixture Model
        print('Bayesian Gaussian Mixture Model: dirichlet process')
        # Bayesian_Gaussian_Mixture(data_all, var, np.int(d[0:-3]), zrange)

        # (e) Covariance Estimator (empirical; maximum likelihood)
        # Covariance_empirical(data_all)

        # (f) Minimum Determinant Covariance

        # (g) Sparse Inverse Covariance

        # (h) Generalized Moment Method (statsmodel)

    '''BIVAR'''
    zmax = len(krange)
    nvar = 2*len(krange)


    return


#----------------------------------------------------------------------
#----------------------------------------------------------------------
def Gaussian_mixture(data, ncomp, var_name):
    # compute Gaussian mixture model for a given number of components
    import itertools
    # (1) Compute Gaussian Mixture Model
    clf = mixture.GaussianMixture(n_components=ncomp, covariance_type='full')
    clf.fit(data)
    Y_ = clf.predict(data)      # predict the labels for the data points

    # (2) plot data and PDF components
    plt.figure()
    color_iter = itertools.cycle(['navy', 'c', 'cornflowerblue', 'gold',
                                  'darkorange'])
    # plot data according to their cluster
    for i, (mean, covar, color) in enumerate(zip(
            clf.means_, clf.covariances_, color_iter)):
        plt.scatter(data[Y_ == i, 0], data[Y_ == i, 1], s=1., alpha=.3, color=color)
    # add PDF contours
    plot_contour_12(data, clf, ncomp, color_iter, var_name, zrange)
    plt.title('Gaussian Mixture Model: ncomp = ' + np.str(ncomp))
    plt.savefig(fullpath_out + 'PDF_nonlocal_alltimes_figures/EM' + np.str(ncomp)+'_PDF_univar_'
                + var_name + '_z' + np.str(np.int(zrange[0])) + 'm_' + np.str(np.int(zrange[1])) + 'm.png')
    plt.close()


    # mean_tot, covar_tot = covariance_estimate_from_multicomp_pdf(clf)
    # plot_covar_matrix(covar_tot, zrange_m, zrange_m, var_name, ncomp, t, z)

    return clf
#----------------------------------------------------------------------
def Gaussian_mixture_ncompselection(data,var_name, t):
    import itertools
    from scipy import linalg
    import matplotlib.mlab as mlab

    # (1) compute GMM and BIC for all number of components in given range
    lowest_bic = np.infty
    bic = []
    n_components_range = range(1, 10)
    cv_types = ['spherical', 'tied', 'diag', 'full']
    for cv_type in cv_types:
        for n_components in n_components_range:
            # Fit a Gaussian mixture with EM
            gmm = mixture.GaussianMixture(n_components=n_components,
                                          covariance_type=cv_type)
            gmm.fit(data)
            # compute BIC = Bayesian Information Criterion
            bic.append(gmm.bic(data))
            if bic[-1] < lowest_bic:
                lowest_bic = bic[-1]
                best_gmm = gmm
    bic = np.array(bic)

    # (2) select best GMM = GMM with lowest BIC
    clf = best_gmm
    n_opt = clf.weights_.size
    print('.....', n_opt, data.shape, clf.covariances_.shape, clf.means_.shape)

    # (3) Plotting
    color_iter = itertools.cycle(['navy', 'turquoise', 'cornflowerblue','darkorange'])
    bars = []
    plt.figure(figsize=(5,12))

    # (3a) Plot the BIC scores
    spl = plt.subplot(3, 1, 1)
    for i, (cv_type, color) in enumerate(zip(cv_types, color_iter)):
        xpos = np.array(n_components_range) + .2 * (i - 2)
        bars.append(plt.bar(xpos, bic[i * len(n_components_range):
                            (i + 1) * len(n_components_range)],
                            width=.2, color=color))
    plt.xticks(n_components_range)
    plt.ylim([bic.min() * 1.01 - .01 * bic.max(), bic.max()])
    plt.title('BIC score per model')
    xpos = np.mod(bic.argmin(), len(n_components_range)) + .65 + \
           .2 * np.floor(bic.argmin() / len(n_components_range))
    plt.text(xpos, bic.min() * 0.97 + .03 * bic.max(), '*', fontsize=14)
    spl.set_xlabel('Number of components')
    spl.legend([b[0] for b in bars], cv_types)

    # (3b) Plot the data in clusters
    spl = plt.subplot(3, 1, 2)
    Y_ = clf.predict(data)
    for i, (mean, cov, color) in enumerate(zip(clf.means_, clf.covariances_,color_iter)):
        if not np.any(Y_ == i):
            continue
        plt.scatter(data[Y_ == i, 0], data[Y_ == i, 1], s=1., color=color, alpha=0.3)

    # (3c) Plot the winner PDF
    spl = plt.subplot(3, 1, 3)
    for i, (mean, cov, color) in enumerate(zip(clf.means_, clf.covariances_,color_iter)):
        if not np.any(Y_ == i):
            continue
        plt.scatter(data[Y_ == i, 0], data[Y_ == i, 1], s=1.2, color=color, alpha=0.3)

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
    means = clf.means_
    sigma = np.sqrt(clf.covariances_)
    for i in range(n_opt):
        sxy = sigma[i,1,0]**2
        Z = mlab.bivariate_normal(X, Y, sigmax=sigma[i,0,0], sigmay=sigma[i,1,1], mux=means[i,0], muy=means[i,1], sigmaxy=sxy)
        ax = plt.contour(X,Y,Z)
    plt.xticks(())
    plt.yticks(())
    plt.title('Selected GMM: '+ np.str(n_opt) + ' components')
    plt.subplots_adjust(hspace=.35, bottom=.02)
    plt.savefig(fullpath_out + 'PDF_nonlocal_alltimes_figures/EM_PDF_univariate_gmm_' + var_name + '_' + str(t) + '.png')
    plt.close()
    return
#----------------------------------------------------------------------
def Bayesian_Gaussian_Mixture(data,var_name, t, zrange):
    # class sklearn.mixture.BayesianGaussianMixture(n_components=1, covariance_type='full', tol=0.001, reg_covar=1e-06,
    #               max_iter=100, n_init=1, init_params='kmeans', weight_concentration_prior_type='dirichlet_process',
    #               weight_concentration_prior=None, mean_precision_prior=None, mean_prior=None,
    #               degrees_of_freedom_prior=None, covariance_prior=None, random_state=None, warm_start=False,
    #               verbose=0, verbose_interval=10
    import itertools
    # Fit a Dirichlet process Gaussian mixture using five components

    dpgmm = mixture.BayesianGaussianMixture(weight_concentration_prior_type="dirichlet_process",n_components=5,
                                            init_params='random',       # setting to 'random' does not help
                                            max_iter=1000,              # setting from 100 to 1000 helped
                                            covariance_type='full').fit(data)
    Y_ = dpgmm.predict(data)
    ncomp = dpgmm.means_.shape[0]
    print('n comp: ', dpgmm.means_.shape[0])

    # Plot for first and second variable
    plt.figure()
    color_iter = itertools.cycle(['navy', 'c', 'cornflowerblue', 'gold',
                                  'darkorange'])
    for i, (mean, covar, color) in enumerate(zip(
            dpgmm.means_, dpgmm.covariances_, color_iter)):
        plt.scatter(data[Y_==i,0], data[Y_==i,1],s=1.,alpha=.5, color=color)

    plot_contour_12(data, dpgmm, ncomp, color_iter, var_name, zrange)
    plt.title('Bayesian Gaussian Mixture Model: ncomp = '+np.str(ncomp) + ' (time: '+np.str(t)+')')
    plt.savefig(fullpath_out + 'PDF_nonlocal_alltimes_figures/EM_PDF_uniariate_dpgmm_'
                + var_name + '_' + str(t) + '_z'+np.str(zrange[0])+'m_'+np.str(zrange[1])+'m.png')
    plt.close()
    return

#----------------------------------------------------------------------
#----------------------------------------------------------------------
def compute_correlation(covar, covar_tot, ncomp):
    # compute correlation matrix = covar_ij / (sigma_i*sigma_j), where covar_ij is computed from the EM model
    global nvar

    sigma = np.ndarray(shape=(ncomp,nvar))
    correlation = np.array(covar, copy=True)
    sigma_tot = np.ndarray(nvar)
    corr_tot = np.array(covar_tot, copy=True)

    print('')
    print('covar: ', covar.shape, covar_tot.shape, ncomp, nvar)
    print('sigma: ', sigma.shape, correlation.shape)
    print(sigma_tot.shape, corr_tot.shape)
    print('')

    for i in range(nvar):
        sigma[:, i] = np.sqrt(covar[:, i, i])
        sigma_tot[i] = np.sqrt(covar_tot[i,i])
    for i in range(nvar):
        for j in range(nvar):
            correlation[:,i,j] = covar[:,i,j] / (sigma[:,i]*sigma[:,j])
            corr_tot[i,j] = covar_tot[i,j] / (sigma_tot[i]*sigma_tot[j])
    return correlation, corr_tot
#----------------------------------------------------------------------
def plot_covar_matrix(covar, corr, x_, y_, var_name, ncomp):
    print('Calling plot covar: '+var_name)
    # plot the correlation and covariance matrix
    global zrange
    from matplotlib import cm

    plt.figure(1,figsize=(18,20))
    # X, Y = np.meshgrid(x_, y_)
    plt.subplot(2,2,1)
    plt.imshow(covar.T, cmap=cm.coolwarm, interpolation='nearest', norm=LogNorm())
    plt.colorbar(shrink=0.75)
    ax = plt.gca()
    labels = [item.get_text() for item in ax.get_xticklabels()]
    if len(labels) == zrange.shape[0]+2:
        labels[1:-1] = x_
        ax.set_xticklabels(labels)
        labels[1:-1] = y_
        ax.set_yticklabels(labels)
    elif zrange.shape[0] == 4:
        labels[1:len(labels):2] = x_
        ax.set_xticklabels(labels)
        labels[1:9:2] = y_
        ax.set_yticklabels(labels)
    plt.grid()
    plt.xlabel('level [m]')
    plt.ylabel('level [m]')
    plt.title('total covariance')

    plt.subplot(2, 2, 2)
    plt.imshow(covar, cmap=cm.coolwarm, interpolation='nearest')
    plt.colorbar(shrink=0.75)
    # labels = [item.get_text() for item in ax.get_xticklabels()]
    # print('....labels', len(labels), zrange)
    ax = plt.gca()
    if len(labels) == zrange.shape[0] + 2:
        labels[1:-1] = x_
        ax.set_xticklabels(labels)
        labels[1:-1] = y_
        ax.set_yticklabels(labels)
    elif zrange.shape[0] == 4:
        labels[1:len(labels):2] = x_
        ax.set_xticklabels(labels)
        labels[1:9:2] = y_
        ax.set_yticklabels(labels)
    plt.grid()
    plt.xlabel('level [m]')
    plt.ylabel('level [m]')
    # ax = plt.gca()
    # ax.set_xticklabels(zrange.tolist())
    # ax.set_yticklabels(zrange.tolist())
    plt.title('total covariance')

    plt.subplot(2, 2, 3)
    plt.imshow(corr, cmap=cm.coolwarm, interpolation='nearest', norm=LogNorm())
    plt.colorbar(shrink=0.75)
    # labels = [item.get_text() for item in ax.get_xticklabels()]
    # plt.xlabel('level [m]')
    # plt.ylabel('level [m]')
    plt.title('total correlation')
    ax = plt.gca()
    if len(labels) == zrange.shape[0] + 2:
        labels[1:-1] = x_
        ax.set_xticklabels(labels)
        labels[1:-1] = y_
        ax.set_yticklabels(labels)
    elif zrange.shape[0] == 4:
        labels[1:len(labels):2] = x_
        ax.set_xticklabels(labels)
        labels[1:9:2] = y_
        ax.set_yticklabels(labels)
    plt.grid()
    plt.xlabel('level [m]')
    plt.ylabel('level [m]')
    plt.subplot(2, 2, 4)
    plt.imshow(corr, cmap=cm.coolwarm, interpolation='nearest')
    plt.colorbar(shrink=0.75)
    # labels = [item.get_text() for item in ax.get_xticklabels()]
    ax = plt.gca()
    if len(labels) == zrange.shape[0] + 2:
        labels[1:-1] = x_
        ax.set_xticklabels(labels)
        labels[1:-1] = y_
        ax.set_yticklabels(labels)
    elif zrange.shape[0] == 4:
        labels[1:len(labels):2] = x_
        ax.set_xticklabels(labels)
        labels[1:9:2] = y_
        ax.set_yticklabels(labels)
    plt.grid()
    plt.xlabel('level [m]')
    plt.ylabel('level [m]')
    plt.title('total correlation')

    # plt.suptitle('Total Covariance & Correlation Matrices: ' + var_name + r'  (z$\in$'+ np.str(zrange) + ')', fontsize=30)
    plt.suptitle('Total Covariance & Correlation Matrices: ' + var_name, fontsize=30)
    plt.savefig(os.path.join(fullpath_out, 'PDF_nonlocal_alltimes_figures', 'EM' + np.str(ncomp) + '_PDF_univar_covariance_'
                + var_name + '.png'))

    plt.close()

    return
#----------------------------------------------------------------------
def plot_contour_12(data, clf, ncomp, colors, var_name, zrange):
    # Description: plot contourplot of computed EM PDF by computing the bivariate PDF from the clf.sigma and clf.means
    colors = ['navy', 'c', 'cornflowerblue', 'gold','darkorange']
    import matplotlib.mlab as mlab
    global nvar     # joint PDF computed for nvar variables
    nvar_ = 2       # only consider two variables here

    # (1) make samples
    n_sample = 300
    x_ = np.linspace(np.amin(data[:, 0]), np.amax(data[:, 0]), n_sample)
    y_ = np.linspace(np.amin(data[:, 1]), np.amax(data[:, 1]), n_sample)
    XX_ = np.ndarray(shape=(n_sample ** nvar_, nvar_))
    X_ = np.ndarray(shape=(n_sample, n_sample))
    Y_ = np.ndarray(shape=(n_sample, n_sample))
    delta_i = n_sample
    for i in range(n_sample):
        for j in range(n_sample):
            shift = i * delta_i + j
            XX_[shift, 0] = x_[i]
            XX_[shift, 1] = y_[j]
            X_[i, j] = x_[j]
            Y_[i, j] = y_[i]
    print('')

    # (2) compute PDF from means and covariance and plot
    means = clf.means_
    try:
        sigma = np.sqrt(clf.covariances_)
        for i in range(ncomp):
            sxy = sigma[i, 1, 0] ** 2
            Z_pred = mlab.bivariate_normal(X_, Y_, sigmax=sigma[i, 0, 0], sigmay=sigma[i, 1, 1], mux=means[i, 0],
                                           muy=means[i, 1], sigmaxy=sxy)
            # Z_pred = mlab.bivariate_normal(X, Y, sigmax=sigma[i, 0, 0], sigmay=sigma[i, 1, 1], mux=means[i, 0],
            #                                muy=means[i, 1],sigmaxy=sxy)
            ax = plt.contour(x_, y_, Z_pred, colors=colors[i], label='new ' + str(i))
    except:
        print('negative value in covariance')

    plt.xlabel(var_name + ' at level z=' + np.str(zrange[0]) + 'm')
    plt.ylabel(var_name + ' at level z=' + np.str(zrange[1]) + 'm')

    # # (3) compute scoring
    # # clf is a multivariate model for nvar=len(zrange) variables --> make matrix of input right dimension
    # n = 2
    # A = np.zeros(shape=(XX_.shape[0],1))
    # while (n < nvar):
    #     XX_ = np.append(XX_, A, axis=1)
    #     n += 1
    # # ZZ_ = clf.score_samples(XX_)
    # # ZZ = np.ndarray(shape=(n_sample, n_sample))
    # # for j in range(n_sample**nvar_):
    # #     j_shift = np.mod(j,delta_i)
    # #     i_shift = (j - j_shift) / delta_i
    # #     ZZ[i_shift, j_shift] = ZZ_[j]

    return


#----------------------------------------------------------------------
def Covariance_empirical(data):
    # class sklearn.covariance.EmpiricalCovariance(store_precision=True, assume_centered=False)
    from sklearn.covariance import EmpiricalCovariance
    global krange

    emp_cov = EmpiricalCovariance().fit(data)
    # print(emp_cov.get_params())
    covar = emp_cov.covariance_
    # print(data.shape, covar.shape, zrange.shape, zrange)
    # print(covar)

    plot_covar_matrix(covar, krange, krange)
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
if __name__ == "__main__":
    main()