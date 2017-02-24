import netCDF4 as nc
import argparse
import os
import pylab as plt
import numpy as np
from sklearn import mixture

from matplotlib.colors import LogNorm

from read_in_files import read_in_nml
from read_in_files import read_in_netcdf_fields

label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size

"""
=======================================
NOTES:
=======================================
- Bayesian Gaussian Mixture Model instable: does not always produce the same output (i.e. the same PDF) for same data
    --> because of problem in sqrt(covariance)??

- ??? Is vertical correlation really Gaussian? Look at EM_PDF_bivar:
    means show more of an exponential / gamma distribution shape;
    no real pattern in covariances
- ??? EM Univar: flat tail of histograms --> might be better represented by a combination of Gaussian and GAMMA DISTRIBUTION
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
    print('fullpath_out:', fullpath_out)
    print('_______________________')
    files = os.listdir(os.path.join(args.path, 'fields'))
    N = len(files)
    print('Found the following directories', files, N, files[0])

    global time
    time = np.zeros((1))
    for d in files:
        time = np.sort(np.append(time, np.int(d[0:-3])))

    print('_______________________')
    '''
        zrange:     z-values for which the PDF is fitted
        var_list:   list of variables that are included in (multi-variate) PDF
    '''
    global zrange, zrange_m
    # zrange = np.arange(5, 21, 5)
    zrange = np.arange(5, 31, 5)
    zrange_m = zrange*dx[2]
    print('zrange: ', zrange, zrange_m)
    if case_name == 'DCBLSoares':
        var_list = ['u', 'w', 's']
    else:
        var_list = ['w', 's', 'qt']
        # var_list = ['w']


    '''UNIVAR'''
    global nvar
    nvar = len(zrange)
    zmax = len(zrange)
    for d in files:
        t = np.int(d[0:-3])
        fullpath_in = os.path.join(args.path, 'fields', d)
        print(fullpath_in)
        # nc_file_name = 'EM_nonlocal_' + str(d)
        # create_statistics_file(fullpath_out, nc_file_name, ncomp, nvar, len(zrange))

        for var in var_list:
        # var = var_list[0]
            data_ = read_in_netcdf_fields(var, fullpath_in).reshape((nx[0] * nx[1], nx[2]))
            data = np.ndarray(shape=(nx[0] * nx[1], zmax))
            for k in range(zmax):
                data[:, k] = data_[:, zrange[k]]

            # (a) Gaussian Mixture Model: ncomp = 2 (univar)
            print('Gaussian Mixture Model: ncomp = 2')
            ncomp = 2
            means, covariance, weights = Gaussian_mixture_univar(data, ncomp, var, np.int(d[0:-3]), zrange[0] * dx[2])
            print('shapes: ', means.shape, covariance.shape, weights.shape, len(zrange))
            correlation = compute_correlation(covariance, ncomp)
            print('')

            # (b) Gaussian Mixture Model: ncomp = 3 (univar)
            print('Gaussian Mixture Model: ncomp = 3')
            ncomp = 3
            means, covariance, weights = Gaussian_mixture_univar(data, ncomp, var, np.int(d[0:-3]), zrange[0] * dx[2])

        # (c) Gaussian Mixture Model with flexible #components
        print('Gaussian Mixture Model: BIC')
        # nvar = len(var_list)*len(zrange)
        zmax = 2
        # Gaussian_mixture_ncompselection(data[:,0:nvar],var, t)

        # (d) Bayesian Gaussian Mixture Model
        print('Bayesian Gaussian Mixture Model: dirichlet process')
        # Bayesian_Gaussian_Mixture(data, var, np.int(d[0:-3]), zrange)

        # (e) Covariance Estimator (empirical; maximum likelihood)
        # Covariance_empirical(data)

        # (f) Minimum Determinant Covariance

        # (g) Sparse Inverse Covariance

        # (h) Generalized Moment Method (statsmodel)

    return


#----------------------------------------------------------------------
def compute_correlation(covariance, ncomp):
    global nvar
    sigma = np.ndarray(shape=(ncomp,nvar))
    correlation = np.array(covariance, copy=True)
    print('nvar, covariance:', nvar, covariance.shape, sigma.shape, correlation.shape)

    for i in range(nvar):
        sigma[:,i] = np.sqrt(covariance[:,i,i])

    for i in range(nvar):
        for j in range(nvar):
            correlation[:,i,j] = covariance[:,i,j] / (sigma[:,i]*sigma[:,j])

    return correlation
#----------------------------------------------------------------------
def plot_covar_matrix(data, x_, y_, var_name, ncomp, t, z):
    from matplotlib import cm
    fig = plt.figure(figsize=(10,5))
    fig, ax = plt.subplots()
    X, Y = np.meshgrid(x_, y_)
    plt.subplot(1,2,1)
    plt.imshow(data, cmap=cm.coolwarm, interpolation='nearest', norm=LogNorm())
    plt.colorbar(shrink=0.5)
    labels = [item.get_text() for item in ax.get_xticklabels()]
    # print('....labels', len(labels), zrange)
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

    plt.subplot(1, 2, 2)
    plt.imshow(data, cmap=cm.coolwarm, interpolation='nearest')
    plt.colorbar(shrink=0.5)
    labels = [item.get_text() for item in ax.get_xticklabels()]
    # print('....labels', len(labels), zrange)
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
    #  plt.xlabel('delta z')
    # plt.ylabel('delta z')

    plt.suptitle('Total Covariance Matrix: ' + np.str(data.shape))
    plt.savefig(fullpath_out + 'figures_PDF_nonlocal/EM' + np.str(ncomp) + '_PDF_univar_covariance_'
                + var_name + '_' + str(t) + '.png')

    plt.close()

    return
#----------------------------------------------------------------------
def plot_contour_12(data, clf, ncomp, colors, var_name, zrange):
    colors = ['navy', 'c', 'cornflowerblue', 'gold','darkorange']
    import itertools
    import matplotlib.mlab as mlab

    # print('plotting: ', data.shape)
    n_sample = 300
    x1_max = np.amax(data[:, 0])
    x1_min = np.amin(data[:, 0])
    x2_max = np.amax(data[:, 1])
    x2_min = np.amin(data[:, 1])
    x = np.linspace(x1_min, x1_max, n_sample)
    y = np.linspace(x2_min, x2_max, n_sample)
    X, Y = np.meshgrid(x, y)
    XX = np.array([X.ravel(), Y.ravel()]).T

    nvar = 2
    A = np.zeros(shape=X.shape)
    while (nvar < data.shape[1]):
        XX = np.append(XX, A.ravel().reshape(X.size, 1), axis=1)
        nvar += 1
    Z = clf.score_samples(XX).reshape(X.shape)
    means = clf.means_
    sigma = np.sqrt(clf.covariances_)
    for i in range(ncomp):
        # pass
        sxy = sigma[i, 1, 0] ** 2
        Z = mlab.bivariate_normal(X, Y, sigmax=sigma[i, 0, 0], sigmay=sigma[i, 1, 1], mux=means[i, 0], muy=means[i, 1],
                                  sigmaxy=sxy)
        ax = plt.contour(X, Y, Z, colors = colors[i])
    plt.xlabel(var_name + ' at level z=' + np.str(zrange[0]) + 'm')
    plt.ylabel(var_name + ' at level z=' + np.str(zrange[1]) + 'm')
    return

#----------------------------------------------------------------------
def Gaussian_mixture_univar(data, ncomp, var_name, t, z):
    import itertools
    # (1) Compute Gaussian Mixture Model
    clf = mixture.GaussianMixture(n_components=ncomp, covariance_type='full')
    clf.fit(data)
    Y_ = clf.predict(data)

    # (2) plot data and PDF components
    plt.figure()
    color_iter = itertools.cycle(['navy', 'c', 'cornflowerblue', 'gold',
                                  'darkorange'])
    for i, (mean, covar, color) in enumerate(zip(
            clf.means_, clf.covariances_, color_iter)):
        plt.scatter(data[Y_ == i, 0], data[Y_ == i, 1], s=1., alpha=.5, color=color)

    plot_contour_12(data, clf, ncomp, color_iter, var_name, zrange)
    plt.title('Gaussian Mixture Model: ncomp = ' + np.str(ncomp) + ' (time: ' + np.str(t) + ')')
    plt.savefig(fullpath_out + 'figures_PDF_nonlocal/EM' + np.str(ncomp)+'_PDF_univar_'
                + var_name + '_' + str(t) + '_z' + np.str(zrange[0]) + 'm_' + np.str(zrange[1]) + 'm.png')



    mean_tot, covar_tot = covariance_estimate_from_multicomp_pdf(clf)
    plot_covar_matrix(covar_tot, zrange_m, zrange_m, var_name, ncomp, t, z)





    # if var_name1 == 'qt' or var_name2 == 'qt':
    #     # data_aux = np.ndarray(shape=((nx * ny), nvar))
    #     # if var_name1 == 'qt':
    #     #     data_aux[:, 0] = data[:, 0] * 1e2
    #     # else:
    #     #     data_aux[:, 0] = data[:, 0]
    #     # if var_name2 == 'qt':
    #     #     data_aux[:, 1] = data[:, 1] * 1e2
    #     # else:
    #     #     data_aux[:, 1] = data[:, 1]
    #     # clf_aux = mixture.GaussianMixture(n_components=2, covariance_type='full')
    #     # clf_aux.fit(data_aux)
    #     # plot_PDF_samples_qt(data, data_aux, var_name1, var_name2, clf, clf_aux, time, z)
    #     plot_PDF_samples_qt(data, var_name1, var_name2, clf, time, z)
    #     print('!!!! qt: factor 100')
    #     # return clf_aux.means_, clf_aux.covariances_
    # else:
    #     plot_PDF_samples(data, var_name1, var_name2, clf, time, z)
    #     # return clf.means_, clf.covariances_
    # plot_PDF_samples(data, var_name1, var_name2, clf, time, z)

    return clf.means_, clf.covariances_, clf.weights_
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
    plt.savefig(fullpath_out + 'figures_PDF_nonlocal/EM_PDF_univariate_gmm_' + var_name + '_' + str(t) + '.png')
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
    plt.savefig(fullpath_out + 'figures_PDF_nonlocal/EM_PDF_uniariate_dpgmm_'
                + var_name + '_' + str(t) + '_z'+np.str(zrange[0])+'m_'+np.str(zrange[1])+'m.png')
    plt.close()

    return


#----------------------------------------------------------------------
def Covariance_empirical(data):
    # class sklearn.covariance.EmpiricalCovariance(store_precision=True, assume_centered=False)
    from sklearn.covariance import EmpiricalCovariance
    global zrange

    emp_cov = EmpiricalCovariance().fit(data)
    # print(emp_cov.get_params())
    covar = emp_cov.covariance_
    # print(data.shape, covar.shape, zrange.shape, zrange)
    # print(covar)

    plot_covar_matrix(covar, zrange, zrange)
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