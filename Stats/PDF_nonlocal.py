import netCDF4 as nc
import argparse
import os
import pylab as plt
import numpy as np
from sklearn import mixture

from read_in_files import read_in_nml
from read_in_files import read_in_netcdf_fields


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
    global zrange
    zrange = np.arange(5, 21, 5)
    print('zrange: ', zrange, zrange*dx[2])
    if case_name == 'DCBLSoares':
        var_list = ['u', 'w', 's']
    else:
        var_list = ['w', 's', 'qt']
        var_list = ['w']

    '''UNIVAR'''
    for d in files:
        t = np.int(d[0:-3])
        fullpath_in = os.path.join(args.path, 'fields', d)
        print(fullpath_in)
        # nc_file_name = 'EM_nonlocal_' + str(d)
        # create_statistics_file(fullpath_out, nc_file_name, ncomp, nvar, len(zrange))

        var = var_list[0]
        data_ = read_in_netcdf_fields(var, fullpath_in).reshape((nx[0] * nx[1], nx[2]))

        # (a) Gaussian Mixture Model: ncomp = 2 (univar)
        print('Gaussian Mixture Model: ncomp = 2')
        ncomp = 2
        nvar = len(zrange)
        data = np.ndarray(shape=(nx[0] * nx[1], nvar))
        for k in range(nvar):
            data[:, k] = data_[:, zrange[k]]
        means, covariance, weights = Gaussian_mixture_univar(data, ncomp, var, np.int(d[0:-3]), zrange[0] * dx[2])

        # # (b) Gaussian Mixture Model: ncomp = 3 (univar)
        print('Gaussian Mixture Model: ncomp = 3')
        ncomp = 3
        nvar = len(zrange)
        data = np.ndarray(shape=(nx[0] * nx[1], nvar))
        for k in range(nvar):
            data[:, k] = data_[:, zrange[k]]
        means, covariance, weights = Gaussian_mixture_univar(data, ncomp, var, np.int(d[0:-3]), zrange[0] * dx[2])

        # (c) Gaussian Mixture Model with flexible #components
        print('Gaussian Mixture Model: BIC')
        # nvar = len(var_list)*len(zrange)
        nvar = 2
        # Gaussian_mixture_ncompselection(data[:,0:nvar],var, t)

        # (d) Bayesian Gaussian Mixture Model
        print('Bayesian Gaussian Mixture Model: dirichlet process')
        Bayesian_Gaussian_Mixture(data, var, np.int(d[0:-3]), zrange)

        # (e) Covariance Estimator (empirical; maximum likelihood)
        Covariance_empirical(data)

        # (f) Minimum Determinant Covariance

        # (g) Sparse Inverse Covariance

        # (h) Generalized Moment Method (statsmodel)

    return


#----------------------------------------------------------------------
# class sklearn.covariance.EmpiricalCovariance(store_precision=True, assume_centered=False)
def Covariance_empirical(data):
    from sklearn.covariance import EmpiricalCovariance
    global zrange

    emp_cov = EmpiricalCovariance().fit(data)
    # print(emp_cov.get_params())
    covar = emp_cov.covariance_
    print(data.shape, covar.shape, zrange.shape, zrange)
    print(covar)

    plot_covar(covar, zrange, zrange)
    return

#----------------------------------------------------------------------
def plot_covar(data, x_, y_):
    from mpl_toolkits.mplot3d import Axes3D
    # import matplotlib.pyplot as plt
    from matplotlib import cm
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    X, Y = np.meshgrid(x_, y_)
    # R = np.sqrt(X ** 2 + Y ** 2)
    # Z = np.sin(R)
    # surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
    #                        linewidth=0, antialiased=False)
    surf = ax.plot_surface(X,Y,data, cmap=cm.coolwarm)
    # plt.show()
    plt.close()

    return

#----------------------------------------------------------------------
def plot_contour_12(data, clf, ncomp, colors, var_name, zrange):
    colors = ['navy', 'c', 'cornflowerblue', 'gold','darkorange']
    import itertools
    import matplotlib.mlab as mlab

    print('plotting: ', data.shape)
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
# class sklearn.mixture.BayesianGaussianMixture(n_components=1, covariance_type='full', tol=0.001, reg_covar=1e-06,
#               max_iter=100, n_init=1, init_params='kmeans', weight_concentration_prior_type='dirichlet_process',
#               weight_concentration_prior=None, mean_precision_prior=None, mean_prior=None,
#               degrees_of_freedom_prior=None, covariance_prior=None, random_state=None, warm_start=False,
#               verbose=0, verbose_interval=10
def Bayesian_Gaussian_Mixture(data,var_name, t, zrange):
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
def Gaussian_mixture_univar(data, ncomp, var_name, t, z):
    import itertools
    clf = mixture.GaussianMixture(n_components=ncomp, covariance_type='full')
    clf.fit(data)
    Y_ = clf.predict(data)

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
if __name__ == "__main__":
    main()