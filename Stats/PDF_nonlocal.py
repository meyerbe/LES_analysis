import netCDF4 as nc
import argparse
import os
import pylab as plt
import numpy as np
from sklearn import mixture

from read_in_files import read_in_nml
from read_in_files import read_in_netcdf_fields

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
        nvar = len(zrange)
        data = np.ndarray(shape=(nx[0] * nx[1], nvar))
        for k in range(nvar):
            data[:, k] = data_[:, zrange[k]]
        means, covariance, weights = Gaussian_mixture_univar(data, var, np.int(d[0:-3]), zrange[0] * dx[2])

        # (b) Gaussian Mixture Model with flexible #components
        print('Gaussian Mixture Model: BIC')
        nvar = 2
        Gaussian_mixture_ncompselection(data[:,0:nvar],var, t)


        # (c) Bayesian Gaussian Mixture Model

        # (d) Covariance Estimater (empirical; maximum likelihood)

        # (e) Sparse Inverse Covariance

        # (f) Generalized Moment Method (statsmodel)

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
    print('.....', n_opt)

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
    Y_ = clf.predict(data)
    for i, (mean, cov, color) in enumerate(zip(clf.means_, clf.covariances_,
                                               color_iter)):

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
    plt.savefig(fullpath_out + 'figures_PDF_nonlocal/EM_PDF_bivariate_gmm_' + var_name + '_' + str(t) + '.png')
    plt.close()
    return
#----------------------------------------------------------------------
def Gaussian_mixture_univar(data, var_name, time, z):
    ncomp = 2
    clf = mixture.GaussianMixture(n_components=ncomp, covariance_type='full')
    clf.fit(data)
    # print('means=' + np.str(clf.means_.shape) + ', ' +np.str(clf.means_))
    # print('covar=' + np.str(clf.covariances_.shape))

    # # if var_name1 == 'qt' or var_name2 == 'qt':
    # #     # data_aux = np.ndarray(shape=((nx * ny), nvar))
    # #     # if var_name1 == 'qt':
    # #     #     data_aux[:, 0] = data[:, 0] * 1e2
    # #     # else:
    # #     #     data_aux[:, 0] = data[:, 0]
    # #     # if var_name2 == 'qt':
    # #     #     data_aux[:, 1] = data[:, 1] * 1e2
    # #     # else:
    # #     #     data_aux[:, 1] = data[:, 1]
    # #     # clf_aux = mixture.GaussianMixture(n_components=2, covariance_type='full')
    # #     # clf_aux.fit(data_aux)
    # #     # plot_PDF_samples_qt(data, data_aux, var_name1, var_name2, clf, clf_aux, time, z)
    # #     plot_PDF_samples_qt(data, var_name1, var_name2, clf, time, z)
    # #     print('!!!! qt: factor 100')
    # #     # return clf_aux.means_, clf_aux.covariances_
    # # else:
    # #     plot_PDF_samples(data, var_name1, var_name2, clf, time, z)
    # #     # return clf.means_, clf.covariances_
    #
    # # plot_PDF_samples(data, var_name1, var_name2, clf, time, z)

    return clf.means_, clf.covariances_, clf.weights_
    # return

#----------------------------------------------------------------------
if __name__ == "__main__":
    main()