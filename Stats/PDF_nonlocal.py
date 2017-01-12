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
        Gaussian_mixture_ncompselection(data)


        # (c) Bayesian Gaussian Mixture Model

        # (d) Covariance Estimater (empirical; maximum likelihood)

        # (e) Sparse Inverse Covariance

        # (f) Generalized Moment Method (statsmodel)

    return
#----------------------------------------------------------------------
def Gaussian_mixture_ncompselection(data):
    lowest_bic = np.infty
    bic = []
    n_components_range = range(1, 7)
    cv_types = ['spherical', 'tied', 'diag', 'full']
    for cv_type in cv_types:
        for n_components in n_components_range:
            # Fit a Gaussian mixture with EM
            gmm = mixture.GaussianMixture(n_components=n_components,
                                          covariance_type=cv_type)
            gmm.fit(data)
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