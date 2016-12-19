import netCDF4 as nc
import argparse
import os
import json as  simplejson
import pylab as plt
from matplotlib.colors import LogNorm

import pickle

import numpy as np
import itertools
from scipy import linalg
from sklearn import mixture


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
    args = parser.parse_args()
    # ______________________
    case_name = 'Bomex'
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
    # for d in files:
    #     time[i] = d[0:-3]
    for d in files:
        time = np.sort(np.append(time, np.int(d[0:-3])))
    # ______________________
    # ______________________
    '''
    zrange:     z-values for which the PDF is fitted
    var_list:   list of variables that are included in (multi-variate) PDF
    '''
    # zrange = map(int,np.linspace(0, 24, 13))
    zrange = map(int,np.linspace(0,10,3))
    print('zrange', zrange)
    # var_list = ['w','s','qt']
    var_list = ['w','s']



    '''
    (1) uni-variate PDF for single variable
    '''
    data = np.ndarray(shape=((nx*ny),1))
    # ---
    ncomp = 2
    nvar = 1
    means_ = np.ndarray(shape=(len(zrange), ncomp, nvar))
    covariance_ = np.zeros(shape=(len(zrange), ncomp, nvar, nvar))
    weights_ = np.zeros(shape=(len(zrange), 2))
    out_dict = {}       # for pickle output
    # ---
    for d in files:
        nc_file_name = 'EM2_univar_' + str(d)
        pkl_file_name = 'EM2_univar_' + str(d[0:-3]) + '.pkl'
        create_statistics_file(fullpath_out, nc_file_name, ncomp, nvar, len(zrange))
        fullpath_in = os.path.join(args.path, 'fields', d)
        print('fullpath_in', fullpath_in)
        for var in var_list:
            # pkl_file_name_var = 'EM2_univar_' + var + '_' + str(d[0:-3]) + '.pkl'
            print('')
            print('...............varvarvar...............', var)
            data_ = read_in_netcdf_fields(var,fullpath_in).reshape((nx*ny),nz)
            # count = 0
            for i in range(len(zrange)):
                iz = zrange[i]
                # print(zrange, i, iz)
                data[:,0] = data_[:,iz]
                means, covariance, weights = Gaussian_mixture_univariate(data, var, np.int(d[0:-3]), iz)
                # ______________
                means_[i,:,0] = means[:,0]
                covariance_[i,0,0,0] = covariance[0]
                covariance_[i,1,0,0] = covariance[1]
                weights_[i,:] = weights[:]
                # count += 1
                print('weights', weights.shape)


            # print(os.path.join(fullpath_out, nc_file_name))
            dump_variable(os.path.join(fullpath_out, nc_file_name),'means', means_, var, ncomp, nvar, len(zrange))
            dump_variable(os.path.join(fullpath_out, nc_file_name), 'covariances', covariance_, var, ncomp, nvar, len(zrange))
            dump_variable(os.path.join(fullpath_out, nc_file_name), 'weights', weights_, var, ncomp, nvar, len(zrange))
            print(means_.shape, len(zrange))
            print('...............varvarvar...............', var)
            print('')

            # Python Dictionary = 'associative memories' or 'associative arrays'
            #   - indexed by keys
            #   - unordered set of 'key:value' pairs
            # out_dict = {}
            var_dict = {}
            var_dict['means'] = np.array(means_, dtype=np.double)
            var_dict['covars'] = np.array(covariance_, dtype=np.double)
            out_dict[var] = var_dict

            # dump_pickle(var_dict, fullpath_out, pkl_file_name_var)
            # test_pickle(fullpath_out, pkl_file_name_var)
        dump_pickle(out_dict, fullpath_out, pkl_file_name)
        test_pickle(fullpath_out,pkl_file_name)


    '''
    (2) bi-variate PDF for (s,qt,w)
    '''
    # ncomp = 2
    # nvar = 2
    # data = np.ndarray(shape=((nx * ny), nvar))
    # means_ = np.ndarray(shape=(len(zrange), ncomp, nvar))
    # covariance_ = np.zeros(shape=(len(zrange), ncomp, nvar, nvar))
    # # ---
    # for d in files:
    #     nc_file_name = 'EM2_bivar_' + str(d)
    #     create_statistics_file(fullpath_out, nc_file_name, ncomp, nvar, len(zrange))
    #
    #     fullpath_in = os.path.join(args.path, 'fields', d)
    #     print(fullpath_in)
    #     for var1 in var_list:
    #         data1_ = read_in_netcdf_fields(var1,fullpath_in).reshape((nx*ny,nz))
    #         for var2 in var_list:
    #             data2_ = read_in_netcdf_fields(var2, fullpath_in).reshape((nx*ny,nz))
    #             count = 0
    #             for i in zrange:
    #                 data[:,0] = data1_[:,i]
    #                 data[:,1] = data2_[:, i]
    #
    #                 means, covariance = Gaussian_mixture_bivariate(data, var1, var2, np.int(d[0:-3]), i*dz)
    #
    #                 means_[count, :, :] = means[:, :]
    #                 covariance_[count,:,:,:] = covariance[:,:,:]
    #                 count += 1
    #
    #             dump_variable(os.path.join(fullpath_out, nc_file_name),'means', means_, var1+var2, ncomp, nvar, len(zrange))
    #             dump_variable(os.path.join(fullpath_out, nc_file_name), 'covariances', covariance_, var1+var2, ncomp, nvar, len(zrange))

    # '''
    # (3) tri - variate PDF for (s, qt, w)
    # '''
    # ncomp = 2
    # nvar = 3
    # data = np.ndarray(shape=((nx * ny), nvar))
    # means_ = np.ndarray(shape=(len(zrange), ncomp, nvar))
    # covariance_ = np.zeros(shape=(len(zrange), ncomp, nvar, nvar))
    # print('means: ', means_.shape)
    # # ---
    # for d in files:
    #     nc_file_name = 'EM2_bivar_' + str(d)
    #     create_statistics_file(fullpath_out, nc_file_name, ncomp, nvar, len(zrange))
    #
    #     fullpath_in = os.path.join(args.path, 'fields', d)
    #     print(fullpath_in)
    #     for var1 in ['w']:
    #         data1_ = read_in_netcdf_fields(var1,fullpath_in).reshape((nx*ny,nz))
    #         for var2 in ['s']:
    #             data2_ = read_in_netcdf_fields(var2, fullpath_in).reshape((nx*ny,nz))
    #             for var3 in ['qt']:
    #                 data3_ = read_in_netcdf_fields(var3, fullpath_in).reshape((nx * ny, nz))
    #                 for i in zrange:
    #                     data[:,0] = data1_[:,i]
    #                     data[:,1] = data2_[:,i]
    #                     data[:,2] = data3_[:,i]
    #                     print('.........', i*dz)
    #                     means, covariance = Gaussian_mixture_trivariate(data, var1, var2, var3, np.int(d[0:-3]), i*dz)
    #


    return


#----------------------------------------------------------------------
def Gaussian_mixture_univariate(data, var_name, time, iz):
    for i in range(1):
        print('Gaussian mixture: '+ var_name + ', height: '+np.str(iz*dz))
        clf = mixture.GaussianMixture(n_components=2,covariance_type='full')

        clf.fit(data[:,i].reshape(nx*ny,1))

        n_sample = 100
        x_max = np.amax(data[:,i])
        x_min = np.amin(data[:,i])
        x = np.linspace(x_min,x_max,n_sample).reshape(n_sample,1)
        score = clf.score_samples(x)

        # print(var_name + ': means=' + np.str(clf.means_))
        # print(var_name + ': covar=' + np.str(clf.covariances_))
        # print(var_name + ': ', clf.means_.shape, clf.covariances_.shape)

        plt.subplot(3,1,1)
        plt.hist(data, bins=30)
        plt.title(var_name + ' (data), t='+str(time)+', z='+str(iz*dz))
        plt.ylabel('samples', fontsize=10)
        plt.subplot(3,1,2)
        plt.plot(x, np.exp(score))
        plt.plot([clf.means_[0], clf.means_[0]], [0, np.amax(np.exp(score))], 'k')
        plt.plot([clf.means_[1],clf.means_[1]],[0,np.amax(np.exp(score))],'k')
        plt.plot([clf.means_[0]-clf.covariances_[0,0], clf.means_[0]+clf.covariances_[0,0]], [np.amax(np.exp(score))/2, np.amax(np.exp(score))/2], 'k')
        plt.plot([clf.means_[1] - clf.covariances_[1, 0], clf.means_[1] + clf.covariances_[1, 0]],[np.amax(np.exp(score)) / 2, np.amax(np.exp(score)) / 2], 'k')
        plt.title('EM fit: likelihood',fontsize=10)
        plt.ylabel('likelihood')
        plt.subplot(3, 1, 3)
        plt.plot(x, score)
        if var_name == 'w':
            plt.xlabel('w [m/s]', fontsize=10)
        elif var_name == 's':
            plt.xlabel('s [J/K]', fontsize=10)
        else:
            plt.xlabel(var_name, fontsize=10)
        plt.ylabel('log likelihood',fontsize=10)
        plt.savefig('./figures_EM/EM2_PDF_univar_'+var_name+'_'+str(time)+'_z'+str(np.int(iz*dz))+'.png')
        plt.close()

    return clf.means_, clf.covariances_, clf.weights_


#----------------------------------------------------------------------
def Gaussian_mixture_bivariate(data, var_name1, var_name2, time, z):
    clf = mixture.GaussianMixture(n_components=2,covariance_type='full')
    clf.fit(data)
    print('means=' + np.str(clf.means_))
    print('covar=' + np.str(clf.covariances_))
    print(clf.means_.shape, clf.covariances_.shape)

    # Plotting
    n_sample = 100
    x1_max = np.amax(data[:,0])
    x1_min = np.amin(data[:,0])
    x2_max = np.amax(data[:,1])
    x2_min = np.amin(data[:,1])
    x = np.linspace(x1_min,x1_max,n_sample)
    y = np.linspace(x2_min,x2_max,n_sample)
    X, Y = np.meshgrid(x,y)
    XX = np.array([X.ravel(),Y.ravel()]).T
    Z = clf.score_samples(XX).reshape(X.shape)

    plt.figure(figsize=(8,16))
    plt.subplot(3,1,1)
    plt.scatter(data[:,0], data[:,1], s=5, alpha=0.5)
    plt.title(var_name1 + var_name2 + ' (data), t='+str(time)+', z='+str(z))
    # ax1 = plt.contour(X,Y,np.exp(Z),levels=np.linspace(10,20,11))
    ax1 = plt.contour(X, Y, Z, levels=np.linspace(10, 20, 11))
    plt.colorbar(ax1,shrink=0.8)
    plt.xlabel('w [m/s]')
    plt.ylabel('entropy s [K]')
    plt.subplot(3,1,2)
    # ax1 = plt.contour(X,Y,np.exp(Z),levels=np.linspace(0,1,11))
    # plt.plot([clf.means_[0,0]],[clf.means_[0,1]], 'o', markersize=10)
    # plt.plot([clf.means_[1, 0]], [clf.means_[1, 1]], 'o', markersize=10)
    ax1 = plt.contourf(X,Y,np.exp(Z))
    plt.colorbar(ax1,shrink=0.8)
    plt.title('EM fit: likelihood')
    plt.xlabel('w [m/s]')
    plt.ylabel('entropy s [K]')
    plt.subplot(3, 1, 3)
    plt.scatter(data[:, 0], data[:, 1], s=5, alpha=0.4)
    ax1 = plt.contour(X, Y, np.exp(Z), levels=np.linspace(0, 1, 11), linewidths=3)
    plt.plot([clf.means_[0, 0]], [clf.means_[0, 1]], 'o', markersize=10)
    plt.plot([clf.means_[1, 0]], [clf.means_[1, 1]], 'o', markersize=10)
    plt.colorbar(ax1, shrink=0.8)
    # plt.title(var_name)
    # plt.title(np.str(clf.means_))
    plt.xlabel('w [m/s]')
    plt.ylabel('entropy s [K]')
    plt.savefig('./figures_EM/EM_PDF_bivariate_'+var_name1+'_'+var_name2+'_'+str(time)+'_z'+str(np.int(z))+'.png')
    plt.close()

    return clf.means_, clf.covariances_




#----------------------------------------------------------------------
def Gaussian_mixture_trivariate(data, var_name1, var_name2, var_name3, time, z):
    clf = mixture.GaussianMixture(n_components=2,covariance_type='full')
    clf.fit(data)
    print('trivar means=' + np.str(clf.means_))
    print('trivar covar=' + np.str(clf.covariances_))
    print(clf.means_.shape, clf.covariances_.shape)

    # Plotting
    n_sample = 100
    x1_max = np.amax(data[:,0])
    x1_min = np.amin(data[:,0])
    x2_max = np.amax(data[:,1])
    x2_min = np.amin(data[:,1])
    x3_max = np.amax(data[:, 0])
    x3_min = np.amin(data[:, 1])
    x = np.linspace(x1_min,x1_max,n_sample)
    y = np.linspace(x2_min,x2_max,n_sample)
    z = np.linspace(x3_min, x3_max, n_sample)
    X, Y, Z = np.meshgrid(x,y,z)
    XXX = np.array([X.ravel(),Y.ravel(),Z.ravel()]).T
    S = clf.score_samples(XXX).reshape(X.shape)
    print('!!!', S.shape, X.shape, Y.shape, Z.shape)

    plt.figure(figsize=(12,16))
    plt.subplot(3,3,1)
    plt.scatter(data[:,0], data[:,1], s=5, alpha=0.5)
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    plt.title(var_name1 + var_name2+ ' (data), t='+ str(time)+ ', z=')
    # plt.title(var_name1 + var_name2 + ' (data), t=' + str(time) + ', z=' + str(z))
    plt.subplot(3, 3, 2)
    plt.scatter(data[:, 1], data[:, 2], s=5, alpha=0.5)
    plt.xlabel(var_name2)
    plt.ylabel(var_name3)
    plt.title(var_name2 + var_name3 + ' (data)')
    plt.subplot(3, 3, 3)
    plt.scatter(data[:, 0], data[:, 2], s=5, alpha=0.5)
    plt.xlabel(var_name1)
    plt.ylabel(var_name3)
    plt.title(var_name1 + var_name3 + ' (data), ')
    # # ax1 = plt.contour(X,Y,np.exp(Z),levels=np.linspace(10,20,11))
    # ax1 = plt.contour(X, Y, Z, levels=np.linspace(10, 20, 11))
    # plt.colorbar(ax1,shrink=0.8)
    # plt.xlabel('w [m/s]')
    # plt.ylabel('entropy s [K]')
    # plt.subplot(3,1,2)
    plt.subplot(3, 3, 4)
    # plt.scatter(data[:, 0], data[:, 1], s=5, alpha=0.5)
    ax1 = plt.contourf(X[:,:,0],Y[:,:,0],np.exp(S[:,:,0]))#,levels=np.linspace(0,100,11), linewidth=3)
    plt.colorbar(ax1)
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    plt.subplot(3, 3, 5)
    # plt.scatter(data[:, 1], data[:, 2], s=5, alpha=0.5)
    ax1 = plt.contourf(X[:, 1, :], Z[:, 1, :], np.exp(S[:, 1, :]))  # ,levels=np.linspace(0,100,11), linewidth=3)
    plt.colorbar(ax1)
    plt.xlabel(var_name2)
    plt.ylabel(var_name3)
    plt.subplot(3, 3, 6)
    plt.scatter(data[:, 0], data[:, 2], s=5, alpha=0.5)
    plt.xlabel(var_name1)
    plt.ylabel(var_name3)
    # # ax1 = plt.contour(X,Y,np.exp(Z),levels=np.linspace(0,1,11))
    # # plt.plot([clf.means_[0,0]],[clf.means_[0,1]], 'o', markersize=10)
    # # plt.plot([clf.means_[1, 0]], [clf.means_[1, 1]], 'o', markersize=10)
    # ax1 = plt.contourf(X,Y,np.exp(Z))
    # plt.colorbar(ax1,shrink=0.8)
    # plt.title('EM fit: likelihood')
    # plt.xlabel('w [m/s]')
    # plt.ylabel('entropy s [K]')
    # plt.subplot(3, 1, 3)
    plt.subplot(3, 3, 7)
    plt.scatter(data[:, 0], data[:, 1], s=5, alpha=0.5)
    # ax1 = plt.contour(X, Y, np.exp(Z), levels=np.linspace(0, 1, 11), linewidths=3)
    # plt.colorbar(ax1, shrink=0.8)
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    plt.subplot(3, 3, 8)
    plt.scatter(data[:, 1], data[:, 2], s=5, alpha=0.5)
    # ax1 = plt.contour(X, Y, np.exp(Z), levels=np.linspace(0, 1, 11), linewidths=3)
    # plt.colorbar(ax1, shrink=0.8)
    plt.xlabel(var_name2)
    plt.ylabel(var_name3)
    plt.subplot(3, 3, 9)
    plt.scatter(data[:, 0], data[:, 2], s=5, alpha=0.5)
    # ax1 = plt.contour(X, Y, np.exp(Z), levels=np.linspace(0, 1, 11), linewidths=3)
    # plt.colorbar(ax1, shrink=0.8)
    plt.xlabel(var_name1)
    plt.ylabel(var_name3)
    # plt.plot([clf.means_[0, 0]], [clf.means_[0, 1]], 'o', markersize=10)
    # plt.plot([clf.means_[1, 0]], [clf.means_[1, 1]], 'o', markersize=10)
    # plt.xlabel('w [m/s]')
    # plt.ylabel('entropy s [K]')
    plt.savefig('./figures_EM/EM_PDF_trivariate_'+var_name1+'_'+var_name2+'_'+var_name3+'_'+str(time)+'_z'+'.png')
    plt.close()
    # plt.show()

    return clf.means_, clf.covariances_


# ____________________
def dump_pickle(data,out_path,file_name):
    data_ = (1.4,42)
    # output = open(os.path.join(out_path,'data.pkl'), 'w')
    output = open(os.path.join(out_path, file_name), 'w')
    pickle.dump(data, output)
    output.close()
    return

def test_pickle(in_path,file_name):
    print ''
    print '------- test pickle ------'
    fullpath_in = os.path.join(in_path,file_name)
    f = open(fullpath_in)
    data = pickle.load(f)
    print(data)
    print ''
    var = data['w']
    print(var)
    print
    ''
    means_ = var['means']
    print(means_)
    print '-------------------------'
    print ''
    return

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
    means_grp.createDimension('nz', nz_)
    means_grp.createDimension('ncomp', ncomp)
    means_grp.createDimension('nvar', nvar)
    cov_grp.createDimension('nz', nz_)
    cov_grp.createDimension('ncomp', ncomp)
    cov_grp.createDimension('nvar', nvar)
    weights_grp.createDimension('nz', nz_)
    weights_grp.createDimension('EM2', 2)
    rootgrp.close()
    print('create file end')
    return

def dump_variable(path, group_name, data_, var_name, ncomp, nvar, nz_):
    print('--------')
    print('dump variable', path, group_name, var_name, data_.shape, ncomp, nvar)
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
        print('dump weights')
        add_weights(path, var_name, ncomp, nvar)
        data = np.empty((nz_, ncomp), dtype=np.double, order='c')
        for i in range(nz_):
            for j in range(ncomp):
                data[i, j] = data_[i, j]
        write_weights(path, group_name, data, var_name)

    # write_field(path, group_name, data, var_name)
    print('--------')
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
    print('add covariance: ', var_name, path)
    # rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
    rootgrp = nc.Dataset(path, 'r+')
    group = rootgrp.groups['covariances']
    var = group.createVariable(var_name, 'f8', ('nz', 'ncomp', 'nvar', 'nvar'))
    rootgrp.close()
    return

def add_weights(path, var_name, ncomp, nvar):
    print('add weights: ', var_name, path)
    # rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
    rootgrp = nc.Dataset(path, 'r+')
    group = rootgrp.groups['weights']
    var = group.createVariable(var_name, 'f8', ('nz', 'EM2'))
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