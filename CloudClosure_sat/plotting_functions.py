import os
import pylab as plt
import numpy as np
import matplotlib.mlab as mlab
from sklearn.preprocessing import StandardScaler
import time


label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['xtick.direction']='out'
plt.rcParams['ytick.direction']='out'
plt.rcParams['legend.fontsize'] = 10
# plt.rcParams['title.fontsize'] = 15


def plot_PDF(data, data_norm, var_name1, var_name2, nvar, clf, scaler, ncomp, error, z, path):
    print('Plot PDF')

    xmin = np.amin(data_norm[:,0])
    xmax = np.amax(data_norm[:,0])
    ymin = np.amin(data_norm[:,1])
    ymax = np.amax(data_norm[:,1])

    n_sample = np.int(1e2)
    x_ = np.linspace(xmin, xmax, n_sample)
    y_ = np.linspace(ymin, ymax, n_sample)
    XX_ = np.ndarray(shape=(n_sample**nvar,nvar))

    # (1) print PDF computed by GMM
    delta_i = n_sample
    for i in range(n_sample):
        for j in range(n_sample):
            shift = i*delta_i + j
            XX_[shift, 0] = x_[i]
            XX_[shift, 1] = y_[j]
    ZZ_ = clf.score_samples(XX_)
    ZZ = np.ndarray(shape=(n_sample,n_sample))
    for j in range(n_sample**nvar):
        jshift = np.mod(j,delta_i)
        ishift = (j - jshift)/delta_i
        ZZ[ishift,jshift] = ZZ_[j]

    # XXii = np.zeros(shape=XX_.shape)
    # print('XX before: ', ZZ.shape, XX_.shape)
    # print(XX_[0:3, 0])
    # print(XX_[0:3, 1])
    # print('XX after inverse')
    # XXi = scaler.inverse_transform(XX_)
    # print(XXi[0:3,0])
    # print(XXi[0:3, 1])
    # XXii[:,0] = XX_[:,0]/scaler.scale_[0]
    # XXii[:, 1] = XX_[:, 1]/ scaler.scale_[1]
    # print(scaler.scale_)
    # print('XX after rescale')
    # print(XXii[0:3, 0])
    # print(XXii[0:3, 1])

    # (2) compute normal distribution from parameters given from GMM
    # X_ = np.ndarray()
    # Z_pred =



    plt.figure(figsize=(12,8))
    plt.subplot(1,3,1)
    plt.scatter(data[:,0], data[:,1], s=5, alpha=0.2)
    plt.title('data')
    labeling(var_name1,var_name2,1,1)
    plt.xlim(np.amin(data[:,0]), np.amax(data[:,0]))
    plt.ylim(np.amin(data[:,1]), np.amax(data[:,1]))

    plt.subplot(1, 3, 2)
    ax = plt.contour(x_, y_, np.exp(ZZ).T)
    plt.scatter(data_norm[:, 0], data_norm[:, 1], s=5, alpha=0.2)
    plt.colorbar(ax)
    plt.title('data')
    labeling(var_name1, var_name2, scaler.scale_[0], scaler.scale_[1])
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

    plt.subplot(1,3,3)
    ax = plt.contourf(x_,y_,np.exp(ZZ).T)
    plt.colorbar(ax)
    # plt.scatter(data_norm[:, 0], data_norm[:, 1], s=5, alpha=0.2)
    plt.title('PDF(thl, qt)')
    labeling(var_name1, var_name2, scaler.scale_[0], scaler.scale_[1])
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

    plt.suptitle('GMM, ncomp='+str(ncomp)+', error: '+str(error)+'    (z='+str(z)+')')

    plt.savefig(os.path.join(path,'CloudClosure_figures','PDF_figures_'+str(z)+'m'+'_ncomp'+str(ncomp)+'.png'))
    # plt.show()
    plt.close()
    return


def plot_samples(type, data, data_ql, samples, samples_ql, var_name1, var_name2, scaler, ncomp, z, path):
    print('Plot samples')
    # type:         normalised data or original data
    # data:         data from LES output
    # data_ql:      ql data from LES output
    # samples:      sample data drawn from GMM PDF
    # samples_ql:   ql data calculated from sample data
    # var_name1:    thl or s
    # var_name :    qt
    # scaler:       scaling factors for normalisation
    # ncomp:        # of components for GMM PDF
    # path:         output path for saving figure

    # print('shapes: ', data.shape, data_ql.shape, samples.shape, samples_ql.shape)
    # print('min/max: ', np.amin(data_ql), np.amax(data_ql), np.amin(samples_ql), np.amax(samples_ql))

    # ql_min = np.min([np.amin(data_ql), np.amin(samples_ql)])
    ql_min = 0.0
    ql_max = np.max([np.amax(data_ql), np.amax(samples_ql)])

    scale_thl = scaler.scale_[0]
    scale_qt = scaler.scale_[1]
    xmin = np.amin(data[:,0])
    xmax = np.amax(data[:,0])
    ymin = np.amin(data[:,1])
    ymax = np.amax(data[:,1])
    xmin = np.min([np.amin(data[:, 0]),np.amin(samples[:, 0])])
    xmax = np.max([np.amax(data[:, 0]),np.amax(samples[:, 0])])
    ymin = np.min([np.amin(data[:, 1]),np.amin(samples[:, 1])])
    ymax = np.max([np.amax(data[:, 1]),np.amax(samples[:, 1])])
    cm = plt.cm.get_cmap('RdYlBu')
    cm = plt.cm.get_cmap('viridis')


    plt.figure(figsize=(8,9))
    plt.subplot(2,2,1)
    plt.scatter(data[:,0], data[:,1], s=5, alpha=0.2)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    if type == 'norm':
        plt.title('data norm')
    else:
        plt.title('data')
    labeling(var_name1, var_name2, scale_thl, scale_qt)

    plt.subplot(2,2,2)
    plt.scatter(samples[:,0], samples[:,1], s=5, alpha=0.2)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.title('samples (n='+str(np.size(samples_ql))+')')
    labeling(var_name1, var_name2, 0, 0)

    plt.subplot(2, 2, 3)
    ax = plt.scatter(data[:, 0], data[:, 1], c=data_ql[:], s=6, alpha=0.5, edgecolors='none', vmin=ql_min, vmax=ql_max)#, cmap = cm)
    if np.amax(data_ql[:]) > 0.0:
        plt.colorbar(ax, shrink=0.6)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    if type == 'norm':
        plt.title('data norm')
    else:
        plt.title('data')
    labeling(var_name1, var_name2, scale_thl, scale_qt)



    plt.subplot(2, 2, 4)
    color_arr = np.zeros(shape = np.size(samples_ql))
    for i in range(np.size(samples_ql)):
        color_arr[i] = samples_ql[i]*1e3
    print('samples_ql ', np.amin(samples_ql), np.amax(samples_ql), np.amin(color_arr), np.amax(color_arr),  np.shape(color_arr))
    try:
        print('try')
        # ax = plt.scatter(samples[:, 0], samples[:, 1], c=color_arr[:], s=6, alpha=0.5, edgecolors='none')
        ax = plt.scatter(samples[:, 0], samples[:, 1], c=samples_ql[:], s=6, alpha=0.5, edgecolors='none', vmin=ql_min, vmax=ql_max)#, cmap = cm)
        if np.amax(data_ql[:]) > 0.0:
            plt.colorbar(ax, shrink=0.6)
    except:
        try:
            print('except-try')
            ax = plt.scatter(samples[:, 0], samples[:, 1], c=color_arr[:], s=6, alpha=0.5, edgecolors='none', vmin=ql_min, vmax=ql_max)#, cmap = cm)
            # ax = plt.scatter(samples[:, 0], samples[:, 1], c=samples_ql[:], s=6, alpha=0.5, edgecolors='none')
            if np.amax(data_ql[:]) > 0.0:
                plt.colorbar(ax, shrink=0.6)
        except:
            print('except-except')
            ax = plt.scatter(samples[:, 0], samples[:, 1], s=6, alpha=0.5, edgecolors='none')#, cmap = cm)
            # pass
        # pass

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.title('samples')
    labeling(var_name1, var_name2, 0, 0)

    # plt.subplot(2, 3, 6)
    # color_arr = np.zeros(shape=np.shape(data_ql))
    # for i in range(np.shape(data_ql)[0]):
    #     color_arr[i] = data_ql[i] * 1e4
    # print('data_ql. ', np.amin(data_ql), np.amax(data_ql), np.amin(color_arr), np.amax(color_arr))
    # plt.scatter(data[:, 0], data[:, 1], c=color_arr[:], s=5, alpha=0.2)
    # # plt.gray()
    # if type == 'norm':
    #     plt.title('data norm')
    # else:
    #     plt.title('data')
    # labeling(var_name1, var_name2, scale_thl, scale_qt)

    plt.suptitle('Data & Samples: ncomp='+str(ncomp)+', z='+str(z)+'m', fontsize=18)
    if type == 'norm':
        savename = 'sample_figure_'+'ncomp'+str(ncomp)+'_norm_'+str(z)+'m.png'
    else:
        savename = 'sample_figure_' + 'ncomp' + str(ncomp) + '_' + str(z) + 'm.png'
    plt.savefig(
        os.path.join(path, 'CloudClosure_figures', savename))
    plt.close()
    return



def plot_hist(data_ql, path):
    plt.figure(figsize=(12,6))
    print('plot hist', data_ql.shape, data_ql[:,0].shape)

    bins = np.linspace(np.amin(data_ql), np.amax(data_ql), 1000)
    # print(np.amin(data_ql), np.amax(data_ql))
    # print('nonzero:', np.count_nonzero(data_ql))
    # # print(np.array(data_ql[:,0]))
    # # print(bins)
    # # print('shape ql: ', data_ql.shape)
    plt.subplot(1,2,1)
    time1 = time.clock()
    plt.hist(data_ql[:,0], bins)
    plt.xlim(0.0, 1e-6)
    # print('time hist:', time.clock() - time1)
    # plt.title('Histogram LES')
    # plt.xlabel('ql')
    plt.subplot(1, 2, 2)

    # bins = np.linspace(np.amin(data_ql), np.amax(data_ql), 10)
    # print(np.amin(data_ql), np.amax(data_ql))
    # plt.subplot(2,3,4)
    # print('shape ql: ', data_ql.shape)
    # time1 = time.clock()
    # plt.hist(data_ql, bins)
    # print('time hist:', time.clock()-time1)
    # plt.title('Histogram LES')
    # plt.xlabel('ql')
    #
    # plt.subplot(2,3,5)
    # plt.hist(samples_ql, bins)
    # plt.title('Histogram samples')
    # plt.xlabel('ql')

    # plt.show()

    savename = 'hist_figure.png'
    plt.savefig(
        os.path.join(path, 'CloudClosure_figures', savename))
    plt.close()
    return


def plot_error_vs_ncomp(error_array, ncomp_array, krange, dz, path):
    # error_array = (nk x ncomp):   array for relative errors
    # ncomp_array

    plt.figure()
    nk = np.shape(error_array)[0]
    for k in range(nk):
        plt.plot(ncomp_array, error_array[k,:], '-o', label='z='+str(krange[k]*dz))
    plt.legend()
    plt.xlabel('# Gaussian components in PDF')
    plt.ylabel('relative error in <ql>')
    plt.title('Relative Error in grid mean liquid water')
    savename = 'error_figure.png'
    plt.savefig(
        os.path.join(path, 'CloudClosure_figures', savename))
    plt.close()
    return





#----------------------------------------------------------------------
def labeling(var_name1, var_name2, scale_1, scale_2):
    # print('labeling', scale_1, scale_2)
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)

    if scale_1 != 0:
        plt.xlabel(var_name1 + ' (' + np.str(np.round(scale_1,3)) + ')')
    if scale_2 != 0:
        plt.ylabel(var_name2 + ' (' + np.str(np.round(scale_2,3)) + ')')

    return