import os
import pylab as plt
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
# from sklearn.preprocessing import StandardScaler
import time
# import traceback


label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['xtick.direction']='out'
plt.rcParams['ytick.direction']='out'
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['figure.titlesize'] = 10
plt.rcParams['lines.linewidth'] = 0.8




def plot_PDF(data, data_norm, ql_all, var_name1, var_name2, clf, dk, ncomp, error, z, Lx, path, save_name):
    print('Plot PDF: ', os.path.join(path+'_figures', 'PDF_figures_' + str(z) + 'm' + '_ncomp' + str(ncomp) + '.png'))

    xmin = np.amin(data[:,0])
    xmax = np.amax(data[:,0])
    ymin = np.amin(data[:,1])
    ymax = np.amax(data[:,1])

    xmin_norm = np.amin(data_norm[:, 0])
    xmax_norm = np.amax(data_norm[:, 0])
    ymin_norm = np.amin(data_norm[:, 1])
    ymax_norm = np.amax(data_norm[:, 1])

    ql_min = np.amin(ql_all)
    ql_max = np.amax(ql_all)

    n_sample = np.int(1e2)
    nvar = 2
    x_ = np.linspace(xmin_norm, xmax_norm, n_sample)
    y_ = np.linspace(ymin_norm, ymax_norm, n_sample)
    XX_ = np.ndarray(shape=(n_sample ** nvar, nvar))

    # (1) print PDF computed by GMM
    delta_i = n_sample
    for i in range(n_sample):
        for j in range(n_sample):
            shift = i * delta_i + j
            XX_[shift, 0] = x_[i]
            XX_[shift, 1] = y_[j]
    ZZ_ = clf.score_samples(XX_)
    ZZ = np.ndarray(shape=(n_sample, n_sample))
    for j in range(n_sample ** nvar):
        jshift = np.mod(j, delta_i)
        ishift = (j - jshift) / delta_i
        ZZ[ishift, jshift] = ZZ_[j]

    plt.figure(figsize=(16,8))
    plt.subplot(1,4,1)
    # plt.scatter(data[:,0], data[:,1], s=5, alpha=0.2)
    plt.scatter(data[:,0], data[:,1], c = ql_all[:], s = 5, alpha = 0.2, edgecolors = 'none', vmin = ql_min, vmax = ql_max)
    plt.title('data (n='+str(data.shape[0])+')' )
    labeling(var_name1, var_name2, xmin, xmax, ymin, ymax)

    plt.subplot(1, 4, 2)
    ax = plt.contour(x_, y_, np.exp(ZZ).T)
    plt.scatter(data_norm[:, 0], data_norm[:, 1], s=5, alpha=0.2)
    plt.colorbar(ax)
    plt.title('data normalised')
    labeling(var_name1, var_name2, xmin_norm, xmax_norm, ymin_norm, ymax_norm)

    plt.subplot(1, 4, 3)
    # ax = plt.contour(x_, y_, ZZ.T)
    ax = plt.contourf(x_, y_, np.exp(ZZ).T)
    plt.colorbar(ax)
    # plt.scatter(data_norm[:, 0], data_norm[:, 1], s=5, alpha=0.2)
    plt.title('PDF(thl, qt)')
    labeling(var_name1, var_name2, xmin_norm, xmax_norm, ymin_norm, ymax_norm)

    plt.subplot(1, 4, 4)
    # # ax = plt.contourf(x_, y_, np.exp(ZZ).T)
    # # plt.colorbar(ax)
    try:
        ax = plt.contourf(x_, y_, np.exp(ZZ).T, norm=LogNorm())
    except:
        print('except in: plot_PDF, subplot(1,4,3)')
    # try:
    #     plt.colorbar(ax)
    # except:
    #     pass
    labeling(var_name1, var_name2, xmin_norm, xmax_norm, ymin_norm, ymax_norm)
    plt.title('PDF(thl, qt)')

    plt.suptitle('GMM, ncomp='+str(ncomp)+ ': error: '+str(error)+'    (Lx=' + str(Lx) + 'z='+str(z)+')')
    # plt.savefig(os.path.join(path+'_figures', 'PDF_figures_' + str(z) + 'm' + '_ncomp' + str(ncomp) + '_dz'+str(dk)+'.png'))
    # plt.savefig(os.path.join(path+'_figures', save_name+'.pdf'))
    plt.savefig(os.path.join(path + '_figures', save_name + '.png'))
    plt.close()
    return



def labeling(var_name1, var_name2, xmin, xmax, ymin, ymax):
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    return