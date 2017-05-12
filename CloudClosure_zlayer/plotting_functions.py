import os
import pylab as plt
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
from sklearn.preprocessing import StandardScaler
import time
import traceback


label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['xtick.direction']='out'
plt.rcParams['ytick.direction']='out'
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['figure.titlesize'] = 15
plt.rcParams['lines.linewidth'] = 3


def plot_PDF(data, data_norm, var_name1, var_name2, clf, dk, ncomp, error, z, path, save_name):
    print('Plot PDF: ', os.path.join(path+'_figures', 'PDF_figures_' + str(z) + 'm' + '_ncomp' + str(ncomp) + '.png'))

    xmin = np.amin(data[:,0])
    xmax = np.amax(data[:,0])
    ymin = np.amin(data[:,1])
    ymax = np.amax(data[:,1])

    xmin_norm = np.amin(data_norm[:, 0])
    xmax_norm = np.amax(data_norm[:, 0])
    ymin_norm = np.amin(data_norm[:, 1])
    ymax_norm = np.amax(data_norm[:, 1])

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
    plt.scatter(data[:,0], data[:,1], s=5, alpha=0.2)
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

    plt.suptitle('GMM, ncomp='+str(ncomp)+', error: '+str(error)+'    (z='+str(z)+')')

    # plt.savefig(
    #     os.path.join(path+'_figures', 'PDF_figures_' + str(z) + 'm' + '_ncomp' + str(ncomp) + '_dz'+str(dk)+'.png'))
    plt.savefig(os.path.join(path+'_figures', save_name+'.pdf'))
    # try:
    #     plt.savefig(os.path.join(path, 'CloudClosure_z_figures','PDF_figures_'+str(z)+'m'+'_ncomp'+str(ncomp)+'.png'))
    # except:
    #     print('!!!!! figure with ncomp='+str(ncomp)+', z='+str(z)+' not saved !!!!!')
    plt.close()
    return


def scatter_data(data0, data1, var_name1, var_name2, dk, ncomp, error, z, path, save_name):
    print('Plot PDF: ', os.path.join(path+'_figures', 'PDF_figures_' + str(z) + 'm' + '_ncomp' + str(ncomp) + '.png'))

    xmin = np.amin(data0)
    xmax = np.amax(data0)
    ymin = np.amin(data1)
    ymax = np.amax(data1)

    plt.figure(figsize=(16,8))
    plt.subplot(1,4,1)
    plt.scatter(data0, 1, s=5, alpha=0.2)
    plt.title('data')
    labeling(var_name1, var_name2, xmin, xmax, ymin, ymax)

    plt.suptitle('GMM, ncomp='+str(ncomp)+', error: '+str(error)+'    (z='+str(z)+')')

    # plt.savefig(
    #     os.path.join(path+'_figures', 'PDF_figures_' + str(z) + 'm' + '_ncomp' + str(ncomp) + '_dz'+str(dk)+'.png'))
    plt.savefig(os.path.join(path+'_figures', save_name+'.pdf'))
    plt.close()
    return



def plot_abs_error(error_ql, ql_mean_ref, error_cf, cf_ref, n_sample, ncomp_array, krange, dz, dk, path):
    cmap = cm.get_cmap('jet')
    nk = error_ql.shape[0]

    plt.figure(figsize=(18,6))
    for k in range(nk):
        col = cmap(np.double(k)/nk)

        plt.subplot(1,2,1)
        # ql_ref = ql_mean_ref[k] * np.ones(len(ncomp_array))
        plt.plot(ncomp_array, np.abs(error_ql[k,:]), '-o', color=col, label=r'$\epsilon(q_l)$, z='+str(krange[k]*dz))
        # plt.plot(ncomp_array, ql_ref, '-', color=col, linewidth=1, label=r'<ql>$_{field}$=' + str(ql_mean_ref[k]))

        plt.subplot(1,2,2)
        # ref = cf_ref[k] * np.ones(len(ncomp_array))
        plt.plot(ncomp_array, np.abs(error_cf[k, :]), '-o', color=col, label=r'$\epsilon$(CF), z=' + str(np.round(krange[k] * dz,4)))
        # plt.plot(ncomp_array, ref, '-', color=col, linewidth=1, label=r'CF$_{field}$=' + str(cf_ref[k]))

    plt.subplot(1, 2, 1)
    plt.legend(loc='best', bbox_to_anchor=(0, 1))
    plt.xlabel('# Gaussian components in PDF')
    plt.ylabel(r'abs(<ql>$_PDF$-<ql>$_{field}$)')
    plt.title('Absolute Error domain mean liquid water (n_sample: ' + str(n_sample) + ')')

    plt.subplot(1, 2, 2)
    plt.legend(loc='best', bbox_to_anchor=(0, 1))
    plt.xlabel('# Gaussian components in PDF')
    plt.ylabel(r'abs(CF$_{PDF}$-CF$_{field}$)')
    plt.title('Absolute Error cloud fraction')

    save_name = 'error_figure_abs_dz'+str(dk)
    plt.savefig(os.path.join(path + '_figures', save_name + '.pdf'))
    plt.close()
    return



def plot_error_vs_ncomp_ql(error_ql, rel_error_ql, n_sample, ql_mean_ref, ql_mean_comp, ncomp_array, krange, dz, dk, path):
    cmap = cm.get_cmap('jet')
    nk = error_ql.shape[0]

    plt.figure(figsize=(18,9))
    for k in range(nk):
        col = cmap(np.double(k)/nk)

        plt.subplot(1, 2, 1)
        plt.plot(ncomp_array, rel_error_ql[k, :], '-o', color=col,
                 label='z=' + str(krange[k] * dz))
        plt.plot(ncomp_array, np.zeros(len(ncomp_array)), 'k-', linewidth=1)

        plt.subplot(1, 2, 2)
        plt.plot(ncomp_array, error_ql[k, :], '-o', color=col,
                 label='z=' + str(krange[k] * dz))

        ql_ref = -ql_mean_ref[k] * np.ones(len(ncomp_array))
        plt.plot(ncomp_array, ql_ref, '-', color=col, linewidth=1, label=r'-<ql>$_{field}$='+str(ql_mean_ref[k]))

    plt.subplot(1, 2, 1)
    plt.legend(loc='best', bbox_to_anchor=(0, 1))
    plt.xlabel('# Gaussian components in PDF')
    plt.ylabel('relative error in <ql>')
    plt.title('Relative Error domain mean liquid water (n_sample: ' + str(n_sample) + ')')

    plt.subplot(1, 2, 2)
    plt.legend(loc='best', bbox_to_anchor=(0, 1))
    plt.xlabel('# Gaussian components in PDF')
    plt.ylabel('absolute error in <ql>')
    plt.title('Absolute Error domain mean liquid water')

    save_name = 'error_figure_ql_dz'+str(dk)
    plt.savefig(os.path.join(path + '_figures', save_name + '.pdf'))
    plt.close()
    return



def plot_error_vs_ncomp_cf(error_cf, rel_error_cf, n_sample, cf_ref, ncomp_array, krange, dz, dk, path):
    cmap = cm.get_cmap('jet')
    nk = np.shape(error_cf)[0]
    ncomp = len(ncomp_array)
    plt.figure(figsize=(12, 6))
    for k in range(nk):
        col = cmap(np.double(k) / nk)

        plt.subplot(1, 2, 1)
        plt.plot(ncomp_array, rel_error_cf[k, :], '-o', color=col, label='z=' + str(krange[k] * dz))

        plt.subplot(1, 2, 2)
        plt.plot(ncomp_array, error_cf[k, :], '-o', color=col, label='z=' + str(krange[k] * dz))

        ref = -cf_ref[k]*np.ones(ncomp)
        plt.plot(ncomp_array, ref, '-', color=col, linewidth=1, label=r'-CF$_{field}$='+str(cf_ref[k]))

    plt.subplot(1, 2, 1)
    plt.legend(loc='best', bbox_to_anchor=(0, 1))
    plt.xlabel('# Gaussian components in PDF')
    plt.ylabel('relative error in CF')
    plt.title('Relative Error Cloud Fraction (n_sample: '+str(n_sample)+')')

    plt.subplot(1, 2, 2)
    # plt.legend(loc=1)
    plt.legend(loc='best', bbox_to_anchor=(0, 1))
    plt.xlabel('# Gaussian components in PDF')
    plt.ylabel('absolute error in CF')
    plt.title('Absolute Error Cloud Fraction')

    save_name = 'error_figure_cf_dz'+str(dk)
    plt.savefig(os.path.join(path + '_figures', save_name + '.pdf'))
    plt.close()
    return




def plot_PDF_components(means_, covars_, weights_, ncomp, krange, dz, dk, path):
    # means_ = (nk, ncomp, nvar)
    # covariance_ = (nk, ncomp, nvar, nvar)
    # weights_ = (nk, ncomp)

    plt.figure(figsize=(12,8))
    # for ncomp in ncomp_range:
    for n in range(ncomp):
        plt.subplot(2, 3, 1)
        plt.plot(means_[:,n,0],krange, '-o', label='mean(th_l), n_comp='+str(n))
        plt.subplot(2, 3, 4)
        plt.plot(means_[:, n, 1], krange, '-o', label='mean(qt), n_comp=' + str(n))

        plt.subplot(2, 3, 2)
        plt.plot(covars_[:, n, 0, 0], krange, '-o', label='var(th_l), n_comp=' + str(n))
        plt.subplot(2, 3, 5)
        plt.plot(covars_[:, n, 1, 1], krange, '-o', label='var(qt), n_comp=' + str(n))

        plt.subplot(2, 3, 3)
        plt.plot(weights_[:, n], krange, '-o', label='weight, n_comp='+str(n))

    plt.subplot(2, 3, 1)
    plt.title(r'<$\theta_l$>')
    plt.legend(loc='best')
    plt.subplot(2, 3, 4)
    plt.title(r'<$q_t$>')
    plt.legend(loc='best')
    plt.subplot(2, 3, 2)
    plt.title(r'Var[$\theta_l$]')
    plt.legend(loc='best')
    plt.subplot(2, 3, 5)
    plt.title(r'Var[$q_t]')
    plt.legend(loc='best')
    plt.subplot(2, 3, 3)
    plt.title('weights')
    plt.legend(loc='best')


    save_name = 'figure_PDFcomps_ncomp'+str(ncomp)+'_dk'+str(dk)
    plt.savefig(os.path.join(path + '_figures', save_name + '.pdf'))

    return







def labeling(var_name1, var_name2, xmin, xmax, ymin, ymax):
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    return