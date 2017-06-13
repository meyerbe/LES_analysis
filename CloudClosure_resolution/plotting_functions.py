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

# # INPUTS:
    # - error_ql = (nk): array with absolute error in mean liquid water ql
    # - ql_mean_ref = (nk): array with mean liquid water ql values from LES field
    # - error_cf = (nk): array with absolute error in cloud fraction
    # - cf_ref = (nk): array with cloud fraction from LES field
    # - n_sample: # of samples in Monte Carlos simulation
    # - ncomp_array: array with number of Gaussian components for which GMM is computed
    # - krange: arrays with z-values at which the PDF model is computed
    # - dz = nml['dz']
    # - Lx: LES domain size for which PDF model is computed
    # - dk: number of layers over which data are accumulated
    # - path, where figures are saved
    #
    # - nk: length krange



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
    # plt.savefig(os.path.join(path+'_figures', save_name+'.pdf'))
    plt.savefig(os.path.join(path + '_figures', save_name + '.png'))
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
    # plt.savefig(os.path.join(path+'_figures', save_name+'.pdf'))
    plt.savefig(os.path.join(path + '_figures', save_name + '.png'))
    plt.close()
    return


def plot_abs_error(error_ql, ql_mean_ref, error_cf, cf_ref, n_sample, ncomp_array, krange, dz, Lx, dk, path):
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

    save_name = 'error_figure_abs_Lx' + str(Lx) + '_dk' + str(dk)
    plt.savefig(os.path.join(path + '_figures', save_name + '.pdf'))
    plt.close()
    return



def plot_error_vs_ncomp_ql(error_ql, rel_error_ql, n_sample, ql_mean_ref, ql_mean_comp, ncomp_array, krange, dz, Lx, dk, path):
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

    save_name = 'error_figure_ql_Lx'+str(Lx)+'_dk'+str(dk)
    plt.savefig(os.path.join(path + '_figures', save_name + '.pdf'))
    plt.close()
    return



def plot_error_vs_ncomp_cf(error_cf, rel_error_cf, n_sample, cf_ref, ncomp_array, krange, dz, Lx, dk, path):
    cmap = cm.get_cmap('jet')
    nk = np.shape(error_cf)[0]
    ncomp = len(ncomp_array)
    plt.figure(figsize=(18, 9))
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

    save_name = 'error_figure_cf_Lx'+str(Lx)+'_dk'+str(dk)
    plt.savefig(os.path.join(path + '_figures', save_name + '.pdf'))
    plt.close()
    return




def plot_PDF_components(means_, covars_, weights_, ncomp, krange, dz, Lx, dk, path):
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


    save_name = 'figure_PDFcomps_ncomp'+str(ncomp)+'_Lx'+str(Lx)+'_dz'+str(dk)
    # plt.savefig(os.path.join(path + '_figures', save_name + '.pdf'))
    plt.savefig(os.path.join(path + '_figures', save_name + '.png'))

    return





def plot_samples(type, data, data_ql, samples, samples_ql, var_name1, var_name2, scaler, ncomp, z, path, save_name):
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
    xmin = np.min([np.amin(data[:, 0]),np.amin(samples[:, 0])])
    xmax = np.max([np.amax(data[:, 0]),np.amax(samples[:, 0])])
    ymin = np.min([np.amin(data[:, 1]),np.amin(samples[:, 1])])
    ymax = np.max([np.amax(data[:, 1]),np.amax(samples[:, 1])])
    cm = plt.cm.get_cmap('RdYlBu')
    cm = plt.cm.get_cmap('viridis')

    plt.figure(figsize=(8,9))
    plt.subplot(2,2,1)
    plt.scatter(data[:,0], data[:,1], s=5, alpha=0.2)
    if type == 'norm':
        plt.title('data norm (n='+str(np.shape(data)[0])+')')
    else:
        plt.title('data (n='+str(np.shape(data)[0])+')')
    labeling(var_name1, var_name2, xmin, xmax, ymin, ymax)


    plt.subplot(2,2,2)
    plt.scatter(samples[:,0], samples[:,1], s=5, alpha=0.2)
    plt.title('samples (n='+str(np.size(samples_ql))+')')
    labeling(var_name1, var_name2, xmin, xmax, ymin, ymax)

    plt.subplot(2, 2, 3)
    try:
        ax = plt.scatter(data[:, 0], data[:, 1], c=data_ql[:], s=6, alpha=0.5, edgecolors='none', vmin=ql_min, vmax=ql_max)#, cmap = cm)
        if np.amax(data_ql[:]) > 0.0:
            plt.colorbar(ax, shrink=0.6)
    except:
        print('except data - color arr: ', data.shape, data_ql.shape)
        # # traceback.print_exc()
        # color_arr = np.zeros(shape=np.size(data_ql))
        # for i in range(np.size(data_ql)):
        #     color_arr[i] = data_ql[i] * 1e3
        # print('data_ql except', np.amin(data_ql), np.amax(data_ql), np.amin(color_arr), np.amax(color_arr),np.shape(color_arr))
        # print(color_arr.shape)
        # print(color_arr)
        # try:
        #     ax = plt.scatter(data[:, 0], data[:, 1], c=color_arr[:], s=6, alpha=0.5, edgecolors='none', vmin=ql_min, vmax=ql_max)  # , cmap = cm)
        #     if np.amax(data_ql[:]) > 0.0:
        #         plt.colorbar(ax, shrink=0.6)
        # except:
        #     # traceback.print_exc()
        #     print('except except (data ql)')
        #     ax = plt.scatter(data[:, 0], data[:, 1], s=6, alpha=0.5, edgecolors='none', vmin=ql_min,vmax=ql_max)  # , cmap = cm)
    if type == 'norm':
        plt.title('data norm (n=' + str(np.shape(data)[0]) + ')')
    else:
        plt.title('data (n=' + str(np.shape(data)[0]) + ')')
    labeling(var_name1, var_name2, xmin, xmax, ymin, ymax)



    plt.subplot(2, 2, 4)
    color_arr = np.zeros(shape=np.size(samples_ql))
    for i in range(np.size(samples_ql)):
        color_arr[i] = samples_ql[i] * 1e3
    print('samples_ql ', np.amin(samples_ql), np.amax(samples_ql), np.amin(color_arr), np.amax(color_arr),
          np.shape(color_arr))
    try:
        print('try')
        ax = plt.scatter(samples[:, 0], samples[:, 1], c=samples_ql[:], s=6, alpha=0.5, edgecolors='none', vmin=ql_min, vmax=ql_max)#, cmap = cm)
        if np.amax(data_ql[:]) > 0.0:
            plt.colorbar(ax, shrink=0.6)
    except:
        pass
        # try:
        #     print('except-try')
        #     ax = plt.scatter(samples[:, 0], samples[:, 1], c=color_arr[:], s=6, alpha=0.5, edgecolors='none', vmin=ql_min, vmax=ql_max)#, cmap = cm)
        #     if np.amax(data_ql[:]) > 0.0:
        #         plt.colorbar(ax, shrink=0.6)
        # except:
        #     print('except-except')
        #     ax = plt.scatter(samples[:, 0], samples[:, 1], s=6, alpha=0.5, edgecolors='none')#, cmap = cm)

    plt.title('samples (n=' + str(np.size(samples_ql)) + ')')
    labeling(var_name1, var_name2, xmin, xmax, ymin, ymax)

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
    # labeling(var_name1, var_name2, xmin, xmax, ymin, ymax)

    plt.suptitle('Data & Samples: ncomp='+str(ncomp)+', z='+str(z)+'m', fontsize=18)
    # if type == 'norm':
    #     savename = 'sample_figure_'+'ncomp'+str(ncomp)+'_norm_'+str(z)+'m.png'
    # else:
    #     savename = 'sample_figure_' + 'ncomp' + str(ncomp) + '_' + str(z) + 'm.png'
    print(os.path.join(path + '_figures', save_name))
    plt.savefig(
        os.path.join(path + '_figures', save_name))
    plt.close()
    return



def labeling(var_name1, var_name2, xmin, xmax, ymin, ymax):
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    return