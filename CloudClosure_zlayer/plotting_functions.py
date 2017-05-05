import os
import pylab as plt
import numpy as np
import matplotlib.mlab as mlab
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
plt.rcParams['legend.fontsize'] = 10
# plt.rcParams['title.fontsize'] = 15


def plot_PDF(data, var_name1, var_name2, dk, ncomp, error, z, path, save_name):
    print('Plot PDF: ', os.path.join(path+'_figures', 'PDF_figures_' + str(z) + 'm' + '_ncomp' + str(ncomp) + '.png'))

    xmin = np.amin(data[:,0])
    xmax = np.amax(data[:,0])
    ymin = np.amin(data[:,1])
    ymax = np.amax(data[:,1])

    # xmin_norm = np.amin(data_norm[:, 0])
    # xmax_norm = np.amax(data_norm[:, 0])
    # ymin_norm = np.amin(data_norm[:, 1])
    # ymax_norm = np.amax(data_norm[:, 1])

    plt.figure(figsize=(16,8))
    plt.subplot(1,4,1)
    plt.scatter(data[:,0], data[:,1], s=5, alpha=0.2)
    plt.title('data (n='+str(data.shape[0])+')' )
    labeling(var_name1, var_name2, xmin, xmax, ymin, ymax)
    # plt.xlim(xmin, xmax)
    # plt.ylim(ymin, ymax)


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





def labeling(var_name1, var_name2, xmin, xmax, ymin, ymax):
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    return