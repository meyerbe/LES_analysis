import os
import pylab as plt
import numpy as np
import matplotlib.mlab as mlab
from sklearn.preprocessing import StandardScaler

def plot_PDF(data, data_norm, var_name1, var_name2, nvar, clf, scaler, ncomp, error, z, path):
    print('plotting')

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

    XXii = np.zeros(shape=XX_.shape)
    print('XX before: ', ZZ.shape, XX_.shape)
    print(XX_[0:3, 0])
    print(XX_[0:3, 1])
    print('XX after inverse')
    XXi = scaler.inverse_transform(XX_)
    print(XXi[0:3,0])
    print(XXi[0:3, 1])
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
    plt.scatter(data[:,0], data[:,1],s=5, alpha=0.2)
    plt.title('data')
    labeling(var_name1,var_name2,1,1)
    plt.xlim(np.amin(data[:,0]), np.amax(data[:,0]))
    plt.ylim(np.amin(data[:,1]), np.amax(data[:,1]))

    plt.subplot(1, 3, 2)
    plt.scatter(data_norm[:, 0], data_norm[:, 1], s=5, alpha=0.2)
    plt.title('data norm')
    labeling(var_name1, var_name2, scaler.scale_[0], scaler.scale_[1])
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

    plt.subplot(1,3,3)
    ax = plt.contour(x_,y_,np.exp(ZZ).T)
    plt.scatter(data_norm[:, 0], data_norm[:, 1], s=5, alpha=0.2)
    plt.title('data norm')
    labeling(var_name1, var_name2, scaler.scale_[0], scaler.scale_[1])
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

    plt.suptitle('GMM, ncomp='+str(ncomp)+', error: '+str(error)+'    (z='+str(z)+')')

    plt.savefig(os.path.join(path,'CloudClosure_figures','test_figures_'+'ncomp'+str(ncomp)+'_'+str(z)+'m.png'))
    # plt.show()
    plt.close()
    return












#----------------------------------------------------------------------
def labeling(var_name1, var_name2, scale_thl, scale_qt):
    print('labeling', scale_thl, scale_qt)
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    if var_name1 == 'qt':
        plt.xlabel(var_name1 + ' (' + np.str(np.round(scale_qt,2)) + ')')
        plt.ylabel(var_name2)
    else:
        plt.xlabel(var_name1)
        plt.ylabel(var_name2 + ' (' + np.str(np.round(scale_thl,2)) + ')')

    return