import pylab as plt

#----------------------------------------------------------------------
# def plot_PDF_samples_qt(data, data_aux, var_name1, var_name2, clf, clf_aux, time, z):
def plot_PDF_samples_qt(data, var_name1, var_name2, clf, time, z):
    global ncomp
    amp = 1e2
    print('')
    print('plot PDF samples qt, factor='+np.str(amp))

    import matplotlib.mlab as mlab
    import matplotlib.cm as cm

    data_aux = np.ndarray(shape=((nx * ny), nvar))
    data_aux[:, 0] = data[:, 0]
    data_aux[:, 1] = data[:, 1] * amp
    clf_aux = mixture.GaussianMixture(n_components=ncomp, covariance_type='full')
    clf_aux.fit(data_aux)

    det_ = np.linalg.det(clf.covariances_[0, :, :])
    fact_ = 1. / np.sqrt((2 * np.pi) ** 2 * det_)
    det_aux_ = np.linalg.det(clf_aux.covariances_[0, :, :])
    fact_aux_ = 1. / np.sqrt((2 * np.pi) ** 2 * det_aux_)

    # Plotting
    n_sample = 100
    x1_max = np.amax(data[:, 0])
    x1_min = np.amin(data[:, 0])
    x2_max = np.amax(data[:, 1])
    x2_min = np.amin(data[:, 1])
    x = np.linspace(x1_min, x1_max, n_sample)
    y = np.linspace(x2_min, x2_max, n_sample)
    X, Y = np.meshgrid(x, y)
    XX = np.array([X.ravel(), Y.ravel()]).T
    Z = clf.score_samples(XX).reshape(X.shape)
    x_aux = np.linspace(np.amin(data_aux[:,0]), np.amax(data_aux[:,0]), n_sample)
    y_aux = np.linspace(np.amin(data_aux[:, 1]), np.amax(data_aux[:, 1]), n_sample)
    X_aux, Y_aux = np.meshgrid(x_aux, y_aux)
    XX_aux = np.array([X_aux.ravel(), Y_aux.ravel()]).T
    Z_aux = clf_aux.score_samples(XX_aux).reshape(X_aux.shape)

    # mx1 = clf.means_[0, 0]
    # my1 = clf.means_[0, 1]
    # sx1 = np.sqrt(clf.covariances_[0, 0, 0])
    # sy1 = np.sqrt(clf.covariances_[0, 1, 1])
    # sxy1 = clf.covariances_[0, 1, 0]
    # Z1 = mlab.bivariate_normal(X, Y, sigmax=sx1, sigmay=sy1, mux=mx1, muy=my1, sigmaxy=sxy1)
    # mx1_aux = clf_aux.means_[0, 0]
    # my1_aux = clf_aux.means_[0, 1]
    # sx1_aux = np.sqrt(clf_aux.covariances_[0, 0, 0])
    # sy1_aux = np.sqrt(clf_aux.covariances_[0, 1, 1])
    # sxy1_aux = clf_aux.covariances_[0, 1, 0]
    # Z1_aux = mlab.bivariate_normal(X_aux, Y_aux, sigmax=sx1_aux, sigmay=sy1_aux, mux=mx1_aux, muy=my1_aux, sigmaxy=sxy1_aux)

    fig = plt.figure(figsize=(12, 16))
    plt.subplot(3, 2, 1)
    plt.scatter(data[:, 0], data[:, 1], s=5, alpha=0.2)
    ax1 = plt.contour(X, Y, Z, colors='w', levels=np.linspace(10, 20, 11))
    plt.colorbar(ax1,shrink=0.8)
    plt.xlim([x1_min,x1_max])
    plt.ylim([x2_min, x2_max])
    plt.title(var_name1 + var_name2 + ' (data), t=' + str(time) + ', z=' + str(z))
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    plt.subplot(3, 2, 2)
    plt.scatter(data[:, 0], data[:, 1], s=5, alpha=0.2)
    ax1 = plt.contour(X_aux, Y_aux, Z_aux, colors='w', levels=np.linspace(10, 20, 11))
    plt.colorbar(ax1, shrink=0.8)
    plt.xlim([x1_min, x1_max])
    plt.ylim([x2_min, x2_max])
    plt.title(var_name1 + var_name2 + ' (data), t=' + str(time) + ', z=' + str(z))
    labeling(var_name1, var_name2, amp)

    plt.subplot(3, 2, 3)
    plt.scatter(data[:, 0], data[:, 1], s=5, alpha=0.2)
    levels = np.linspace(0, fact_, 10)
    ax1 = plt.contour(X, Y, np.exp(Z), levels=levels, linewidths=2)
    plt.plot([clf.means_[0, 0]], [clf.means_[0, 1]], 'o', markersize=8)
    plt.colorbar(ax1, shrink=0.8)
    plt.title('EM PDF')
    plt.xlim([x1_min, x1_max])
    plt.ylim([x2_min, x2_max])
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    plt.subplot(3, 2, 4)
    plt.scatter(data_aux[:, 0], data_aux[:, 1], s=5, alpha=0.2)
    # ax1 = plt.contour(X_aux, Y_aux, np.exp(Z_aux), levels=levels, linewidths=2)
    ax1 = plt.contour(X_aux, Y_aux, np.exp(Z_aux), linewidths=2)
    plt.plot([clf_aux.means_[0, 0]], [clf_aux.means_[0, 1]], 'o', markersize=8)
    plt.colorbar(ax1, shrink=0.8)
    plt.title('EM PDF')
    labeling(var_name1, var_name2, amp)

    plt.subplot(3, 2, 5)
    ax1 = plt.contourf(X, Y, np.exp(Z))
    # plt.scatter(X, Y, s=2, alpha=0.5)
    plt.colorbar(ax1, shrink=0.8)
    plt.title('EM PDF')
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    plt.subplot(3, 2, 6)
    ax1 = plt.contourf(X_aux, Y_aux, np.exp(Z_aux))
    plt.colorbar(ax1, shrink=0.8)
    plt.title('EM PDF')
    labeling(var_name1, var_name2, amp)

    fig.suptitle('Cloud Closure: Univariate Gaussian PDF fit', fontsize=20)
    plt.savefig(
        os.path.join(
            fullpath_out,'CloudClosure_alltimes_figures/CC_' + var_name1 + '_' + var_name2 + '_z'
                         + str(np.int(z)) + 'm_bivariate_alltime.png')
    )

    plt.close()
    return
#----------------------------------------------------------------------
def plot_PDF_samples(data, var_name1, var_name2, clf, time, z):
    import matplotlib.mlab as mlab
    import matplotlib.cm as cm

    det_ = np.linalg.det(clf.covariances_[0, :, :])
    fact_ = 1. / np.sqrt((2 * np.pi) ** 2 * det_)

    # Plotting
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
    mx1 = clf.means_[0, 0]
    my1 = clf.means_[0, 1]
    sx1 = np.sqrt(clf.covariances_[0, 0, 0])
    sy1 = np.sqrt(clf.covariances_[0, 1, 1])
    sxy1 = clf.covariances_[0, 1, 0]
    Z1 = mlab.bivariate_normal(X, Y, sigmax=sx1, sigmay=sy1, mux=mx1, muy=my1, sigmaxy=sxy1)

    plt.figure(figsize=(12, 12))
    levels_tot = np.linspace(0, fact_, 10)
    if fact_ <= 2:
        levels_cont = np.arange(0, fact_, 0.2)
        levels_contf = np.arange(0, fact_, 0.2)
    elif fact_ <= 10:
        levels_cont = np.arange(0,fact_,0.5)
        levels_contf = np.arange(0,fact_,0.5)
    elif fact_ <= 20:
        levels_cont = np.arange(0,fact_,2)
        levels_contf = np.arange(0,fact_,2)
    elif fact_ <= 50:
        levels_cont = np.arange(0,fact_,5)
        levels_contf = np.arange(0,fact_,5)
    else:
        levels_cont = np.arange(0,fact_,20)
        levels_contf = np.arange(0,fact_,20)
    levels_comp = np.linspace(0, fact_, 7)

    plt.subplot(3, 2, 1)
    plt.scatter(data[:, 0], data[:, 1], s=2, alpha=0.05)
    ax1 = plt.contour(X, Y, Z, levels=np.linspace(10, 20, 2))
    plt.colorbar(ax1,shrink=0.8)
    plt.xlim([x1_min,x1_max])
    plt.ylim([x2_min, x2_max])
    plt.title(var_name1 + var_name2 + ' (data), t=' + str(time) + ', z=' + str(z))
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    plt.colorbar(ax1,shrink=0.8)
    plt.xlim([x1_min,x1_max])
    plt.ylim([x2_min, x2_max])
    plt.title(var_name1 + var_name2 + ' (data), t=' + str(time) + ', z=' + str(z))
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)

    plt.subplot(3, 2, 3)
    ax1 = plt.hist2d(data[:, 0], data[:, 1], bins=30, normed=True)
    plt.colorbar(shrink=0.8)
    ax2 = plt.contour(X, Y, np.exp(Z), levels=levels_cont, linewidths=1, colors='w')
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    plt.title('data histogram')
    plt.subplot(3, 2, 4)
    ax1 = plt.contourf(X, Y, np.exp(Z),levels=levels_contf)
    plt.colorbar(ax1, shrink=0.8)
    plt.title('EM PDF')
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    plt.subplot(3, 2, 5)
    plt.scatter(data[:, 0], data[:, 1], s=2, alpha=0.05)
    ax1 = plt.contour(X, Y, Z1, levels=levels_comp, linewidths=1.5)
    plt.plot([clf.means_[0, 0]], [clf.means_[0, 1]], 'wo', markersize=6)
    plt.colorbar(ax1, shrink=0.8)
    plt.title('f1, f2')
    plt.xlim([x1_min, x1_max])
    plt.ylim([x2_min, x2_max])
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    plt.subplot(3, 2, 6)
    plt.scatter(data[:, 0], data[:, 1], s=2, alpha=0.05)
    ax1 = plt.contour(X, Y, np.exp(Z), levels=levels_cont, linewidths=1.5)
    # ax1 = plt.contour(X, Y, Z1+Z2, linewidths=1.5)
    plt.plot([clf.means_[0, 0]], [clf.means_[0, 1]], 'wo', markersize=6)
    plt.colorbar(ax1, shrink=0.8)
    plt.xlim([x1_min, x1_max])
    plt.ylim([x2_min, x2_max])
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    plt.title('f = f1 + f2')

    plt.savefig(fullpath_out+'CloudClosure_alltimes_figures/CC_bivariate_' + var_name1 + '_' + var_name2 + '_z' + str(
        np.int(z)) + 'm_alltime.png')

    plt.close()
    return