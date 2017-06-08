import os
import json as simplejson
import time
import netCDF4 as nc
import numpy as np
cimport numpy as np
from sklearn import mixture
from sklearn.preprocessing import StandardScaler
import pickle


import CC_thermodynamics_c
from CC_thermodynamics_c cimport LatentHeat, ClausiusClapeyron
from CC_thermodynamics_c import sat_adj_fromentropy, sat_adj_fromthetali
from plotting_functions import plot_labels_comparison, plot_labels, plot_labels_test, plot_parameters_pdf

# NEW:
#   - no more accumlation over times (files) or levels for PDF model
#   - no more looping over components

# class Updrafts:
# (1) PDF model
#   (a) Updrafts.update: read in fields, compute PDF, label points
#       (i) read in 3D files
#       (ii) self.predict_PDF(...) --> output: clf (PDF) & labels (2D array for all data points); error for <ql> and CF
#               (A) Set up files
#               (B) Compute PDF
#                   - read in fields (and accumulate over several time steps)
#                   for each k...
#                       - normalise data
#                       - compute PDF for two components
#                       - compute labels
#                       - sort labels and PDF parameters according to <qt>
#                       - rearrange labels into 3D array
#               (C) Compute Error
#               (D) IO:
#                       - output updrafts file (labels)
#                       - output PDF parameters (sorted)
#                       - output errors
#           calls:
#               self.sort_PDF(...)              --> sorts PDF parameters
#               self.create_updrafts_file()
#               self.write_updrafts_file()      --> output of labels
#               self.create_statistics_file()   --> output of sorted PDF parameters
#       (iii) plotting functions
#
#   (b) read in labelled fields directly
#
# (2) Tracer model: read in pkl-files from Colleen
#       self.read_in_updrafts_colleen(...) --> output: python arrays with labels
#
# (3) scatter plot (a) PDF-labelling, (b) Colleen labeling (4 types), (c) compare labels --> green if both, red if only PDF, blue if only Colleen (separately for each label)


# TODO:
# - read in z-profile from Colleens files (..._updrafts.pkl)
# - find same z-value as in PDF_labels
# - send this 2d section to plotting
# - !!!! need same data for mine and for Colleens sampling!! (her data = 192x192, 3-times as large iin horizontal)


cdef class Updrafts:
    def __init__(self):
        self.path_ref = ' '
        self.p_ref = None
        self.z_ref = None

        self.krange = None
        self.zrange = None
        self.times_tracers = None
        return



    cpdef initialize(self, krange, path, case_name):
        print('----- Updraft Clustering: initialize -----')

        # (A) PDF Model
        nml = simplejson.loads(open(os.path.join(path, case_name+'.in')).read())
        cdef:
            int nz = nml['grid']['nz']
            int dz = nml['grid']['dz']
        self.path_out = os.path.join(path, 'Updrafts')

        if case_name[0:8] == 'ZGILS_S6':
            self.path_ref = os.path.join(path, 'Stats.ZGILS_S6_1xCO2_SST_FixSub.nc')
        else:
            self.path_ref = os.path.join(path, 'Stats.'+case_name+'.nc')
        self.p_ref = np.zeros([nz], dtype=np.double, order='c')
        try:
            self.p_ref = read_in_netcdf('p0_half', 'reference', self.path_ref)
        except:
            print('no p0_half profile')
            self.p_ref = read_in_netcdf('p0', 'reference', self.path_ref)[:]

        self.z_ref = read_in_netcdf('z', 'reference', self.path_ref)
        self.zrange = krange * dz
        # print('zrange', self.zrange.shape, krange.shape, type(krange), type(self.zrange))

        '''Initialize Latent Heat and ClausiusClapeyron'''

        self.LH = CC_thermodynamics_c.LatentHeat(nml)
        self.CC = CC_thermodynamics_c.ClausiusClapeyron()
        self.CC.initialize(nml, self.LH)


        # (B) Tracer Model
        self.times_tracers = np.asarray([10800, 11700, 12600, 13500, 14400, 15300, 16200, 17100, 18000, 18900, 19800, 20700, 21600], dtype=np.int32)
        self.krange = np.array(krange, dtype=np.int32)


        print('')

        return



    # cpdef update(self, files, ncomp_range, dz_range, krange, nml, path, path_tr):
    cpdef update(self, files, ncomp_range, krange, nml, path, path_tr):
        print('')
        print('----- Updraft Clustering: update -----')

        # (A) read in 3D fields for one time step and one level
        cdef:
            int nvar = 2
            int ncomp = 2
            int nx = nml['grid']['nx']
            int ny = nml['grid']['ny']
            int nz = nml['grid']['nz']
            int dz = nml['grid']['dz']
            double [:,:,:] s_
            double [:,:,:] T_
            double [:,:,:] qt_
            double [:,:,:] ql_

        var_list = ['s', 'qt', 'temperature', 'ql']
        data = np.ndarray(shape=((nx * ny), nvar))
        for d in files:
            if d[-1] == 'c':
                tt = d[0:-3]
            else:
                tt = np.int(d)
            print('- time: d='+str(d) + ', t='+str(tt)+' -')
            nc_file_name = d
            path_fields = os.path.join(path, 'fields', nc_file_name)
            s_, qt_, T_, ql_ = read_in_fields('fields', var_list, path_fields)
            w_ = read_in_netcdf('w', 'fields', path_fields)



            # (B) PDF Model
            # labels_pdf = self.predict_PDF(s_, qt_, T_, ql_, path, ncomp, dz_range, krange, tt, nml)
            labels_pdf, clf = self.predict_PDF(s_, qt_, T_, ql_, path, ncomp, krange, tt, nml)
            up_type = 'PDF'
            plot_labels(qt_, ql_, w_, labels_pdf, tt, krange, dz, up_type, path)

            # (C) Tracer Model
            type_list = ['Couvreux', 'Cloud', 'Coherent', 'Core']
            for up_type in type_list:
                labels_tr = self.read_in_updrafts_colleen(up_type, tt, path_tr)
                plot_labels(qt_, ql_, w_, labels_tr, tt, krange, dz, up_type, path)
                # plot_labels_test(qt_, ql_, w_, labels_tr, tt, krange, dz, up_type, path)
                plot_labels_comparison(qt_, ql_, labels_pdf, labels_tr, up_type, tt, krange, dz, path)


        return





    # ----------------------------------------------------------------------
    #               PDF Model
    # ----------------------------------------------------------------------
    # cpdef predict_PDF(self, s_, qt_, T_, ql_, path, int ncomp, dz_range_, krange_, tt, nml):
    cpdef predict_PDF(self, s_, qt_, T_, ql_, path, int ncomp, krange_, tt, nml):
        # input:
        #     - 3D fields at given time tt
        # output:
        #     - labels = (nx,ny,nz):      3D array with labels from PDF clustering
        #
        # (A) Set up files
        # (B) Compute PDF
        #       - read in fields (and accumulate over several time steps)
        #       for each k...
        #           - normalise data
        #           - compute PDF for two components
        #           - compute labels
        #           - sort labels according to <qt>
        #           - rearrange labels into 3D array
        # (C) Compute Error
        # (D) IO:
        #       - output updrafts file (labels)
        #       - output PDF parameters (sorted)
        #       - output errors


        print('')
        print('--- PDF Prediction ---')
        print('')
        cdef extern from "thermodynamic_functions.h":
            inline double thetali_c(const double p0, const double T, const double qt, const double ql, const double qi, const double L)

        # ________________________________________________________________________________________
        '''(A) Set up files etc.'''
        cdef:
            # int ncomp
            int nvar = 2
            int [:] krange = krange_
            int nk = len(krange)
            # int dz_range = dz_range_
            Py_ssize_t k, iz
            str d
            int dz = nml['grid']['dz']
            int nx = nml['grid']['nx']
            int ny = nml['grid']['ny']
            int nz = nml['grid']['nz']
            double [:] p_ref = np.zeros([nz],dtype=np.double,order='c')       # <type 'CloudClosure._memoryviewslice'>

        for k in range(nz):
            p_ref[k] = self.p_ref[k]

        # " (a) from Profiles "
        # time_profile = read_in_netcdf('t', 'timeseries', self.path_ref)
        # nt = time_profile.shape[0]
        # print('time in stats file: ', time_profile.shape, nt, dt_stat)
        # ql_profile = read_in_netcdf('ql_mean', 'profiles', path_ref)

        # N = len(files)

        cdef:
            # LatentHeat LH = CC_thermodynamics_c.LatentHeat(nml)
            # ClausiusClapeyron CC = CC_thermodynamics_c.ClausiusClapeyron()
            LatentHeat LH = self.LH
            ClausiusClapeyron CC = self.CC
        # CC.initialize(nml, LH)
        print('')


        # ________________________________________________________________________________________
        '''(B) Compute PDF f(s,qt) from LES data'''
        #       - read in fields
        #       - compute theta_l
        #       - compute PDFs f(s,qt), g(th_l, qt)

        '''(1) Read in Fields'''
        # for PDF construction
        cdef:
            int i, j, ij
            int ishift = ny
            double [:,:] qt = np.zeros([nx*ny,nk],dtype=np.double,order='c')         # <type 'CloudClosure._memoryviewslice'>
            double [:,:] ql = np.zeros([nx*ny,nk],dtype=np.double,order='c')         # <type 'CloudClosure._memoryviewslice'>
            double qi_ = 0.0
            double [:,:] theta_l = np.zeros([nx*ny,nk],dtype=np.double,order='c')         # <type 'CloudClosure._memoryviewslice'>
            double Lv

        # for PDF sampling
        cdef:
            int n_sample = np.int(1e6)
            double [:] T_comp_thl = np.zeros([n_sample],dtype=np.double,order='c')
            double [:] ql_comp_thl = np.zeros([n_sample],dtype=np.double,order='c')

        # for Error Computation
        cdef:
            double [:] ql_mean_field = np.zeros(nk, dtype=np.double)        # computation from 3D LES field
            double [:] ql_mean_comp = np.zeros(nk, dtype=np.double)         # computation from PDF sampling
            double [:] cf_field = np.zeros(shape=(nk))                      # computation from 3D LES field
            double [:] cf_comp = np.zeros(nk, dtype=np.double)              # computation from PDF sampling
            double [:] error_ql = np.zeros(shape=(nk))
            double [:] rel_error_ql = np.zeros(shape=(nk))
            double [:] error_cf = np.zeros(shape=(nk))
            double [:] rel_error_cf = np.zeros(shape=(nk))

        # for Labeling
        # cdef:
        #     double [:,:,:] labels = np.zeros(shape=(nx, ny, nk))
        labels = np.zeros(shape=(nx, ny, nk))

        data = np.ndarray(shape=((nx * ny), nvar))
        means_ = np.ndarray(shape=(nk, ncomp, nvar))
        covariances_ = np.zeros(shape=(nk, ncomp, nvar, nvar))
        weights_ = np.zeros(shape=(nk, ncomp))

        '''(1) Statistics File'''
        nc_file_name_CC = 'CC_updrafts_time'+str(tt)+'.nc'
        self.create_statistics_file(self.path_out, nc_file_name_CC, tt, ncomp, nvar, nk)
        nc_file_name_labels = 'Labeling_t'+str(tt)+'.nc'
        self.create_updrafts_file(self.path_out, nc_file_name_labels, tt, ncomp, nvar, nk, nml)

        for k in range(nk):
            iz = krange[k]
            print('- z = '+str(iz*dz)+ ', ncomp = '+str(ncomp)+' -')
            '''(2) Compute liquid potential temperature from temperature and moisture'''
            for i in range(nx):
                for j in range(ny):
                    ij = i*ishift + j
                    Lv = LH.L(T_[i,j,iz],LH.Lambda_fp(T_[i,j,iz]))
                    theta_l[ij,k] = thetali_c(p_ref[iz], T_[i,j,iz], qt_[i,j,iz], ql_[i,j,iz], qi_, Lv)
                    qt[ij,k] = qt_[i,j,iz]
                    ql[ij,k] = ql_[i,j,iz]
                    ql_mean_field[k] += ql_[i,j,iz]
                    if ql[ij,k] > 0.0:
                        cf_field[k] += 1.0
            data[:, 0] = theta_l[:, k]
            data[:, 1] = qt[:, k]

            ql_mean_field[k] /= (nx*ny)
            cf_field[k] /= (nx*ny)

            '''(3) Normalise Data'''
            scaler = StandardScaler()
            data_norm = scaler.fit_transform(data)

            '''(4) Compute bivariate PDF'''
            #   (a) for (s,qt)
            #   (b) for (th_l,qt)
            clf_thl_norm = mixture.GaussianMixture(n_components=ncomp,covariance_type='full')
            clf_thl_norm.fit(data_norm)

            '''(5) Find Labels, sort and rearrange'''
            # sort PDF components, s.t. 0/1 always correspond to the same component (environment vs. updrafts)
            # Note: it is either (# of zeros in labels_new) = (# of zeros in labels_) or (# of zeros in labels_new) = (nx*ny - (# of zeros in labels_))
            labels_ = clf_thl_norm.predict(data_norm)

            means_[k, :, :], covariances_[k, :, :, :], weights_[k, :], labels_new = self.sort_PDF(clf_thl_norm.means_,
                                                              clf_thl_norm.covariances_, clf_thl_norm.weights_, labels_)
            print('Labels: ', np.count_nonzero(labels_), np.count_nonzero(labels_new), labels_.shape[0]-np.count_nonzero(labels_))
            if ( np.count_nonzero(labels_new) != labels_.shape[0]-np.count_nonzero(labels_)
                 and np.count_nonzero(labels_new) != np.count_nonzero(labels_) ):
                print('!!!!! Labels Problem !!!!')
            print('')
            # means_[k, :, :] = clf_thl_norm.means_
            # covariances_[k, :,:,:] = clf_thl_norm.covariances_
            # weights_[k, :] = clf_thl_norm.weights_
            # labels_new = labels_

            # rearrange into 2D array (for only one data file)
            for i in range(nx):
                for j in range(ny):
                    ij = i*ishift + j
                    labels[i,j,k] = labels_new[ij]
                    # '''(D) Compute mean liquid water <ql> from PDF f(s,qt)'''

            '''(6) Compute Error in <ql> and CF'''
            #   a. sample (th_l, qt) from PDF (Monte Carlo ???
            #   b. compute ql for samples
            #   c. consider ensemble average <ql> = domain mean representation???
            '''     (a) Draw samples'''
            Th_l_norm, y_norm = clf_thl_norm.sample(n_samples=n_sample)
            '''     (b) Rescale theta_l and qt'''
            Th_l = scaler.inverse_transform(Th_l_norm)      # Inverse Normalisation

            '''     (c) Compute ql (saturation adjustment) & Cloud Fraction '''
            for i in range(n_sample-2):
                T_comp_thl[i], ql_comp_thl[i] = sat_adj_fromthetali(p_ref[iz], Th_l[i, 0], Th_l[i, 1], CC, LH)
                ql_mean_comp[k] = ql_mean_comp[k] + ql_comp_thl[i]
                if ql_comp_thl[i] > 0:
                    cf_comp[k] += 1
            ql_mean_comp[k] = ql_mean_comp[k] / n_sample
            cf_comp[k] = cf_comp[k] / n_sample
            error_ql[k] = ql_mean_comp[k] - ql_mean_field[k]
            error_cf[k] = cf_comp[k] - cf_field[k]
            if ql_mean_field[k] > 0.0:
                rel_error_ql[k] = (ql_mean_comp[k] - ql_mean_field[k]) / ql_mean_field[k]
            if cf_field[k] > 0.0:
                rel_error_cf[k] = (cf_comp[k] - cf_field[k]) / cf_field[k]

            print('')
            print('<ql> from CloudClosure Scheme: ', ql_mean_comp[k])
            print('<ql> from ql fields:           ', ql_mean_field[k])
            print('error (<ql>_CC - <ql>_field): '+str(error_ql[k]))
            print('rel err: '+ str(rel_error_ql[k]))
            print('')
            print('CF from Cloud Closure Scheme: ', cf_comp[k])
            print('CF from ql fields: ', cf_field[k])
            print('error: '+str(error_cf[k]))
            print('rel error: ', rel_error_cf[k])
            print('')


        plot_parameters_pdf(labels, means_, covariances_, weights_, krange_, self.zrange, tt, dz, path)

        '''(6) Save Updraft Labels & Gaussian Mixture PDFs '''
        self.write_updrafts_file(self.path_out, nc_file_name_labels, labels, nml)

        print('')
        print('Dumping files: '+ self.path_out)
        self.dump_variable(os.path.join(self.path_out, nc_file_name_CC), 'means', means_, 'qtT', ncomp, nvar, nk)
        self.dump_variable(os.path.join(self.path_out, nc_file_name_CC), 'covariances', covariances_, 'qtT', ncomp, nvar, nk)
        self.dump_variable(os.path.join(self.path_out, nc_file_name_CC), 'weights', weights_, 'qtT', ncomp, nvar, nk)
        self.dump_variable(os.path.join(self.path_out, nc_file_name_CC), 'error', np.asarray(error_ql[:]), 'error_ql', ncomp, nvar, nk)
        self.dump_variable(os.path.join(self.path_out, nc_file_name_CC), 'error', np.asarray(rel_error_ql[:]), 'rel_error_ql', ncomp, nvar, nk)
        self.dump_variable(os.path.join(self.path_out, nc_file_name_CC), 'error', np.asarray(error_cf[:]), 'error_cf', ncomp, nvar, nk)
        self.dump_variable(os.path.join(self.path_out, nc_file_name_CC), 'error', np.asarray(rel_error_cf[:]), 'rel_error_cf', ncomp, nvar, nk)

        return labels, clf_thl_norm



    cpdef sort_PDF(self, means_, covariances_, weights_, labels_):
        '''sort PDF components according to <qt>'''
        # sort such that PDF-component 1 has larger mean qt than PDF-component 0, i.e. PDF-component 1 represents the updrafts

        # means_ = (ncomp x nvar)
        # covariances_ = (ncomp, nvar, nvar)
        # weights_ = (ncomp)
        # labels_ = (nx * ny)
        cdef:
            double [:,:] means = means_
            double [:,:,:] covars = covariances_
            double [:] weights = weights_
            # double [:] labels = np.ndarray(shape=labels_.shape, ndtype=np.double)
            # double [:] labels = np.ndarray(shape=labels_.shape)
            int nk = len(self.zrange)
            int nij = labels_.shape[0]
            int nvar = 2

        labels = np.zeros(shape = labels_.shape)

        # change PDF-component label if mean qt of component 0 smaller than mean qt of PDF-component 1
        if means[0, 1] < means[1, 1]:
            # print('')
            print('sorting')
            aux = weights[1]
            weights[1] = weights[0]
            weights[0] = aux
            for i1 in range(nvar):  # loop over variables
                aux = means[1, i1]
                means[1, i1] = means[0, i1]
                means[0, i1] = aux
                for i2 in range(nvar):
                    aux = covars[1, i1, i2]
                    covars[1, i1, i2] = covars[0, i1, i2]
                    covars[0, i1, i2] = aux
            # labels = [int(not labels_[i]) for i in range(nij)]
            labels = labels_
        else:
            # labels = labels_
            labels = [int(not labels_[i]) for i in range(nij)]
            # # trivar_plot_means(var, weights, means, covars, mean_tot, covars_tot, time, 'sortB')
            # # trivar_plot_covars(var, weights, means, covars, mean_tot, covars_tot, time, 'sortB')
            # # trivar_plot_weights(var, weights, means, covars, mean_tot, covars_tot, time, 'sortB')
            # # print('')
        return means, covars, weights, labels




    # cpdef sort_PDF_allk(self, means_, covariances_, weights_, labels_):
    #     '''sort PDF components according to <qt>'''
    #     # sort such that PDF-component 0 has larger mean qt than PDF-component 1, i.e. PDF-component 0 represents the updrafts
    #
    #     # means_ = (nk x ncomp x nvar)
    #     # covariances_ = (nk x ncomp x nvar x nvar)
    #     # weights_ = (nk x ncomp)
    #     # labels_ = (nx x ny x nk)
    #
    #     cdef:
    #         double [:,:,:] means = means_[:,:,:]
    #         double [:,:,:,:] covars = covariances_
    #         double [:,:] weights = weights_
    #         double [:,:,:] labels = np.ndarray(shape=labels_.shape,ndtype=np.double)
    #         int nk = len(self.zrange)
    #         # Py_ssize_t nij = labels.shape
    #         int nij = labels_.shape
    #         int nvar = 2
    #     # means = means_
    #     # covars = covariances_
    #     # weights = weights_
    #
    #     print('')
    #     print('sorting', nij)
    #     print(means.shape, covars.shape, weights.shape)
    #     for k in range(nk):
    #         # change PDF-component label if mean qt of component 0 smaller than mean qt of PDF-component 1
    #         if means[k, 0, 1] < means[k, 1, 1]:
    #             aux = weights[k, 1]
    #             weights[k, 1] = weights[k, 0]
    #             weights[k, 0] = aux
    #             for i1 in range(nvar):  # loop over variables
    #                 aux = means[k, 1, i1]
    #                 means[k, 1, i1] = means[k, 0, i1]
    #                 means[k, 0, i1] = aux
    #                 for i2 in range(nvar):
    #                     aux = covars[k, 1, i1, i2]
    #                     covars[k, 1, i1, i2] = covars[k, 0, i1, i2]
    #                     covars[k, 0, i1, i2] = aux
    #     #         # for i in range(nx):
    #     #         #     for j in range(ny):
    #     #         a = [int(not labels[i]) for i in range(nij)]
    #     #
    #     # print(a[0:3])
    #     # print(labels[0:3])
    #
    #     # trivar_plot_means(var, weights, means, covars, mean_tot, covars_tot, time, 'sortB')
    #     # trivar_plot_covars(var, weights, means, covars, mean_tot, covars_tot, time, 'sortB')
    #     # trivar_plot_weights(var, weights, means, covars, mean_tot, covars_tot, time, 'sortB')
    #     # print('')
    #     return means, covars, weights, labels


    #----------------------------------------------------------------------
    #----------------------------------------------------------------------
    cpdef create_updrafts_file(self, path, file_name, time, ncomp, nvar, nz_, nml):
        print('create updrafts file: '+ path +', '+ file_name)
        # ncomp: number of Gaussian components in EM
        # nvar: number of variables of multi-variate Gaussian components
        rootgrp = nc.Dataset(os.path.join(path, file_name), 'w', format='NETCDF4')
        # dimgrp = rootgrp.createGroup('dims')

        z_grp = rootgrp.createGroup('profiles')
        nk = len(self.zrange)
        z_grp.createDimension('nz', nk)
        var = z_grp.createVariable('z', 'f8', 'nz')
        for k in range(nk):
            var[k] = self.zrange[k]
        # var[:] = self.zrange[:]

        field_grp = rootgrp.createGroup('fields')
        field_grp.createDimension('nx', nml['grid']['nx'])
        field_grp.createDimension('ny', nml['grid']['ny'])
        # field_grp.createDimension('nz', nml['grid']['nz'])
        field_grp.createDimension('nz', nk)
        labels = field_grp.createVariable('labels', 'f8', ('nx','ny','nz'))
        labels.units = ' '
        rootgrp.close()
        return



    cpdef write_updrafts_file(self, path, file_name, data, nml):
        print('write updrafts: ' + os.path.join(path, file_name))
        root = nc.Dataset(os.path.join(path, file_name), 'r+', format='NETCDF4')
        nk = len(self.zrange)
        labels = root.groups['fields'].variables['labels']
        labels[:,:,:] = data
        root.close()
        # # field_grp = root.groups['fields']
        # # var = field_grp.variables['labels']
        # # var = root.groups['fields'].variables['labels']
        # # var = fieldgrp.variables[name]
        # # var[:,:,:] = np.array(data)
        return

    #----------------------------------------------------------------------
    cpdef create_statistics_file(self, path, file_name, time, ncomp, nvar, nz_):
        cdef:
            int i
            int nz = nz_

        print('create statistics file: '+ path+', '+ file_name)
        # ncomp: number of Gaussian components in EM
        # nvar: number of variables of multi-variate Gaussian components
        rootgrp = nc.Dataset(os.path.join(path,file_name), 'w', format='NETCDF4')
        dimgrp = rootgrp.createGroup('dims')

        ts_grp = rootgrp.createGroup('time')
        ts_grp.createDimension('nt',len(time))
        var = ts_grp.createVariable('t','f8',('nt'))
        for i in range(len(time)-1):
            var[i] = time[i+1]
        z_grp = rootgrp.createGroup('z-profile')
        z_grp.createDimension('nz', nz)
        var = z_grp.createVariable('height', 'f8', ('nz'))
        for i in range(nz):
            var[i] = self.zrange[i]

        means_grp = rootgrp.createGroup('means')
        means_grp.createDimension('nz', nz)
        means_grp.createDimension('ncomp', ncomp)
        means_grp.createDimension('nvar', nvar)
        cov_grp = rootgrp.createGroup('covariances')
        cov_grp.createDimension('nz', nz)
        cov_grp.createDimension('ncomp', ncomp)
        cov_grp.createDimension('nvar', nvar)
        weights_grp = rootgrp.createGroup('weights')
        weights_grp.createDimension('nz', nz)
        weights_grp.createDimension('ncomp', ncomp)

        error_grp = rootgrp.createGroup('error')
        error_grp.createDimension('nz', nz)
        rootgrp.close()
        return



    cpdef dump_variable(self, path, group_name, data_, var_name, ncomp, nvar, nz_):
        print('-------- dump variable --------', var_name, group_name, path)
        # # print('dump variable', path, group_name, var_name, data_.shape, ncomp, nvar)
        rootgrp = nc.Dataset(path, 'r+')
        if group_name == 'means':
            # rootgrp = nc.Dataset(path, 'r+')
            var = rootgrp.groups['means'].createVariable(var_name, 'f8', ('nz', 'ncomp', 'nvar'))
            # var = nc.Dataset(path, 'r+').groups['means'].createVariable(var_name, 'f8', ('nz', 'ncomp', 'nvar'))
            # var = nc.Dataset(path, 'r+').groups['means'].createVariable(var_name, 'f8', ('nz', 'ncomp', 'nvar'))[:,:,:]
            var[:,:,:] = data_[:,:,:]

        elif group_name == 'covariances':
            var = rootgrp.groups['covariances'].createVariable(var_name, 'f8', ('nz', 'ncomp', 'nvar', 'nvar'))
            var[:,:,:,:] = data_[:,:,:,:]

        elif group_name == 'weights':
            var = rootgrp.groups['weights'].createVariable(var_name, 'f8', ('nz', 'ncomp'))
            var[:,:] = data_[:,:]

        elif group_name == 'error':
            var = rootgrp.groups['error'].createVariable(var_name, 'f8', ('nz'))
            var[:] = data_[:]

        # # write_field(path, group_name, data, var_name)
        # # print('--------')
        rootgrp.close()
        print('')
        return


    # ----------------------------------------------------------------------
    #                Tracer Model
    # ----------------------------------------------------------------------
    cpdef read_in_updrafts_colleen(self, str type, t, path_):
        print('')
        print('--- Updraft Colleen: read in ---')
        import pickle
        print('')

        # path = os.path.join(path_, 'updrafts_colleen')
        path = path_
        files = os.listdir(path)
        # times = ['10800', '11700', '12600', '13500', '14400', '15300', '16200', '17100', '18000', '18900', '19800', '20700', '21600']
        # times = ['10800']
        print(path_)
        print('time: ', t)
        # print(files)
        print('')


        if type == 'Cloud':
            root = 'Bomex_Cloud_'
        elif type == 'Coherent':
            root = 'Bomex_Coherent_'
        elif type == 'Couvreux':
            root = 'Bomex_Couvreux_'
        elif type == 'Core':
            root = 'Bomex_Core_'
        print(root)



        print('')
        path = os.path.join(path_, root + 'updraft.pkl')
        data = pickle.load(open(path))
        print(path +': ', data.keys())
        # print(data.items())
        # dict.items(): Return a copy of the dictionary’s list of (key, value) pairs.

        print('')
        path = os.path.join(path_, root + 'environment.pkl')
        data = pickle.load(open(path))
        print(path +': ', data.keys())


        print('')
        path = os.path.join(path_, root + 'time_' + t + '_Grid.pkl')
        labels = pickle.load(open(path))
        # data = pickle.load(open(path))
        print(path +': ', labels.shape)


        # # (a)
        # root1 = 'Bomex_Cloud_'
        # print('')
        # for t in times:
        #     path = os.path.join(path_, root1 + 'time_' + t + '_Grid.pkl')
        #     # data = pickle.load(open(path))
        #     print(path)
        # path = os.path.join(path_, root1 + 'updraft.pkl')
        # path = os.path.join(path_, root1 + 'environment.pkl')
        #
        # # (b)
        # root2 = 'Bomex_Coherent_'
        # print('')
        # for t in times:
        #     path = os.path.join(path_, root2 + 'time_' + t + '_Grid.pkl')
        #     # data = pickle.load(open(path))
        #     print(path)
        # path = os.path.join(path_, root2 + 'updraft.pkl')
        # path = os.path.join(path_, root2 + 'environment.pkl')
        #
        # # (c)
        # root3 = 'Bomex_Couvreux_'
        # print('')
        # for t in times:
        #     path = os.path.join(path_, root3 + 'time_' + t + '_Grid.pkl')
        #     # data = pickle.load(open(path))
        #     print(path)
        # path = os.path.join(path_, root3 + 'updraft.pkl')
        # path = os.path.join(path_, root3 + 'environment.pkl')
        #
        # # (d)
        # root4 = 'Bomex_Core_'
        # print('')
        # for t in times:
        #     path = os.path.join(path_, root4 + 'time_' + t + '_Grid.pkl')
        #     # data = pickle.load(open(path))
        #     print(path)

        return labels


    # ----------------------------------------------------------------------
    #                Comparison
    # ----------------------------------------------------------------------
    # cpdef read_in_data()




# ----------------------------------------------------------------------
def read_in_netcdf(variable_name, group_name, fullpath_in):
    print('io_read_in_files: read in netcdf', variable_name, group_name, fullpath_in)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    # grp = rootgrp.groups[group_name]
    if group_name == 'reference' or group_name == 'timeseries':
        var = rootgrp.groups[group_name].variables[variable_name][:]
        rootgrp.close()
    elif group_name == 'profiles':
        var = rootgrp.groups[group_name].variables[variable_name][:,:]
        rootgrp.close()
    elif group_name == 'fields':
        var = rootgrp.groups[group_name].variables[variable_name][:,:,:]
        rootgrp.close()
    else:
        var_ = rootgrp.groups[group_name].variables[variable_name]
        shape_ = var_.shape
        dims_ = np.size(shape_)
        var = np.ndarray(shape=shape_)
        var[:] = var_[:]
    return var


def read_in_fields(group_name, var_list, path):
    rootgrp = nc.Dataset(path, 'r')
    grp = rootgrp.groups[group_name]
    if group_name == 'reference':
        var1 = grp.variables[var_list[0]][:]
        var2 = grp.variables[var_list[1]][:]
        var3 = grp.variables[var_list[2]][:]
        var4 = grp.variables[var_list[3]][:]
        rootgrp.close()
    elif group_name == 'fields':
        var1 = grp.variables[var_list[0]][:,:,:]
        var2 = grp.variables[var_list[1]][:,:,:]
        var3 = grp.variables[var_list[2]][:,:,:]
        var4 = grp.variables[var_list[3]][:,:,:]
        rootgrp.close()
    return var1, var2, var3, var4