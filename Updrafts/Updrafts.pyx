# (1)   a) read in 3D output files from LES
#
#       b) compute double-Gaussian PDF from 3D output file per level
#           - for given level make 1D array out of 2D array
#           -
# (2) find label (to which cluster each data point belongs) for each data point
#           Attribute: clf.predict(X[, y]), X.shape = (# data, # features)
#           - returns array with 0/1 of size #data
#           --> return into 2D array
#           --> sort PDF components, s.t. 0/1 always correspond to the same component (environment vs. updrafts)
# (3) read in Colleens output data (HDF5-data with 0/1 encoding of updraft or environment)
# (4) compare both clustering algorithms:
#       - plotting: coloring scheme
#           blue = environment in both
#           green = updraft in both
#           red = differences (could be further distinguised into which algorithm classifies it as an updraft point)
#       - compute conditional statistics, e.g. cloud fraction and compare


# Option: read in fields and PDF parameters from Cloud Closure and use clf.predict(field_data) for producing labels


import os
import json as simplejson
import time
import netCDF4 as nc
import numpy as np
from sklearn import mixture
from sklearn.preprocessing import StandardScaler


import CC_thermodynamics_c
from CC_thermodynamics_c cimport LatentHeat, ClausiusClapeyron
from CC_thermodynamics_c import sat_adj_fromentropy, sat_adj_fromthetali

cdef class CloudClosure:
    def __init__(self):
        return
    cpdef initialize(self):
        return


cdef class Updrafts:
    def __init__(self):
        self.path_ref = ' '
        self.p_ref = None
        self.z_ref = None
        return

    cpdef initialize(self, krange, path, case_name):
        print('--- Updraft Clustering ---')
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

        return




    cpdef predict_pdf(self, files, path, ncomp_range, dz_range_, krange_, nml):
        print('')
        print('--- PDF Prediction ---')
        print('')
        cdef extern from "thermodynamic_functions.h":
            inline double thetali_c(const double p0, const double T, const double qt, const double ql, const double qi, const double L)

        time1 = time.clock()

        # ________________________________________________________________________________________
        '''(A) Set up files etc.'''
        cdef:
            int ncomp
            int nvar = 2
            int [:] krange = krange_
            int nk = len(krange)
            int dz_range = dz_range_
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

        '''(B) Initialize Latent Heat and ClausiusClapeyron'''
        cdef:
            LatentHeat LH = CC_thermodynamics_c.LatentHeat(nml)
            ClausiusClapeyron CC = CC_thermodynamics_c.ClausiusClapeyron()
        CC.initialize(nml, LH)
        print('')


        # ________________________________________________________________________________________
        '''(C) Compute PDF f(s,qt) from LES data'''
        #       - read in fields
        #       - compute theta_l
        #       - compute PDFs f(s,qt), g(th_l, qt)

        '''(1) Read in Fields'''
        # for PDF construction
        cdef:
            double [:,:,:] s_
            double [:,:,:] T_
            double [:,:,:] qt_
            double [:,:] qt = np.zeros([nx*ny,nk],dtype=np.double,order='c')         # <type 'CloudClosure._memoryviewslice'>
            double [:,:,:] ql_
            double [:,:] ql = np.zeros([nx*ny,nk],dtype=np.double,order='c')         # <type 'CloudClosure._memoryviewslice'>
            double [:] ql_all = np.ndarray(shape=(0))
            double qi_ = 0.0
            double [:,:] theta_l = np.zeros([nx*ny,nk],dtype=np.double,order='c')         # <type 'CloudClosure._memoryviewslice'>
            double [:,:] data_all
            int ishift = ny
            int i, j, ij
            double Lv

        # for PDF sampling
        cdef:
            int n_sample = np.int(1e6)
            double [:] cf_field = np.zeros(shape=(nk))
            double [:] T_comp_thl = np.zeros([n_sample],dtype=np.double,order='c')
            double [:] ql_comp_thl = np.zeros([n_sample],dtype=np.double,order='c')
            double [:] alpha_comp_thl = np.zeros(n_sample)

        # for Error Computation
        cdef:
            double ql_mean_field
            double ql_mean = 0.0
            double cf_comp = 0.0
            int count_ncomp = 0
            double [:,:] error_ql = np.zeros(shape=(nk,len(ncomp_range)))
            double [:,:] rel_error_ql = np.zeros(shape=(nk,len(ncomp_range)))
            double [:,:] error_cf = np.zeros(shape=(nk,len(ncomp_range)))
            double [:,:] rel_error_cf = np.zeros(shape=(nk,len(ncomp_range)))

        # for Labeling
        # cdef:
        #     double [:,:,:] labels = np.zeros(shape=(nx, ny, nk))
        labels = np.zeros(shape=(nx, ny, nk))


        var_list = ['s', 'qt', 'temperature', 'ql']
        data = np.ndarray(shape=((nx * ny), nvar))
        for ncomp in ncomp_range:
            means_ = np.ndarray(shape=(nk, ncomp, nvar))
            covariance_ = np.zeros(shape=(nk, ncomp, nvar, nvar))
            weights_ = np.zeros(shape=(nk, ncomp))
            '''(1) Statistics File'''
            tt = files[0][0:-3]
            nc_file_name_CC = 'CC_updrafts_time'+str(tt)+'.nc'
            self.create_statistics_file(self.path_out, nc_file_name_CC, tt, ncomp, nvar, nk)
            nc_file_name_labels = 'Labeling_time'+str(tt)+'.nc'
            self.create_updrafts_file(self.path_out, nc_file_name_labels, tt, ncomp, nvar, nk, nml)

            for k in range(nk):
                iz = krange[k]
                print('- z = '+str(iz*dz)+ ', ncomp = '+str(ncomp)+' -')
                data_all = np.ndarray(shape=(0, nvar))
                ql_all = np.ndarray(shape=(0))
                ql_mean_field = 0.0

                for d in files:
                    nc_file_name = d
                    path_fields = os.path.join(path, 'fields', nc_file_name)
                    # print('fields: ', path_fields)
                    s_, qt_, T_, ql_ = read_in_fields('fields', var_list, path_fields)
                    '''(2) Compute liquid potential temperature from temperature and moisture'''
                    for i in range(nx):
                        for j in range(ny):
                            ij = i*ishift + j
                            Lv = LH.L(T_[i,j,iz],LH.Lambda_fp(T_[i,j,iz]))
                            theta_l[ij,k] = thetali_c(p_ref[iz], T_[i,j,iz], qt_[i,j,iz], ql_[i,j,iz], qi_, Lv)
                            qt[ij,k] = qt_[i,j,iz]
                            ql[ij,k] = ql_[i,j,iz]
                            ql_mean_field += ql_[i,j,iz]
                            if ql[ij,k] > 0.0:
                                cf_field[k] += 1.0
                    del s_, ql_, qt_, T_
                    data[:, 0] = theta_l[:, k]
                    data[:, 1] = qt[:, k]
                    data_all = np.append(data_all, data, axis=0)
                    ql_all = np.append(ql_all, ql[:,k], axis=0)

                ql_mean_field /= len(files)*(nx*ny)
                cf_field[k] /= len(files)*(nx*ny)

                '''(3) Normalise Data'''
                scaler = StandardScaler()
                data_all_norm = scaler.fit_transform(data_all)

                '''(4) Compute bivariate PDF'''
                #   (a) for (s,qt)
                #           ...
                #   (b) for (th_l,qt)
                # clf_thl_norm = Gaussian_bivariate(ncomp, data_all_norm, 'T', 'qt', np.int(d[0:-3]), iz * dz, path)
                clf_thl_norm = mixture.GaussianMixture(n_components=ncomp,covariance_type='full')
                clf_thl_norm.fit(data_all_norm)
                means_[k, :, :] = clf_thl_norm.means_[:, :]
                covariance_[k,:,:,:] = clf_thl_norm.covariances_[:,:,:]
                weights_[k,:] = clf_thl_norm.weights_[:]


                '''(5) Sort and look for labels'''
                labels_ = clf_thl_norm.predict(data_all_norm)
                # labels_ = np.int(clf_thl_norm.predict(data_all_norm))
                # print('LABELS: ', type(labels))
                print('')
                print('Labels: ', np.count_nonzero(labels_), labels_.shape, data_all_norm.shape)
                # sort PDF components, s.t. 0/1 always correspond to the same component (environment vs. updrafts)
                means, covars, weight, labels_new = self.sort_PDF(clf_thl_norm.means_, clf_thl_norm.covariances_, clf_thl_norm.weights_, labels_)
                print('Labels new: ', np.count_nonzero(labels_new))

                # rearrange into 2D array (for only one data file)
                for i in range(nx):
                    for j in range(ny):
                        ij = i*ishift + j
                        labels[i,j,k] = labels_new[ij]


            #     '''(D) Compute mean liquid water <ql> from PDF f(s,qt)'''
            #     #       1. sample (th_l, qt) from PDF (Monte Carlo ???
            #     #       2. compute ql for samples
            #     #       3. consider ensemble average <ql> = domain mean representation???
            #     '''(1) Draw samples'''
            #     Th_l_norm, y_norm = clf_thl_norm.sample(n_samples=n_sample)
            #     '''(2) Rescale theta_l and qt'''
            #     Th_l = scaler.inverse_transform(Th_l_norm)      # Inverse Normalisation
            #
            #     '''(3) Compute ql (saturation adjustment) & Cloud Fraction '''
            #     for i in range(n_sample-2):
            #         # T_comp[i], ql_comp[i], alpha_comp[i] = sat_adj_fromentropy(p_ref[iz-1], S[i,0],S[i,1])
            #         T_comp_thl[i], ql_comp_thl[i], alpha_comp_thl[i] = sat_adj_fromthetali(p_ref[iz], Th_l[i, 0], Th_l[i, 1], CC, LH)
            #         ql_mean = ql_mean + ql_comp_thl[i]
            #         if ql_comp_thl[i] > 0:
            #             cf_comp += 1
            #     ql_mean = ql_mean / n_sample
            #     cf_comp = cf_comp / n_sample
            #
            #     error_ql[k,count_ncomp] = ql_mean - ql_mean_field
            #     error_cf[k,count_ncomp] = cf_comp - cf_field[k]
            #     if ql_mean_field > 0.0:
            #         rel_error_ql[k,count_ncomp] = (ql_mean - ql_mean_field) / ql_mean_field
            #     if cf_field[k] > 0.0:
            #         rel_error_cf[k,count_ncomp] = (cf_comp - cf_field[k]) / cf_field[k]
            #
            #     print('<ql> from CloudClosure Scheme: ', ql_mean)
            #     print('<ql> from ql fields: ', ql_mean_field)
            #     print('error (<ql>_CC - <ql>_field): '+str(error_ql[k,count_ncomp]))
            #     print('rel err: '+ str(rel_error_ql[k,count_ncomp]))
            #
            # count_ncomp += 1

            # sort PDF components, s.t. 0/1 always correspond to the same component (environment vs. updrafts)
            # means, covars, weight, labels = self.sort_PDF(means_, covariance_, weights_, labels_)

            # self.write_updrafts_field(self.path_out, nc_file_name_labels, labels)
            self.write_updrafts_field(self.path_out, nc_file_name_labels, labels, nml)

        return



    cpdef sort_PDF(self, means_, covariances_, weights_, labels_):
        '''sort PDF components according to <qt>'''
        # sort such that PDF-component 0 has larger mean qt than PDF-component 1, i.e. PDF-component 0 represents the updrafts

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

        print(means.shape, covars.shape, weights.shape)
        # change PDF-component label if mean qt of component 0 smaller than mean qt of PDF-component 1
        if means[0, 1] < means[1, 1]:
            print('')
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
            labels = [int(not labels_[i]) for i in range(nij)]

            # # trivar_plot_means(var, weights, means, covars, mean_tot, covars_tot, time, 'sortB')
            # # trivar_plot_covars(var, weights, means, covars, mean_tot, covars_tot, time, 'sortB')
            # # trivar_plot_weights(var, weights, means, covars, mean_tot, covars_tot, time, 'sortB')
            # # print('')
        return means, covars, weights, labels




    cpdef sort_PDF_allk(self, means_, covariances_, weights_, labels_):
        '''sort PDF components according to <qt>'''
        # sort such that PDF-component 0 has larger mean qt than PDF-component 1, i.e. PDF-component 0 represents the updrafts

        # means_ = (nk x ncomp x nvar)
        # covariances_ = (nk x ncomp x nvar x nvar)
        # weights_ = (nk x ncomp)
        # labels_ = (nx x ny x nk)

        cdef:
            double [:,:,:] means = means_[:,:,:]
            double [:,:,:,:] covars = covariances_
            double [:,:] weights = weights_
            double [:,:,:] labels = np.ndarray(shape=labels_.shape,ndtype=np.double)
            int nk = len(self.zrange)
            # Py_ssize_t nij = labels.shape
            int nij = labels_.shape
            int nvar = 2
        # means = means_
        # covars = covariances_
        # weights = weights_

        print('')
        print('sorting', nij)
        print(means.shape, covars.shape, weights.shape)
        for k in range(nk):
            # change PDF-component label if mean qt of component 0 smaller than mean qt of PDF-component 1
            if means[k, 0, 1] < means[k, 1, 1]:
                aux = weights[k, 1]
                weights[k, 1] = weights[k, 0]
                weights[k, 0] = aux
                for i1 in range(nvar):  # loop over variables
                    aux = means[k, 1, i1]
                    means[k, 1, i1] = means[k, 0, i1]
                    means[k, 0, i1] = aux
                    for i2 in range(nvar):
                        aux = covars[k, 1, i1, i2]
                        covars[k, 1, i1, i2] = covars[k, 0, i1, i2]
                        covars[k, 0, i1, i2] = aux
        #         # for i in range(nx):
        #         #     for j in range(ny):
        #         a = [int(not labels[i]) for i in range(nij)]
        #
        # print(a[0:3])
        # print(labels[0:3])

        # trivar_plot_means(var, weights, means, covars, mean_tot, covars_tot, time, 'sortB')
        # trivar_plot_covars(var, weights, means, covars, mean_tot, covars_tot, time, 'sortB')
        # trivar_plot_weights(var, weights, means, covars, mean_tot, covars_tot, time, 'sortB')
        # print('')
        return means, covars, weights, labels


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
        # print('....', self.zrange.shape, type(self.zrange), type(self.zrange[:]))
        for k in range(nk):
            var[k] = self.zrange[k]
        # var[:] = self.zrange[:]

        # field_grp = rootgrp.createGroup('fields')
        # field_grp.createDimension('nx', nml['grid']['nx'])
        # field_grp.createDimension('ny', nml['grid']['ny'])
        # # field_grp.createDimension('nz', nml['grid']['nz'])
        # field_grp.createDimension('nz', nk)
        # labels = field_grp.createVariable('labels', 'f8', ('nx','ny','nz',))
        # labels.units = ' '
        return

    # cpdef add_field(self, name):
    #     rootgrp = nc.Dataset(self.path_plus_file, 'r+', format='NETCDF4')
    #     fieldgrp = rootgrp.groups['labels']
    #     fieldgrp.createVariable(name, 'f8', ('nx, ny, nz'))
    #     rootgrp.close()
    #     return

    # cpdef write_updrafts_field(self, path, file_name, double[:,:,:] data):
    cpdef write_updrafts_field(self, path, file_name, data, nml):
        print('write updrafts: ' + os.path.join(path, file_name))
        root = nc.Dataset(os.path.join(path, file_name), 'r+', format='NETCDF4')
        field_grp = root.createGroup('fields')
        nk = len(self.zrange)
        field_grp.createDimension('nx', nml['grid']['nx'])
        field_grp.createDimension('ny', nml['grid']['ny'])
        # field_grp.createDimension('nz', nml['grid']['nz'])
        field_grp.createDimension('nz', nk)
        labels = field_grp.createVariable('labels', 'f8', ('nx','ny','nz',))
        labels.units = ' '
        labels[:,:,:] = data
        # # field_grp = root.groups['fields']
        # # var = field_grp.variables['labels']
        # # var = root.groups['fields'].variables['labels']
        # # var = fieldgrp.variables[name]
        # # var[:,:,:] = np.array(data)
        # root.close()
        return


    cpdef create_statistics_file(self, path, file_name, time, ncomp, nvar, nz_):
        print('create statistics file: '+ path+', '+ file_name)
        # ncomp: number of Gaussian components in EM
        # nvar: number of variables of multi-variate Gaussian components
        rootgrp = nc.Dataset(os.path.join(path,file_name), 'w', format='NETCDF4')
        dimgrp = rootgrp.createGroup('dims')
        means_grp = rootgrp.createGroup('means')
        means_grp.createDimension('nz', nz_)
        means_grp.createDimension('ncomp', ncomp)
        means_grp.createDimension('nvar', nvar)
        cov_grp = rootgrp.createGroup('covariances')
        cov_grp.createDimension('nz', nz_)
        cov_grp.createDimension('ncomp', ncomp)
        cov_grp.createDimension('nvar', nvar)
        weights_grp = rootgrp.createGroup('weights')
        weights_grp.createDimension('nz', nz_)
        weights_grp.createDimension('EM2', 2)
        ts_grp = rootgrp.createGroup('time')
        ts_grp.createDimension('nt',len(time)-1)
        var = ts_grp.createVariable('t','f8',('nt'))
        for i in range(len(time)-1):
            var[i] = time[i+1]
        z_grp = rootgrp.createGroup('z-profile')
        z_grp.createDimension('nz', nz_)
        var = z_grp.createVariable('height', 'f8', ('nz'))
        for i in range(nz_):
            var[i] = self.zrange[i]
        rootgrp.close()
        return


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