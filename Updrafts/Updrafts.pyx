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

# !!!! sort PDF components, s.t. 0/1 always correspond to the same component (environment vs. updrafts)

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
        cdef:
            double [:,:,:] labels = np.zeros(shape=(nx, ny, nk))


        var_list = ['s', 'qt', 'temperature', 'ql']
        data = np.ndarray(shape=((nx * ny), nvar))
        for ncomp in ncomp_range:
            means_ = np.ndarray(shape=(nk, ncomp, nvar))
            covariance_ = np.zeros(shape=(nk, ncomp, nvar, nvar))
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

                '''(4) Normalise Data'''
                scaler = StandardScaler()
                data_all_norm = scaler.fit_transform(data_all)

                '''(5) Compute bivariate PDF'''
                #   (a) for (s,qt)
                #           ...
                #   (b) for (th_l,qt)
                # clf_thl_norm = Gaussian_bivariate(ncomp, data_all_norm, 'T', 'qt', np.int(d[0:-3]), iz * dz, path)
                # cpdef Gaussian_bivariate(ncomp_, data, var_name1, var_name2, time, z, path):
                clf_thl_norm = mixture.GaussianMixture(n_components=ncomp,covariance_type='full')
                clf_thl_norm.fit(data_all_norm)
                means_[k, :, :] = clf_thl_norm.means_[:, :]
                covariance_[k,:,:,:] = clf_thl_norm.covariances_[:,:,:]


                '''(6) Sort and look for labels'''
                labels_ = clf_thl_norm.predict(data_all_norm)
                # print('LABELS: ', type(labels))
                print('')
                print('Labels: ', np.count_nonzero(labels_), labels_.shape, data_all_norm.shape)
                # rearrange into 2D array (for only one data file)
                for i in range(nx):
                    for j in range(ny):
                        ij = i*ishift + j
                        labels[i,j,k] = labels_[ij]
            # !!!! sort PDF components, s.t. 0/1 always correspond to the same component (environment vs. updrafts)
            # self.write_updrafts_field(self.path_out, nc_file_name_labels, labels)
            self.write_updrafts_field(self.path_out, nc_file_name_labels, labels, nml)



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
        var[:] = self.zrange

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
    cpdef write_updrafts_field(self, path, file_name, double[:,:,:] data, nml):
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
        # field_grp = root.groups['fields']
        # var = field_grp.variables['labels']
        # var = root.groups['fields'].variables['labels']
        # var = fieldgrp.variables[name]
        # var[:,:,:] = np.array(data)
        root.close()
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