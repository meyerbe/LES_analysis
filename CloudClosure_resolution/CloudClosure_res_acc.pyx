import os, sys
#include 'parameters.pxi'
import json as simplejson
import time
import pylab as plt
import netCDF4 as nc
import numpy as np

from sklearn import mixture
from sklearn.preprocessing import StandardScaler

# TODO
# - compute ql_mean_field only once for all components (do in intialisation)
# - make nml --> self.nml
# - output  / dump statistics file
# - include labeling


import CC_thermodynamics_c
from CC_thermodynamics_c cimport LatentHeat, ClausiusClapeyron
from CC_thermodynamics_c import sat_adj_fromentropy, sat_adj_fromthetali

from plotting_functions import plot_PDF, plot_PDF_components
from plotting_functions import plot_error_vs_ncomp_ql, plot_error_vs_ncomp_cf, plot_abs_error


cdef class CloudClosure:
    def __init__(self):
        self.path_ref = ' '
        self.path_out = ' '
        self.p_ref = None
        self.z_ref = None
        self.zrange = None
        return



    cpdef initialize(self, krange, path, case_name):
        print('')
        print('--- Cloud Closure Scheme ---')
        print('nml: ', os.path.join(path, case_name+'.in'))
        print('')
        self.nml = simplejson.loads(open(os.path.join(path, case_name+'.in')).read())
        cdef:
            int nz = self.nml['grid']['nz']
            int dz = self.nml['grid']['dz']
        self.path_out = os.path.join(path, 'CloudClosure_res')


        '''Initialize Reference Pressure'''
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
        self.zrange = np.double(krange) * dz
        print('')
        print('zrange', np.double(krange) * dz)
        print('')


        '''Initialize Latent Heat and ClausiusClapeyron'''
        self.LH = CC_thermodynamics_c.LatentHeat(self.nml)
        self.CC = CC_thermodynamics_c.ClausiusClapeyron()
        # cdef:
        #     LatentHeat LH = CC_thermodynamics_c.LatentHeat(nml)
        #     ClausiusClapeyron CC = CC_thermodynamics_c.ClausiusClapeyron()
        self.CC.initialize(self.nml, self.LH)
        print('')


        return




    cpdef predict_pdf(self, files, path, int n_sample_, ncomp_range, Lx_, Ly_, dk_, int [:] krange_, nml):
        print('')
        print('--- PDF Prediction ---')
        print('')
        cdef extern from "thermodynamic_functions.h":
            inline double thetali_c(const double p0, const double T, const double qt, const double ql, const double qi, const double L)
        # time1 = time.clock()

        # ________________________________________________________________________________________
        '''(A) Set up files etc.'''
        cdef:
            int ncomp
            int nvar = 2
            int [:] krange = krange_
            int nk = len(krange)
            # double [:] dk_range = dk_range_
            int dk = dk_ + 1
            Py_ssize_t k, iz
            str d
            int dx = nml['grid']['dx']
            int dy = nml['grid']['dy']
            int dz = nml['grid']['dz']
            int nx = nml['grid']['nx']
            int ny = nml['grid']['ny']
            int nz = nml['grid']['nz']
            double [:] p_ref = np.zeros([nz],dtype=np.double,order='c')       # <type 'CloudClosure._memoryviewslice'>

        ''' Compute horizontal length scale'''
        if Lx_ > nx*dx or Ly_ > ny*dy:
            print('!!! WARNING !!!: Lx or Ly larger than domain size!')
            print('Program exit')
            sys.exit()
        cdef:
            int nx_ = np.floor(Lx_ / dx)
            int ny_ = np.floor(Ly_ / dy)
            int N_ = nx_*ny_
            int n_col_x = np.floor(nx*dx/Lx_)
            int n_col_y = np.floor(ny*dy/Ly_)
            int n_col_xy = n_col_x*n_col_y
        print('')
        print('Sampling domain size: ', Lx_, Ly_, '; nx_, ny_: ', nx_, ny_, 'N_:', N_)
        print('Total domain size: ', nx*dx, ny*dy, '; nx, ny, dx, dy: ', nx, ny, dx, dy)
        print('Number of columns (x, y, total): ', n_col_x, n_col_y, n_col_xy)
        print('Number of levels: ', dk_, dk)
        print('')

        # time_profile = read_in_netcdf('t', 'timeseries', self.path_ref)
        # nt_stats = time_profile.shape[0]
        # print('time in stats file: ', time_profile.shape, nt_stats, dt_stats)
        # ql_profile = read_in_netcdf('ql_mean', 'profiles', path_ref)

        # N = len(files)

        for k in range(nz):
            p_ref[k] = self.p_ref[k]


        cdef:
            LatentHeat LH = self.LH
            ClausiusClapeyron CC = self.CC

        # ________________________________________________________________________________________
        '''(B) Compute PDF f(s,qt) from LES data'''
        #       - read in fields
        #       - compute theta_l
        #       - compute PDFs f(s,qt), g(th_l, qt)

        # for PDF construction
        cdef:
            int i, j, ij
            int i_shift, j_shift, k_shift
            double [:,:,:] s_
            double [:,:,:] T_
            double [:,:,:] qt_
            double [:,:] qt = np.zeros([dk*n_col_xy*nx_*ny_,nk],dtype=np.double,order='c')         # <type 'CloudClosure._memoryviewslice'>
            double [:,:,:] ql_
            double [:,:] ql = np.zeros([dk*n_col_xy*nx_*ny_,nk],dtype=np.double,order='c')         # <type 'CloudClosure._memoryviewslice'>
            double [:] ql_all
            double qi_ = 0.0
            double [:,:] theta_l = np.zeros([dk*n_col_xy*nx_*ny_,nk],dtype=np.double,order='c')    # <type 'CloudClosure._memoryviewslice'>
            # double [:,:] data_all
            double Lv
            int count_ncol = 0
            int col_shift = 0

        # for PDF sampling
        cdef:
            int n_sample = n_sample_
            double [:] T_comp_thl = np.zeros([n_sample],dtype=np.double,order='c')
            double [:] ql_comp_thl = np.zeros([n_sample],dtype=np.double,order='c')
            # double [:,:] Th_l = np.zeros([n_sample, nvar], dtype=np.double, order='c')
            double [:] alpha_comp_thl = np.zeros(n_sample)

        # for Error Computation / <ql> intercomparison
        cdef:
            double [:] ql_mean_field = np.zeros(nk, dtype=np.double)        # computation from 3D LES field
            double [:] ql_mean_comp = np.zeros(nk, dtype=np.double)              # computation from PDF sampling
            double [:] cf_field = np.zeros(shape=(nk))                      # computation from 3D LES field
            double [:] cf_comp = np.zeros(nk, dtype=np.double)              # computation from PDF sampling
            int count_ncomp = 0
            double [:,:] error_ql = np.zeros(shape=(nk,len(ncomp_range)))
            double [:,:] rel_error_ql = np.zeros(shape=(nk,len(ncomp_range)))
            double [:,:] error_cf = np.zeros(shape=(nk,len(ncomp_range)))
            double [:,:] rel_error_cf = np.zeros(shape=(nk,len(ncomp_range)))



        '''(1) Read in Fields'''
        var_list = ['s', 'qt', 'temperature', 'ql']
        data = np.ndarray(shape=((n_col_xy * nx_ * ny_ * dk), nvar))
        for ncomp in ncomp_range:
            means_ = np.zeros(shape=(nk, ncomp, nvar))
            covariances_ = np.zeros(shape=(nk, ncomp, nvar, nvar))
            weights_ = np.zeros(shape=(nk, ncomp))
            '''(1) Statistics File'''
            tt = files[0][0:-3]
            nc_file_name_out = 'CC_alltime_ncomp'+str(ncomp)+'_Lx'+str(Lx_)+'Ly'+str(Ly_)+'_dz'+str((dk-1)*dz)\
                               +'_time'+str(tt)+'.nc'
            self.create_statistics_file(self.path_out, nc_file_name_out, str(tt), ncomp, nvar, nk)

            for k in range(nk):
                iz = krange[k]
                print('')
                print('- z = '+str(iz*dz)+ ', ncomp = '+str(ncomp)+' -')
                ql_all = np.ndarray(shape=(0))          # data averaged over all levels and time step
                data_all = np.ndarray(shape=(0, nvar))  # data averaged over all levels and time step

                for d in files:
                    print('time: d='+str(d))
                    nc_file_name = d
                    path_fields = os.path.join(path, 'fields', nc_file_name)
                    s_, qt_, T_, ql_ = read_in_fields('fields', var_list, path_fields)

                    '''(2) Compute liquid potential temperature from temperature and moisture'''
                    # for k_ in range(-dk+1,dk):
                    for k_ in range(0,dk):
                        k_shift = k_*n_col_xy*N_            # accumulation over multiple layers
                        print('layer: k_='+str(k_), str(iz), str(iz+k_), 'dk:', dk, dk_, 'k_shift: ', k_shift)
                        count_ncol = 0
                        for i_col in range(n_col_x):
                            for j_col in range(n_col_y):
                                col_shift = count_ncol*N_
                                for i in range(nx_):
                                    for j in range(ny_):
                                        i_shift = i_col*nx_ + i
                                        j_shift = j_col*ny_ + j
                                        ij = i*ny_ + j              # 2D --> 1D
                                        Lv = LH.L(T_[i_shift,j_shift,iz+k_],LH.Lambda_fp(T_[i_shift,j_shift,iz+k_]))
                                        theta_l[k_shift+col_shift+ij,k] = thetali_c(p_ref[iz+k_],
                                                    T_[i_shift,j_shift,iz+k_], qt_[i_shift,j_shift,iz+k_],
                                                                                    ql_[i_shift,j_shift,iz+k_], qi_, Lv)
                                        qt[k_shift+col_shift+ij,k] = qt_[i_shift,j_shift,iz+k_]
                                        ql[k_shift+col_shift+ij,k] = ql_[i_shift,j_shift,iz+k_]
                                        ql_mean_field[k] += ql_[i_shift,j_shift,iz+k_]
                                        if ql_[i_shift,j_shift,iz+k_] > 0.0:
                                            cf_field[k] += 1.0
                                count_ncol += 1
                                            # # test test test test
                                            # theta_l[k_*nx_*ny_+ij,k] = thetali_c(p_ref[iz], T_[i,j,iz], qt_[i,j,iz], ql_[i,j,iz], qi_, Lv)
                                            # qt[k_*nx_*ny_+ij,k] = qt_[i,j,iz]
                                            # ql[k_*nx_*ny_+ij,k] = ql_[i,j,iz]
                                            # ql_mean_field[k] += ql_[i,j,iz]
                                            # if ql_[i,j,iz] > 0.0:
                                            #     cf_field[k] += 1.0
                                            # save_name = 'ncomp'+str(ncomp)+'_dk'+str(dk) + '_z'+str((iz+k_)*dz)+'m'
                                            # scatter_data(theta_l[:,k], ql[:,k], 'thl', 'qt', dk, ncomp, err_ql, iz*dz, self.path_out, save_name)
                                            # print('')
                    del s_, ql_, qt_, T_
                    data[:, 0] = theta_l[:, k]
                    data[:, 1] = qt[:, k]
                    data_all = np.append(data_all, data, axis=0)
                    ql_all = np.append(ql_all, ql[:,k], axis=0)

                ql_mean_field[k] /= ( len(files)*(nx_*ny_)*dk )
                cf_field[k] /= ( len(files)*(nx_*ny_)*dk )

                '''(3) Normalise Data'''
                scaler = StandardScaler()
                data_all_norm = scaler.fit_transform(data_all)

                '''(4) Compute bivariate PDF'''
                #   (a) for (s,qt)
                #   (b) for (th_l,qt)
                clf_thl_norm = mixture.GaussianMixture(n_components=ncomp,covariance_type='full')
                clf_thl_norm.fit(data_all_norm)
                means_[k, :, :] = clf_thl_norm.means_[:, :]
                covariances_[k,:,:,:] = clf_thl_norm.covariances_[:,:,:]
                weights_[k,:] = clf_thl_norm.weights_[:]

                # '''(5) Compute Kernel-Estimate PDF '''
                # kde, kde_aux = Kernel_density_estimate(data, 'T', 'qt', np.int(d[0:-3]), iz * dz, path_in)
                # relative_entropy(data, clf, kde)
                # print('')

                # -------------------------------------------


                '''(D) Compute mean liquid water <ql> from PDF f(s,qt)'''
                #       1. sample (th_l, qt) from PDF (Monte Carlo ???
                #       2. compute ql for samples
                #       3. consider ensemble average <ql> = domain mean representation???
                '''(1) Draw samples'''
                Th_l_norm, y_norm = clf_thl_norm.sample(n_samples=n_sample)
                '''(2) Rescale theta_l and qt'''
                Th_l = scaler.inverse_transform(Th_l_norm)      # Inverse Normalisation

                '''(3) Compute ql (saturation adjustment) & Cloud Fraction '''
                # print('ql_mean_comp[k], k', k, ql_mean_comp[k])
                # print('Th_l_norm[0]', np.amin(Th_l_norm[:,0]), np.amax(Th_l_norm[:,0]), 'Th_l_norm[1]', np.amin(Th_l_norm[:,1]), np.amax(Th_l_norm[:,1]))
                # print('Th_l[0]', np.amin(Th_l[:,0]), np.amax(Th_l[:,0]), 'Th_l[1]', np.amin(Th_l[:,1]), np.amax(Th_l[:,1]))
                # print('theta_l', np.amin(theta_l[:,k]), np.amax(theta_l[:,k]), 'qt', np.amin(qt[:,k]), np.amax(qt[:,k]))
                # print('data_all', np.amin(data_all[:,0]), np.amax(data_all[:,0]), np.amin(data_all[:,0]), np.amax(data_all[:,1]))
                # print('data_all_norm', np.amin(data_all_norm[:,0]), np.amax(data_all_norm[:,0]), np.amin(data_all_norm[:,0]), np.amax(data_all_norm[:,1]))
                # print('')
                for i in range(n_sample-2):
                    # ??? ok to use same reference pressure for all ik+k_ points?
                    # T_comp_thl, ql_comp_thl, alpha_comp_thl = sat_adj_fromthetali(p_ref[iz], Th_l[i, 0], Th_l[i, 1], CC, LH)
                    T_comp_thl[i], ql_comp_thl[i], alpha_comp_thl[i] = sat_adj_fromthetali(p_ref[iz], Th_l[i, 0], Th_l[i, 1], CC, LH)
                    ql_mean_comp[k] = ql_mean_comp[k] + ql_comp_thl[i]
                    if ql_comp_thl[i] > 0:
                        cf_comp[k] += 1
            #     print('ql_mean_comp[k], k', k, ql_mean_comp[k], T_comp_thl[i], ql_comp_thl[i])
                ql_mean_comp[k] = ql_mean_comp[k] / n_sample
                cf_comp[k] = cf_comp[k] / n_sample
                error_ql[k,count_ncomp] = ql_mean_comp[k] - ql_mean_field[k]
                error_cf[k,count_ncomp] = cf_comp[k] - cf_field[k]
                if ql_mean_field[k] > 0.0:
                    rel_error_ql[k,count_ncomp] = (ql_mean_comp[k] - ql_mean_field[k]) / ql_mean_field[k]
                if cf_field[k] > 0.0:
                    rel_error_cf[k,count_ncomp] = (cf_comp[k] - cf_field[k]) / cf_field[k]

                print('')
                print('<ql> from CloudClosure Scheme: ', ql_mean_comp[k])
                print('<ql> from ql fields: ', ql_mean_field[k])
                print('error (<ql>_CC - <ql>_field): '+str(error_ql[k,count_ncomp]))
                print('rel err: '+ str(rel_error_ql[k,count_ncomp]))
                print('')
                print('CF from Cloud Closure Scheme: ', cf_comp[k])
                print('CF from ql fields: ', cf_field[k])
                print('error: '+str(error_cf[k,count_ncomp]))
                print('rel error: ', rel_error_cf[k,count_ncomp])
                print('')

                '''(E) Plotting'''
                save_name = 'PDF_figures_'+str(iz*dz)+'m'+'_ncomp'+str(ncomp)+'_Lx'+str(Lx_)+'_dk'+str(dk-1)
                plot_PDF(data_all, data_all_norm, 'thl', 'qt', clf_thl_norm, dk, ncomp, error_ql[k,count_ncomp], iz*dz, self.path_out, save_name)
                # plot_samples('norm', data_all_norm, ql_all[:], Th_l_norm, ql_comp_thl, 'thl', 'qt', scaler, ncomp, iz*dz, path_out)
                # plot_samples('original', data_all, ql_all[:], Th_l, ql_comp_thl, 'thl', 'qt', scaler, ncomp, iz*dz, path_out)
                # plot_hist(ql, path_out)
                print('')
            count_ncomp += 1
            plot_error_vs_ncomp_ql(error_ql, rel_error_ql, n_sample, ql_mean_field, ql_mean_comp, ncomp_range, krange, dz, Lx_, dk-1, self.path_out)
            plot_error_vs_ncomp_cf(error_cf, rel_error_cf, n_sample, cf_field, ncomp_range, krange, dz, Lx_, dk-1, self.path_out)
            plot_abs_error(error_ql, ql_mean_field, error_cf, cf_field, n_sample, ncomp_range, krange, dz, Lx_, dk-1, self.path_out)

            plot_PDF_components(means_, covariances_, weights_, ncomp, krange, dz, Lx_, dk-1, self.path_out)
        return




#     cpdef sample_pdf(self, data, clf, double ql_mean_ref, double cf_ref, double pref,
#                             ClausiusClapeyron CC, LatentHeat LH, ncomp_range, n_sample, nml):
#         # 1. sample (th_l, qt) from PDF (Monte Carlo ???
#         # 2. compute ql for samples
#         # 3. consider ensemble average <ql> = domain mean representation???
#
#         cdef:
#             int i, k
#             int nvar = 2
#             int nz = nml['grid']['nz']
#             int nk = ql_mean_ref.shape[0]
#
#         # for PDF sampling
#         cdef:
#             double [:] T_comp_thl = np.zeros([n_sample],dtype=np.double,order='c')
#             double [:] ql_comp_thl = np.zeros([n_sample],dtype=np.double,order='c')
#             double [:,:] Th_l = np.zeros([n_sample, nvar], dtype=np.double, order='c')
#             double [:] alpha_comp_thl = np.zeros(n_sample)
#
#         # for Error Computation / <ql> intercomparison
#         cdef:
#             int count_ncomp = 0
#             # double [:] ql_mean_field = np.zeros(nk, dtype=np.double)        # computation from 3D LES field
#             # double [:] cf_field = np.zeros(shape=(nk))                      # computation from 3D LES field
#             double ql_mean_comp = 0.0
#             double cf_comp = 0.0
#             double error_ql = 0.0
#             double rel_error_ql = 0.0
#             double error_cf = 0.0
#             double rel_error_cf = 0.0
#
#
#
#         '''(1) Draw samples'''
#         scaler = StandardScaler()
#         data_nrom = scaler.fit_transform(data)
#         Th_l_norm, y_norm = clf.sample(n_samples=n_sample)
#         '''(2) Rescale theta_l and qt'''
#         Th_l = scaler.inverse_transform(Th_l_norm)      # Inverse Normalisation
#
#         '''(3) Compute ql (saturation adjustment) & Cloud Fraction '''
#         for i in range(n_sample-2):
#             T_comp_thl[i], ql_comp_thl[i], alpha_comp_thl[i] = sat_adj_fromthetali(pref, Th_l[i, 0], Th_l[i, 1], CC, LH)
#             ql_mean_comp = ql_mean_comp + ql_comp_thl[i]
#             if ql_comp_thl[i] > 0:
#                 cf_comp += 1
#         ql_mean_comp = ql_mean_comp / n_sample
#         cf_comp = cf_comp / n_sample
#         error_ql = ql_mean_comp - ql_mean_ref
#         error_cf = cf_comp- cf_ref
#         if ql_mean_ref > 0.0:
#             rel_error_ql = (ql_mean_comp - ql_mean_ref) / ql_mean_ref
#         if cf_ref > 0.0:
#             rel_error_cf = (cf_comp - cf_ref) / cf_ref
#
#         print('')
#         print('<ql> from CloudClosure Scheme: ', ql_mean_comp)
#         print('<ql> from ql fields: ', ql_mean_ref)
#         print('error (<ql>_CC - <ql>_field): '+ str(error_ql))
#         print('rel err: '+ str(rel_error_ql))
#         print('')
#         print('CF from Cloud Closure Scheme: ', cf_comp)
#         print('CF from ql fields: ', cf_ref)
#         print('error: '+str(error_cf))
#         print('rel error: ', rel_error_cf)
#         print('')
#
#         # '''(E) Plotting'''
#         # save_name = 'PDF_figures_'+str(iz*dz)+'m'+'_ncomp'+str(ncomp)+'_dz'+str(dk-1)
#         # plot_PDF(data_all, data_all_norm, 'thl', 'qt', clf_thl_norm, dk, ncomp, error_ql[k,count_ncomp], iz*dz, self.path_out, save_name)
#         # # plot_samples('norm', data_all_norm, ql_all[:], Th_l_norm, ql_comp_thl, 'thl', 'qt', scaler, ncomp, iz*dz, path_out)
#         # # plot_samples('original', data_all, ql_all[:], Th_l, ql_comp_thl, 'thl', 'qt', scaler, ncomp, iz*dz, path_out)
#         # # plot_hist(ql, path_out)
#         # print('')
#
#         return ql_mean_comp, error_ql, rel_error_ql, cf_comp, error_cf, rel_error_cf








    #----------------------------------------------------------------------
    #----------------------------------------------------------------------
    def create_statistics_file(self, path, file_name, time, ncomp, nvar, nz_):
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



#----------------------------------------------------------------------
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