import os, sys
#include 'parameters.pxi'
import json as simplejson
import time
import pylab as plt
import netCDF4 as nc
import numpy as np
cimport numpy as np

from sklearn import mixture
from sklearn.preprocessing import StandardScaler

''' OUTPUT '''
'(1) Output per ncomp: PDF_cond_<type>_ncomp<x>_Lx<Lx>Ly<Ly>_dz<dz>_time<tt>.nc'
#     - PDF parameters:
#         - means
#         - covariances
#         - weights
#
#     - Profiles:
#         nk = len(krange)
#         Ncomp = len(ncomp_range)
#         - n_environment = (nk):     # of 'environmental' points according to tracers by Colleen
#         - n_updraft = (nk):         # of 'updraft' points according to tracers by Colleen
#
#         - ql_mean_comp = (nk):            <ql>_{env, comp}environment computed from PDF sampling
#         - ql_mean_environment = (nk):     <ql>_environment from conditional sampling of LES data (i.e. updraft points are masked)
#         - ql_mean_updraft = (nk):         <ql>_updraft from conditional sampling of LES data (i.e. only updraft points)
#         - ql_mean_domain = (nk):          <ql> from total LES domain (no conditional sampling)
#         - cf_comp = (nk):                 CF_environment computed from PDF sampling
#         - cf_environment = (nk):          CF_environment from conditional sampling of LES data (i.e. updraft points are masked)
#         - cf_updraft = (nk):              CF_updraft from conditional sampling of LES (i.e. only updraft points)
#         - cf_domain = (nk):               CF from total LES domain (no conditional sampling)
#
#         - error_ql_env = (nk):          e = ql_mean_comp - ql_mean_env
#         - error_ql_domain = (nk):       e = ql_mean_comp - ql_mean_domain
#         - rel_error_ql_env = (nk):      e_ = (ql_mean_comp - ql_mean_env) / ql_mean_env
#         - rel_error_ql_domain = (nk):   e_ = (ql_mean_comp - ql_mean_domain) / ql_mean_domain
#
#         - error_cf_env = (nk):          e = cf_comp - cf_env
#         - error_cf_domain = (nk):       e = cf_comp - cf_domain
#         - rel_error_cf_env = (nk):      e_ = (cf_comp - cf_env) / cf_env
#         - rel_error_cf_domain = (nk):   e_ = (cf_comp - cf_domain) / cf_domain


'(2) Output for all: PDF_cond_<type>_allncomp_Lx<Lx>Ly<Ly>_dz<dz>_time<tt>.nc'
#         nk = len(krange)
#         Ncomp = len(ncomp_range)

#     - Profiles:
#         - krange: k-values [-] where PDF-model computed
#         - zrange: z-values [m] where PDF-model computed
#
#         - n_environment = (nk):     # of 'environmental' points according to tracers by Colleen
#         - n_updraft = (nk):         # of 'updraft' points according to tracers by Colleen
#
#         - ql_mean_environment = (nk):     <ql>_environment from conditional sampling of LES data (i.e. updraft points are masked)
#         - ql_mean_updraft = (nk):         <ql>_updraft from conditional sampling of LES data (i.e. only updraft points)
#         - ql_mean_domain = (nk):          <ql> from total LES domain (no conditional sampling)
#         - cf_comp = (nk):                 CF_environment computed from PDF sampling
#         - cf_environment = (nk):          CF_environment from conditional sampling of LES data (i.e. updraft points are masked)
#         - cf_updraft = (nk):              CF_updraft from conditional sampling of LES (i.e. only updraft points)
#         - cf_domain = (nk):               CF from total LES domain (no conditional sampling)
#
#         - error_ql_env = (nk, Ncomp):          e = ql_mean_comp - ql_mean_env
#         - error_ql_domain = (nk, Ncomp):       e = ql_mean_comp - ql_mean_domain
#         - rel_error_ql_env = (nk, Ncomp):      e_ = (ql_mean_comp - ql_mean_env) / ql_mean_env
#         - rel_error_ql_domain = (nk, Ncomp):   e_ = (ql_mean_comp - ql_mean_domain) / ql_mean_domain
#
#         - error_cf_env = (nk, Ncomp):          e = cf_comp - cf_env
#         - error_cf_domain = (nk, Ncomp):       e = cf_comp - cf_domain
#         - rel_error_cf_env = (nk, Ncomp):      e_ = (cf_comp - cf_env) / cf_env
#         - rel_error_cf_domain = (nk, Ncomp):   e_ = (cf_comp - cf_domain) / cf_domain




import CC_thermodynamics_c
from CC_thermodynamics_c cimport LatentHeat, ClausiusClapeyron
from CC_thermodynamics_c import sat_adj_fromentropy, sat_adj_fromthetali

from plotting_functions import plot_PDF#, plot_PDF_components
# from plotting_functions import plot_error_vs_ncomp_ql, plot_error_vs_ncomp_cf, plot_abs_error


cdef class PDF_conditional:
    def __init__(self):
        self.path_ref = ' '
        self.path_out = ' '
        self.p_ref = None
        self.z_ref = None
        self.zrange = None
        self.krange = None
        return




    cpdef initialize(self, krange, zrange, path, case_name):
        print('')
        print('--- PDF conditional (not anomaly) ---')
        print('nml: ', os.path.join(path, case_name+'.in'))
        print('')
        self.nml = simplejson.loads(open(os.path.join(path, case_name+'.in')).read())
        cdef:
            int nz = self.nml['grid']['nz']
            int dz = self.nml['grid']['dz']
        self.path_out = os.path.join(path, 'PDF_cond')

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

        self.krange = krange
        try:
            self.z_ref = read_in_netcdf('z_half', 'profiles', self.path_ref)
        except:
            self.z_ref = read_in_netcdf('z', 'reference', self.path_ref)
        self.zrange = zrange
        print('')
        print('zrange: ' + str(np.double(krange) * dz))
        print('z_ref: ' + str(self.z_ref[:]))
        print('')


        '''Initialize Latent Heat and ClausiusClapeyron'''
        self.LH = CC_thermodynamics_c.LatentHeat(self.nml)
        self.CC = CC_thermodynamics_c.ClausiusClapeyron()
        self.CC.initialize(self.nml, self.LH)
        print('')

        return




    cpdef predict_pdf(self, files, path, int n_sample_, ncomp_range, Lx_, Ly_, dk_, int [:] krange_, nml):
        print('')
        print('--- PDF Prediction ---')
        print('dk='+str(dk_), 'Lx='+str(Lx_))
        print('______________________')
        print('')
        cdef extern from "thermodynamic_functions.h":
            inline double thetali_c(const double p0, const double T, const double qt, const double ql, const double qi, const double L)


        # ________________________________________________________________________________________
        '''(A) Set up files etc.'''
        cdef:
            int ncomp
            int nvar = 2
            int [:] krange = krange_
            int nk = len(krange)
            int dk = dk_ + 1
            Py_ssize_t k, iz, k_
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
            print('Domain Size: ', nx*dx, ny*dy, ' Lx: ', Lx_, Ly_)
            print('Program exit')
            sys.exit()
        cdef:
            int nx_ = np.floor(Lx_ / dx)
            int ny_ = np.floor(Ly_ / dy)
        print('')
        print('Sampling domain size: ', Lx_, Ly_, '; nx_, ny_: ', nx_, ny_)
        print('Total domain size: ', nx*dx, ny*dy, '; nx, ny, dx, dy: ', nx, ny, dx, dy)
        print('Number of levels: ', dk_, dk)
        print('')

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
            int i, j
            double [:,:,:] s_
            double [:,:,:] T_
            double [:,:,:] qt_
            double [:] qt = np.ndarray(shape=(0), dtype=np.double)
            double [:,:,:] ql_
            double [:] ql = np.ndarray(shape=(0), dtype=np.double)
            double [:] ql_all
            double qi_ = 0.0
            double [:] theta_l = np.ndarray(shape=(0), dtype=np.double)
            double [:,:] data_all
            double Lv

        # for PDF sampling
        cdef:
            int n_sample = n_sample_
            double [:] T_comp_thl = np.zeros([n_sample],dtype=np.double,order='c')
            double [:] ql_comp_thl = np.zeros([n_sample],dtype=np.double,order='c')
            double [:,:] Th_l = np.zeros([n_sample, nvar], dtype=np.double, order='c')

        # for Error Computation / <ql> intercomparison
        cdef:
            double [:, :] ql_mean_comp = np.zeros(shape=(nk,len(ncomp_range)), dtype=np.double)                 # computation from PDF sampling
            double [:, :] cf_comp = np.zeros(shape=(nk,len(ncomp_range)), dtype=np.double)                      # computation from PDF sampling
            int count_ncomp = 0
            double [:,:] error_ql_domain = np.zeros(shape=(nk,len(ncomp_range)))
            double [:,:] error_ql_env = np.zeros(shape=(nk,len(ncomp_range)))
            # double [:,:] error_ql_up = np.zeros(shape=(nk,len(ncomp_range)))
            double [:,:] rel_error_ql_domain = np.zeros(shape=(nk,len(ncomp_range)))
            double [:,:] rel_error_ql_env = np.zeros(shape=(nk,len(ncomp_range)))
            # double [:,:] rel_error_ql_up = np.zeros(shape=(nk,len(ncomp_range)))
            double [:,:] error_cf_domain = np.zeros(shape=(nk,len(ncomp_range)))
            double [:,:] error_cf_env = np.zeros(shape=(nk,len(ncomp_range)))
            # double [:,:] error_cf_up = np.zeros(shape=(nk,len(ncomp_range)))
            double [:,:] rel_error_cf_domain = np.zeros(shape=(nk,len(ncomp_range)))
            double [:,:] rel_error_cf_env = np.zeros(shape=(nk,len(ncomp_range)))
            # double [:,:] rel_error_cf_up = np.zeros(shape=(nk,len(ncomp_range)))

        # for Updraft / Environment Decompositoin
        cdef:
            # int n_env = 0
            # int n_updraft = 0
            int [:] n_env = np.zeros(nk, dtype=np.int32)
            int [:] n_updraft = np.zeros(nk, dtype=np.int32)
            double [:] ql_mean_domain = np.zeros(nk, dtype=np.double)        # computation from 3D LES field
            double [:] ql_mean_env = np.zeros(nk, dtype=np.double)
            double [:] ql_mean_updraft = np.zeros(nk, dtype=np.double)
            double [:] cf_domain = np.zeros(shape=(nk))                      # computation from 3D LES field
            double [:] cf_env = np.zeros(nk, dtype=np.double)
            double [:] cf_updraft = np.zeros(nk, dtype=np.double)



        '''(1) Read in Fields'''
        var_list = ['s', 'qt', 'temperature', 'ql']
        type_list = ['Couvreux']
        # type_list = ['Couvreux', 'Coherent']
        # ________________________________________________________________________________________
        # Read in Tracer Labels
        path_tracers = os.path.join(path, 'tracer_fields')
        tt = files[0][0:-3]

        ''' (2) Compute PDF Model '''
        for type_ in type_list:
            labels_tracers = read_in_updrafts_colleen(type_, tt, path_tracers)
            print('labels_tr', labels_tracers.shape, np.count_nonzero(labels_tracers))
            # ncomp = ncomp_range[0]
            for ncomp in ncomp_range:
                means_ = np.zeros(shape=(nk, ncomp, nvar))
                covariances_ = np.zeros(shape=(nk, ncomp, nvar, nvar))
                weights_ = np.zeros(shape=(nk, ncomp))
                n_env = np.zeros(nk, dtype=np.int32)
                n_updraft = np.zeros(nk, dtype=np.int32)
                ql_mean_domain = np.zeros(nk, dtype=np.double)
                ql_mean_env = np.zeros(nk, dtype=np.double)
                ql_mean_updraft = np.zeros(nk, dtype=np.double)
                cf_domain = np.zeros(nk, dtype=np.double)
                cf_env = np.zeros(nk, dtype=np.double)
                cf_updraft = np.zeros(nk, dtype=np.double)

                ''' (b) initialize Statistics File '''
                nc_file_name_out = 'PDF_cond_' + type_ + '_ncomp'+str(ncomp)+'_Lx'+str(Lx_)+'Ly'+str(Ly_)+'_dz'+str((dk)*dz) \
                                   +'_time'+str(tt)+'.nc'
                self.create_statistics_file(self.path_out, nc_file_name_out, str(tt), ncomp, nvar, nk)




                '''(2) Compute liquid potential temperature from temperature and moisture'''
                for k in range(nk):
                    iz = krange[k]
                    # n_env = 0
                    # n_updraft = 0
                    print('')
                    print('-- z = '+str(iz*dz)+ ', ncomp = '+str(ncomp)+' --')
                    theta_l = np.ndarray(shape=(0), dtype=np.double)
                    qt = np.ndarray(shape=(0), dtype=np.double)
                    ql = np.ndarray(shape=(0), dtype=np.double)
                    ql_all = np.ndarray(shape=(0))          # data averaged over all levels and time step
                    data_all = np.ndarray(shape=(0, nvar))  # data averaged over all levels and time step


                    for d in files:
                        # pass
                        # d = files[0]
                        print('time: d='+str(d))
                        nc_file_name = d
                        path_fields = os.path.join(path, 'fields', nc_file_name)
                        print(path_fields)
                        s_, qt_, T_, ql_ = read_in_fields('fields', var_list, path_fields)
                        for k_ in range(0,dk):
                            print('----- k_='+str(k_), ' dk='+str(dk))
                            # k_ = 0
                            for i in range(nx_):
                                for j in range(ny_):
                                    ql_mean_domain[k] += ql_[i,j,iz+k_]
                                    if ql_[i,j,iz+k_] > 0.0:
                                        cf_domain[k] += 1.0
                                    if labels_tracers[i,j,iz+k_] == 0:
                                        n_env[k] += 1
                                        Lv = LH.L(T_[i,j,iz+k_],LH.Lambda_fp(T_[i,j,iz+k_]))
                                        aux = thetali_c(p_ref[iz+k_], T_[i,j,iz+k_], qt_[i,j,iz+k_], ql_[i,j,iz+k_], qi_, Lv)
                                        theta_l = np.append(theta_l, aux)
                                        qt = np.append(qt, qt_[i,j,iz+k_])
                                        ql = np.append(ql, ql_[i,j,iz+k_])
                                        ql_mean_env[k] += ql_[i,j,iz+k_]
                                        if ql_[i,j,iz+k_] > 0.0:
                                            cf_env[k] += 1.0
                                    else:
                                        n_updraft[k] += 1
                                        ql_mean_updraft[k] += ql_[i,j,iz+k_]
                                        if ql_[i,j,iz+k_] > 0.0:
                                            cf_updraft[k] += 1.0

                    # print('n_env='+str(n_env[k]), theta_l.shape)
                    # del s_, ql_, qt_, T_
                    data = np.ndarray(shape=(n_env[k], nvar))
                    # print('..')
                    # print('n_up, n_env, dk', n_updraft[k], n_env[k], dk)
                    # print(n_env)
                    # print(theta_l.shape, qt.shape, data.shape)
                    # print('..', i, nx_, j, ny_, nx_*ny_, np.count_nonzero(labels_tracers[0:nx_,0:ny_,iz]))
                    data[:, 0] = theta_l[:]
                    data[:, 1] = qt[:]
                    data_all = np.append(data_all, data, axis=0)
                    ql_all = np.append(ql_all, ql[:], axis=0)

                    # ql_mean_domain[k] /= ( len(files)*dk*(nx_*ny_) )
                    # cf_domain[k] /= ( len(files)*dk*(nx_*ny_) )
                    ql_mean_domain[k] /= ( 1*dk*(nx_*ny_) )
                    cf_domain[k] /= ( 1*dk*(nx_*ny_) )
                    ql_mean_env[k] /= n_env[k]
                    cf_env[k] /= n_env[k]
                    if n_updraft[k] > 0.0:
                        ql_mean_updraft[k] /= n_updraft[k]
                        cf_updraft[k] /= n_updraft[k]

                    '''(3) Normalise Data'''
                    scaler = StandardScaler()
                    data_all_norm = scaler.fit_transform(data_all)

                    '''(4) Compute bivariate PDF for (th_l,qt)'''
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
                    print('ql_mean_comp[k], k', k, ql_mean_comp[k, count_ncomp])
                    for i in range(n_sample-2):
                        # ??? ok to use same reference pressure for all ik+k_ points?
                        T_comp_thl[i], ql_comp_thl[i] = sat_adj_fromthetali(p_ref[iz], Th_l[i, 0], Th_l[i, 1], CC, LH)
                        ql_mean_comp[k, count_ncomp] += ql_comp_thl[i]
                        if ql_comp_thl[i] > 0:
                            cf_comp[k, count_ncomp] += 1
                    print('ql_mean_comp[k], k', k, ql_mean_comp[k, count_ncomp], T_comp_thl[i], ql_comp_thl[i])
                    ql_mean_comp[k, count_ncomp] = ql_mean_comp[k, count_ncomp] / (n_sample-2)
                    cf_comp[k, count_ncomp] = cf_comp[k, count_ncomp] / (n_sample-2)
                    error_ql_domain[k,count_ncomp] = ql_mean_comp[k, count_ncomp] - ql_mean_domain[k]
                    error_cf_domain[k,count_ncomp] = cf_comp[k, count_ncomp] - cf_domain[k]
                    error_ql_env[k,count_ncomp] = ql_mean_comp[k, count_ncomp] - ql_mean_env[k]
                    error_cf_env[k,count_ncomp] = cf_comp[k, count_ncomp] - cf_env[k]
                    if ql_mean_domain[k] > 0.0:
                        rel_error_ql_domain[k,count_ncomp] = error_ql_domain[k, count_ncomp] / ql_mean_domain[k]
                    if cf_domain[k] > 0.0:
                        rel_error_cf_domain[k,count_ncomp] = error_cf_domain[k, count_ncomp] / cf_domain[k]
                    if ql_mean_env[k] > 0.0:
                        rel_error_ql_env[k,count_ncomp] = error_ql_env[k, count_ncomp] / ql_mean_env[k]
                    if cf_env[k] > 0.0:
                        rel_error_cf_env[k,count_ncomp] = error_cf_env[k, count_ncomp] / cf_env[k]

                    print('')
                    print('<ql> from CloudClosure Scheme: ', ql_mean_comp[k, count_ncomp])
                    print('<ql> from ql fields (env):     ', ql_mean_env[k])
                    print('<ql> from ql fields (domain):  ', ql_mean_domain[k])
                    print('error env (<ql>_CC - <ql>_env):       '+str(error_ql_env[k,count_ncomp]))
                    print('error domain (<ql>_CC - <ql>_domain): '+str(error_ql_domain[k,count_ncomp]))
                    print('rel err env:    '+ str(rel_error_ql_env[k,count_ncomp]))
                    print('rel err domain: '+ str(rel_error_ql_domain[k,count_ncomp]))
                    print('')
                    print('CF from Cloud Closure Scheme: ', cf_comp[k, count_ncomp])
                    print('CF from ql fields (env):      ', cf_env[k])
                    print('CF from ql fields (domain):   ', cf_domain[k])
                    print('error env:    '+str(error_cf_env[k,count_ncomp]))
                    print('error domain: '+str(error_cf_domain[k,count_ncomp]))
                    print('rel error env:    ', rel_error_cf_env[k,count_ncomp])
                    print('rel error domain: ', rel_error_cf_domain[k,count_ncomp])
                    print('')

                    '''(E) Plotting'''
                    save_name = 'PDF_figures_'+str(iz*dz)+'m'+ type_ + '_ncomp'+str(ncomp)+'_Lx'+str(Lx_)+'_dk'+str(dk-1)
                    plot_PDF(data_all, data_all_norm, ql_all, 'thl', 'qt', clf_thl_norm, dk, ncomp, error_ql_env[k,count_ncomp], iz*dz, Lx_, self.path_out, save_name)
                    # # save_name = 'sample_figure_'+'ncomp'+str(ncomp)+'_Lx' + str(Lx_) +'_dk'+str(dk-1) + '_z' + str(iz*dz)+'m.png'
                    # # plot_samples('norm', data_all_norm, ql_all[:], Th_l_norm, ql_comp_thl, 'thl', 'qt', scaler, ncomp, iz*dz, path_out)
                    # # plot_samples('original', data_all, ql_all[:], Th_l, ql_comp_thl, 'thl', 'qt', scaler, ncomp, iz*dz, path_out)
                    # # plot_hist(ql, path_out)
                    # print('')

                    # plot_error_vs_ncomp_ql(error_ql, rel_error_ql, n_sample, ql_mean_field, ql_mean_comp, ncomp_range, krange, dz, Lx_, dk-1, self.path_out)
                    # plot_error_vs_ncomp_cf(error_cf, rel_error_cf, n_sample, cf_field, ncomp_range, krange, dz, Lx_, dk-1, self.path_out)
                    # plot_abs_error(error_ql, ql_mean_field, error_cf, cf_field, n_sample, ncomp_range, krange, dz, Lx_, dk-1, self.path_out)
                    #
                    # plot_PDF_components(means_, covariances_, weights_, ncomp, krange, dz, Lx_, dk-1, self.path_out)

                del s_, ql_, qt_, T_

                '''(F) Save Gaussian Mixture PDFs '''
                print('')
                print('Dumping files: '+ self.path_out)
                #     print('ncomp', ncomp, ncomp_range)
                dump_variable(os.path.join(self.path_out, nc_file_name_out), 'means', means_, 'qtT', ncomp, nvar, nk)
                dump_variable(os.path.join(self.path_out, nc_file_name_out), 'covariances', covariances_, 'qtT', ncomp, nvar, nk)
                dump_variable(os.path.join(self.path_out, nc_file_name_out), 'weights', weights_, 'qtT', ncomp, nvar, nk)
                dump_variable(os.path.join(self.path_out, nc_file_name_out), 'profiles', np.asarray(n_env[:]), 'n_environment', ncomp, nvar, nk)
                dump_variable(os.path.join(self.path_out, nc_file_name_out), 'profiles', np.asarray(n_updraft[:]), 'n_updraf', ncomp, nvar, nk)
                dump_variable(os.path.join(self.path_out, nc_file_name_out), 'profiles', np.asarray(ql_mean_comp[:, count_ncomp]), 'ql_mean_pdf', ncomp, nvar, nk)
                dump_variable(os.path.join(self.path_out, nc_file_name_out), 'profiles', np.asarray(ql_mean_env[:]), 'ql_mean_environment', ncomp, nvar, nk)
                dump_variable(os.path.join(self.path_out, nc_file_name_out), 'profiles', np.asarray(ql_mean_updraft[:]), 'ql_mean_updraft', ncomp, nvar, nk)
                dump_variable(os.path.join(self.path_out, nc_file_name_out), 'profiles', np.asarray(ql_mean_domain[:]), 'ql_mean_domain', ncomp, nvar, nk)
                dump_variable(os.path.join(self.path_out, nc_file_name_out), 'profiles', np.asarray(cf_comp[:, count_ncomp]), 'cf_pdf', ncomp, nvar, nk)
                dump_variable(os.path.join(self.path_out, nc_file_name_out), 'profiles', np.asarray(cf_env[:]), 'cf_environment', ncomp, nvar, nk)
                dump_variable(os.path.join(self.path_out, nc_file_name_out), 'profiles', np.asarray(cf_updraft[:]), 'cf_updraft', ncomp, nvar, nk)
                dump_variable(os.path.join(self.path_out, nc_file_name_out), 'profiles', np.asarray(cf_domain[:]), 'cf_domain', ncomp, nvar, nk)
                dump_variable(os.path.join(self.path_out, nc_file_name_out), 'error', np.asarray(error_ql_env[:,count_ncomp]), 'error_ql_env', ncomp, nvar, nk)
                dump_variable(os.path.join(self.path_out, nc_file_name_out), 'error', np.asarray(error_ql_domain[:,count_ncomp]), 'error_ql_domain', ncomp, nvar, nk)
                dump_variable(os.path.join(self.path_out, nc_file_name_out), 'error', np.asarray(rel_error_ql_env[:,count_ncomp]), 'rel_error_ql_env', ncomp, nvar, nk)
                dump_variable(os.path.join(self.path_out, nc_file_name_out), 'error', np.asarray(rel_error_ql_domain[:,count_ncomp]), 'rel_error_ql_domain', ncomp, nvar, nk)
                dump_variable(os.path.join(self.path_out, nc_file_name_out), 'error', np.asarray(error_cf_env[:,count_ncomp]), 'error_cf_env', ncomp, nvar, nk)
                dump_variable(os.path.join(self.path_out, nc_file_name_out), 'error', np.asarray(error_cf_domain[:,count_ncomp]), 'error_cf_domain', ncomp, nvar, nk)
                dump_variable(os.path.join(self.path_out, nc_file_name_out), 'error', np.asarray(rel_error_cf_env[:,count_ncomp]), 'rel_error_cf_env', ncomp, nvar, nk)
                dump_variable(os.path.join(self.path_out, nc_file_name_out), 'error', np.asarray(rel_error_cf_domain[:,count_ncomp]), 'rel_error_cf_domain', ncomp, nvar, nk)

                count_ncomp += 1

            error_file_name = 'PDF_cond_' + type_ + '_allcomp_Lx'+str(np.floor(Lx_))+'Ly'+str(np.floor(Ly_))+'_dz'+str((dk)*dz)+'_time'+str(tt)+'.nc'
            self.dump_error_file(self.path_out, error_file_name, str(tt), ncomp_range, nvar, nk,
                                 np.asarray(ql_mean_env), np.asarray(ql_mean_updraft), np.asarray(ql_mean_domain),
                                 np.asarray(cf_env), np.asarray(cf_updraft), np.asarray(cf_domain),
                                 np.asarray(ql_mean_comp), np.asarray(cf_comp),
                                 np.asarray(error_ql_env), np.asarray(error_ql_domain),
                                 np.asarray(rel_error_ql_env), np.asarray(rel_error_ql_domain),
                                 np.asarray(error_cf_env), np.asarray(error_cf_domain),
                                 np.asarray(rel_error_cf_env), np.asarray(rel_error_cf_domain),
                                 np.asarray(n_env), np.asarray(n_updraft))
        return






    #----------------------------------------------------------------------
    #----------------------------------------------------------------------

    def dump_error_file(self, path, file_name, time, comp_range, nvar, nz_,
                        ql_mean_env, ql_mean_up, ql_mean_domain,
                        cf_env, cf_up, cf_domain,
                        ql_mean_comp, cf_comp,
                        error_ql_env, error_ql_domain, rel_error_ql_env, rel_error_ql_domain,
                        error_cf_env, error_cf_domain, rel_error_cf_env, rel_error_cf_domain,
                        n_env, n_updraft):
        print('---------- dump error ---------- ')
        print(os.path.join(path, file_name))
        N_comp = max(comp_range)

        rootgrp = nc.Dataset(os.path.join(path, file_name), 'w', format = 'NETCDF4')
        prof_grp = rootgrp.createGroup('profiles')
        prof_grp.createDimension('nz', nz_)
        prof_grp.createDimension('Ncomp', N_comp)
        var = prof_grp.createVariable('zrange', 'f8', ('nz'))
        var[:] = np.asarray(self.zrange)[:]
        var = prof_grp.createVariable('krange', 'f8', ('nz'))
        var[:] = np.asarray(self.krange)[:]

        var = prof_grp.createVariable('n_environment', 'f8', ('nz'))
        var[:] = n_env[:]
        var = prof_grp.createVariable('n_updraft', 'f8', ('nz'))
        var[:] = n_updraft[:]

        var = prof_grp.createVariable('ql_mean_environment', 'f8', ('nz'))
        var[:] = ql_mean_env[:]
        var = prof_grp.createVariable('ql_mean_updraft', 'f8', ('nz'))
        var[:] = ql_mean_up[:]
        var = prof_grp.createVariable('ql_mean_domain', 'f8', ('nz'))
        var[:] = ql_mean_domain[:]
        var = prof_grp.createVariable('cf_environment', 'f8', ('nz'))
        var[:] = cf_env[:]
        var = prof_grp.createVariable('cf_updraft', 'f8', ('nz'))
        var[:] = cf_up[:]
        var = prof_grp.createVariable('cf_domain', 'f8', ('nz'))
        var[:] = cf_domain[:]

        var = prof_grp.createVariable('ql_mean_pdf', 'f8', ('nz', 'Ncomp'))
        var[:] = ql_mean_comp[:,:]
        var = prof_grp.createVariable('cf_pdf', 'f8', ('nz', 'Ncomp'))
        var[:] = cf_comp[:,:]

        error_grp = rootgrp.createGroup('error')
        error_grp.createDimension('nz', nz_)
        error_grp.createDimension('Ncomp', N_comp)
        var = error_grp.createVariable('error_ql_env', 'f8', ('nz', 'Ncomp'))
        var[:,:] = error_ql_env[:,:]
        var = error_grp.createVariable('error_ql_domain', 'f8', ('nz', 'Ncomp'))
        var[:,:] = error_ql_domain[:,:]
        var = error_grp.createVariable('rel_error_ql_env', 'f8', ('nz', 'Ncomp'))
        var[:,:] = rel_error_ql_env[:,:]
        var = error_grp.createVariable('rel_error_ql_domain', 'f8', ('nz', 'Ncomp'))
        var[:,:] = rel_error_ql_domain[:,:]
        var = error_grp.createVariable('error_cf_env', 'f8', ('nz', 'Ncomp'))
        var[:,:] = error_cf_env[:,:]
        var = error_grp.createVariable('error_cf_domain', 'f8', ('nz', 'Ncomp'))
        var[:,:] = error_cf_domain[:,:]
        var = error_grp.createVariable('rel_error_cf_env', 'f8', ('nz', 'Ncomp'))
        var[:,:] = rel_error_cf_env[:,:]
        var = error_grp.createVariable('rel_error_cf_domain', 'f8', ('nz', 'Ncomp'))
        var[:,:] = rel_error_cf_domain[:,:]

        rootgrp.close()
        return


    def create_statistics_file(self, path, file_name, time, ncomp, nvar, nz_):
        print('create statistics file: '+ path+', '+ file_name)
        # ncomp: number of Gaussian components in EM
        # nvar: number of variables of multi-variate Gaussian components
        rootgrp = nc.Dataset(os.path.join(path,file_name), 'w', format='NETCDF4')
        dimgrp = rootgrp.createGroup('dims')
        ts_grp = rootgrp.createGroup('time')
        ts_grp.createDimension('nt',len(time)-1)
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
        weights_grp.createDimension('ncomp', ncomp)
        error_grp = rootgrp.createGroup('error')
        error_grp.createDimension('nz', nz_)

        var = ts_grp.createVariable('t','f8',('nt'))
        for i in range(len(time)-1):
            var[i] = time[i+1]
        z_grp = rootgrp.createGroup('profiles')
        z_grp.createDimension('nz', nz_)
        var = z_grp.createVariable('z', 'f8', ('nz'))
        for i in range(nz_):
            var[i] = self.zrange[i]
        var = z_grp.createVariable('k', 'f8', ('nz'))
        for i in range(nz_):
            var[i] = self.krange[i]
        rootgrp.close()
        return


#----------------------------------------------------------------------
def dump_variable(path, group_name, data_, var_name, ncomp, nvar, nz_):
    print('-------- dump variable --------', var_name, group_name, path)
    # print('dump variable', path, group_name, var_name, data_.shape, ncomp, nvar)
    rootgrp = nc.Dataset(path, 'r+')
    if group_name == 'means':
        # rootgrp = nc.Dataset(path, 'r+')
        var = rootgrp.groups['means'].createVariable(var_name, 'f8', ('nz', 'ncomp', 'nvar'))
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

    elif group_name == 'profiles':
        var = rootgrp.groups['profiles'].createVariable(var_name, 'f8', ('nz'))
        var[:] = data_[:]

    # # write_field(path, group_name, data, var_name)
    # # print('--------')
    rootgrp.close()
    print('')
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





def read_in_updrafts_colleen(type, t, path_):
    print('')
    print('- Updraft Colleen: read in ' + type + ' -')
    import pickle
    print('')


    path = path_
    files = os.listdir(path)
    print(path_)
    print('time: ', t)
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

    # print('')
    # path = os.path.join(path_, root + 'updraft.pkl')
    # data = pickle.load(open(path))
    # print(path + ': ', data.keys())
    #
    # print('')
    # path = os.path.join(path_, root + 'environment.pkl')
    # data = pickle.load(open(path))
    # print(path + ': ', data.keys())

    print('')
    path = os.path.join(path_, root + 'time_'+ str(t) + '_Grid.pkl')
    print(path)
    print('')
    labels = pickle.load(open(path))
    # print(type(labels))
    # print(path + ': ', labels.shape())

    return labels


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