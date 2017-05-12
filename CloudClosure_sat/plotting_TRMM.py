import numpy as np
import pylab as plt
import os
import sys

sys.path.append("..")
from parameters import *
from io_read_in_files import read_in_netcdf
from thermodynamic_functions import CC_Magnus, theta_c

label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['figure.titlesize'] = 15

# from check_ql_file import plot_mean_profile


def main():
    case_name = 'TRMM_LBA'

    nml = {}
    nml['grid'] = {}
    nml['grid']['nx'] = 20
    nml['grid']['ny'] = 20
    nml['grid']['nz'] = 200

    InitTRMM_LBA(nml)

    path = '/Volumes/Data/ClimatePhysics/LES/output/TRMM_LBA'
    path_ref = os.path.join(path, 'Stats.'+case_name+'.nc')
    # out_path = '/Volumes'

    global dt_stats, k_hours, zrange
    dt_stats, k_hours, time_prof = get_dt_stats(path, path_ref)
    zrange = get_zrange(path, path_ref)
    plot_profile('qt', time_prof, path_ref, path)
    plot_profile('ql', time_prof, path_ref, path)
    plot_profile('qi', time_prof, path_ref, path)
    return


# def InitTRMM_LBA(namelist, Grid.Grid Gr,PrognosticVariables.PrognosticVariables PV,
#                        ReferenceState.ReferenceState RS, Th, NetCDFIO_Stats NS, ParallelMPI.ParallelMPI Pa , LatentHeat LH):
def InitTRMM_LBA(nml):
    # reference_profiles = AdjustedMoistAdiabat(namelist, LH, Pa)
    #
    # RS.Tg  = 296.85   # surface values for reference state (RS) which outputs p0 rho0 alpha0
    # RS.Pg  = 991.3*100
    # pvg = Th.get_pv_star(RS.Tg)
    # RS.qtg = eps_v * pvg/(RS.Pg - pvg)


    ## TRMM_LBA inputs

    z_in = np.array([0.130,  0.464,  0.573,  1.100,  1.653,  2.216,  2.760,
                     3.297,  3.824,  4.327,  4.787,  5.242,  5.686,  6.131,
                     6.578,  6.996,  7.431,  7.881,  8.300,  8.718,  9.149,
                     9.611, 10.084, 10.573, 11.008, 11.460, 11.966, 12.472,
                    12.971, 13.478, 13.971, 14.443, 14.956, 15.458, 16.019,
                    16.491, 16.961, 17.442, 17.934, 18.397, 18.851, 19.331,
                    19.809, 20.321, 20.813, 21.329, 30.000]) * 1000 - 130.0 #LES z is in meters

    p_in = np.array([991.3, 954.2, 942.0, 886.9, 831.5, 778.9, 729.8,
                     684.0, 641.7, 603.2, 570.1, 538.6, 509.1, 480.4,
                     454.0, 429.6, 405.7, 382.5, 361.1, 340.9, 321.2,
                     301.2, 281.8, 263.1, 246.1, 230.1, 213.2, 197.0,
                     182.3, 167.9, 154.9, 143.0, 131.1, 119.7, 108.9,
                     100.1,  92.1,  84.6,  77.5,  71.4,  65.9,  60.7,
                      55.9,  51.3,  47.2,  43.3,  10.3]) * 100 # LES pres is in pasc

    T_in = np.array([23.70,  23.30,  22.57,  19.90,  16.91,  14.09,  11.13,
                      8.29,   5.38,   2.29,  -0.66,  -3.02,  -5.28,  -7.42,
                    -10.34, -12.69, -15.70, -19.21, -21.81, -24.73, -27.76,
                    -30.93, -34.62, -38.58, -42.30, -46.07, -50.03, -54.67,
                    -59.16, -63.60, -67.68, -70.77, -74.41, -77.51, -80.64,
                    -80.69, -80.00, -81.38, -81.17, -78.32, -74.77, -74.52,
                    -72.62, -70.87, -69.19, -66.90, -66.90]) + 273.15 # LES T is in deg K

    RH_in = np.array([98.00,  86.00,  88.56,  87.44,  86.67,  83.67,  79.56,
                      84.78,  84.78,  89.33,  94.33,  92.00,  85.22,  77.33,
                      80.11,  66.11,  72.11,  72.67,  52.22,  54.67,  51.00,
                      43.78,  40.56,  43.11,  54.78,  46.11,  42.33,  43.22,
                      45.33,  39.78,  33.78,  28.78,  24.67,  20.67,  17.67,
                      17.11,  16.22,  14.22,  13.00,  13.00,  12.22,   9.56,
                       7.78,   5.89,   4.33,   3.00,   3.00])

    u_in = np.array([0.00,   0.81,   1.17,   3.44,   3.53,   3.88,   4.09,
                     3.97,   1.22,   0.16,  -1.22,  -1.72,  -2.77,  -2.65,
                    -0.64,  -0.07,  -1.90,  -2.70,  -2.99,  -3.66,  -5.05,
                    -6.64,  -4.74,  -5.30,  -6.07,  -4.26,  -7.52,  -8.88,
                    -9.00,  -7.77,  -5.37,  -3.88,  -1.15,  -2.36,  -9.20,
                    -8.01,  -5.68,  -8.83, -14.51, -15.55, -15.36, -17.67,
                   -17.82, -18.94, -15.92, -15.32, -15.32])

    v_in = np.array([-0.40,  -3.51,  -3.88,  -4.77,  -5.28,  -5.85,  -5.60,
                     -2.67,  -1.47,   0.57,   0.89,  -0.08,   1.11,   2.15,
                      3.12,   3.22,   3.34,   1.91,   1.15,   1.01,  -0.57,
                     -0.67,   0.31,   2.97,   2.32,   2.66,   4.79,   3.40,
                      3.14,   3.93,   7.57,   2.58,   2.50,   6.44,   6.84,
                      0.19,  -2.20,  -3.60,   0.56,   6.68,   9.41,   7.03,
                      5.32,   1.14,  -0.65,   5.27,   5.27])

    print('')
    print('nz='+str(v_in.shape), nml['grid']['nz'])
    print(z_in.shape)
    print(T_in.shape)

    # compute qt, ql from RH and reference density profile?!
    # RH = pv / pv*
    # specific humidity: q_w = \rho_w / \rho_t = \rho_w / (\rho_w + \rho_d)
    # mixing ratio: r_w = \rho_w / \rho_d

    # qv_star = pv_star * epsi / (
    #       p[k] - pv_star + epsi * pv_star * RH[k] / 100.0)  # eq. 37 in pressel et al and the def of RH
    # qt[k] = qv_star * RH[k] / 100.0

    # --
    Th = theta_c(p_in, T_in)

    pv_star = CC_Magnus(T_in)
    pv = RH_in * 0.01 * pv_star

    rho_v_star = pv_star / (Rv * T_in)
    rho_v = RH_in * rho_v_star

    qv_star = pv_star * eps_vi / (
        p_in - pv_star + eps_vi * pv_star * RH_in * 0.01)  # eq. 37 in pressel et al and the def of RH
    qt = qv_star * RH_in * 0.01
    # ql =
    # T_comp_thl[i], ql_comp_thl[i], alpha_comp_thl[i] = sat_adj_fromthetali(p_ref[iz], Th_l[i, 0], Th_l[i, 1], CC, LH)
    # theta_l = thetali_c(p, T, qt, )
    # theta_l[ij, k] = thetali_c(p_ref[iz], T_[i, j, iz], qt_[i, j, iz], ql_[i, j, iz], qi_, Lv)

    # plotting
    plt.figure(figsize=(15,15))
    plt.subplot(3, 3, 1)
    plt.plot(T_in, z_in)
    plt.xlabel('T [K]')
    plt.ylabel('z')
    plt.subplot(3, 3, 2)
    plt.plot(T_in[0:6], z_in[0:6])
    plt.xlabel('T [K]')
    plt.ylabel('z')
    plt.subplot(3, 3, 3)
    plt.plot(Th[0:6], z_in[0:6])
    plt.xlabel(r'$\theta$ [K]')
    plt.ylabel('z')


    plt.subplot(3, 3, 4)
    plt.plot(RH_in, z_in)
    plt.xlabel('RH')
    plt.ylabel('z')
    plt.subplot(3, 3, 5)
    plt.plot(rho_v, z_in)
    plt.xlabel('rho_v')
    plt.ylabel('z')

    plt.subplot(3, 3, 7)
    plt.plot(qt, z_in)
    plt.xlabel('qt')
    plt.ylabel('z')
    plt.subplot(3, 3, 8)
    plt.plot(qt[0:6], z_in[0:6])
    plt.xlabel('qt')
    plt.ylabel('z')
    # plt.subplot(3, 3, 9)
    # # plt.plot(rho_v, z_in)
    # plt.xlabel('')
    # plt.ylabel('z')
    plt.suptitle('Initial Thermodynamic Soundings')
    save_name = 'initial_thermod'
    plt.savefig(os.path.join('../TRMM_figs_stats', save_name + '.pdf'))
    plt.close()

    plt.figure(figsize=(15, 8))
    plt.subplot(1, 3, 1)
    plt.plot(p_in, z_in, '-', label='p')
    plt.plot(p_in - pv, z_in, '--', label='p_d')
    plt.legend()
    plt.xlabel('p')
    plt.ylabel('z')
    plt.subplot(1, 3, 2)
    plt.plot(pv_star, z_in, label='pv_star')
    plt.legend()
    plt.xlabel('pv* (CC Magnus)')
    plt.ylabel('z')
    plt.subplot(1, 3, 3)
    plt.plot(pv, z_in, '--', label='pv')
    # plt.plot(p_in, z_in, '--', label='p')
    plt.legend()
    plt.xlabel('pv')
    plt.ylabel('z')
    # plt.show()
    plt.suptitle('Initial Pressure Soundings')
    save_name = 'initial_pressure'
    plt.savefig(os.path.join('../TRMM_figs_stats', save_name + '.pdf'))
    plt.close()

    plt.figure(figsize=(9, 6))
    plt.subplot(1, 3, 1)
    plt.plot(u_in, z_in)
    plt.xlabel('u')
    plt.ylabel('z')
    plt.subplot(1, 3, 2)
    plt.plot(v_in, z_in)
    plt.xlabel('v')
    plt.ylabel('z')
    # plt.show()
    plt.suptitle('Initial Velocity Soundings')
    save_name = 'initial_uv'

    plt.savefig(os.path.join('../TRMM_figs_stats', save_name+'.pdf'))
    plt.close()

    return



def interpolate(nml):

    nz = nml['grid']['nz']

    # cdef:
    #     Py_ssize_t i, j, k, ijk, ishift, jshift
    #     Py_ssize_t istride = Gr.dims.nlg[1] * Gr.dims.nlg[2]
    #     Py_ssize_t jstride = Gr.dims.nlg[2]
    #     Py_ssize_t u_varshift = PV.get_varshift(Gr, 'u')
    #     Py_ssize_t v_varshift = PV.get_varshift(Gr, 'v')
    #     Py_ssize_t w_varshift = PV.get_varshift(Gr, 'w')
    #     Py_ssize_t s_varshift = PV.get_varshift(Gr, 's')
    #     Py_ssize_t qt_varshift = PV.get_varshift(Gr, 'qt')
    #     double[:] T = np.zeros((Gr.dims.nlg[2],), dtype=np.double,order='c')  # change to temp interp to zp_hlaf (LES press is in pasc)
    #     double[:] qt = np.zeros((Gr.dims.nlg[2],), dtype=np.double, order='c')
    #     double[:] u = np.zeros((Gr.dims.nlg[2],), dtype=np.double, order='c')
    #     double[:] v = np.zeros((Gr.dims.nlg[2],), dtype=np.double, order='c')

    T = np.zeros((nz,), dtype=np.double,order='c')  # change to temp interp to zp_hlaf (LES press is in pasc)
    qt = np.zeros((nz,), dtype=np.double,order='c')
    u = np.zeros((nz,), dtype=np.double,order='c')
    v = np.zeros((nz,), dtype=np.double,order='c')

    pv_star = np.zeros((nz,), dtype=np.double, order='c')

    T = np.interp(Gr.zp_half, z_in, T_in)
    p = np.interp(Gr.zp_half, z_in, p_in)
    RH = np.interp(Gr.zp_half, z_in, RH_in)
    u = np.interp(Gr.zp_half, z_in, u_in)
    v = np.interp(Gr.zp_half, z_in, v_in)

    # Set velocities for Galilean transformation
    RS.u0 = 0.5 * (np.amax(u) + np.amin(u))
    RS.v0 = 0.5 * (np.amax(v) + np.amin(v))

    # Generate initial perturbations (here we are generating more than we need)
    # random fluctuations
    # I need to perturbed the temperature and only later calculate the entropy
    np.random.seed(Pa.rank)
    # cdef:
    #     double[:] T_pert = np.random.random_sample(Gr.dims.npg)
    #     double T_pert_
    #     double pv_star
    #     double qv_star

    epsi = 287.1 / 461.5
    # Here we fill in the 3D arrays
    # We perform saturation adjustment on the S6 data, although this should not actually be necessary (but doesn't hurt)
    for k in xrange(nz):
        pv_star[k] = CC_Magnus(T[k])
        qv_star[k] = pv_star[k] * epsi / (p[k] - pv_star[k] + epsi * pv_star[k] * RH[k] / 100.0)
    # for i in xrange(Gr.dims.nlg[0]):
    #     ishift = istride * i
    #     for j in xrange(Gr.dims.nlg[1]):
    #         jshift = jstride * j
    #         for k in xrange(Gr.dims.nlg[2]):
    #             ijk = ishift + jshift + k
    #             pv_star = Th.get_pv_star(T[k])
    #             qv_star = pv_star * epsi / (
    #             p[k] - pv_star + epsi * pv_star * RH[k] / 100.0)  # eq. 37 in pressel et al and the def of RH
    #             qt[k] = qv_star * RH[k] / 100.0
    #             PV.values[ijk + u_varshift] = u[k] - RS.u0
    #             PV.values[ijk + v_varshift] = v[k] - RS.v0
    #             PV.values[ijk + w_varshift] = 0.0
    #             PV.values[ijk + qt_varshift] = qt[k]
    #
    #             if Gr.zp_half[k] < 1000.0:
    #                 T_pert_ = (T_pert[ijk] - 0.5) * 0.1
    #                 PV.values[ijk + s_varshift] = Th.entropy(RS.p0_half[k], T[k] + T_pert_, qt[k], 0.0, 0.0)
    #             else:
    #                 PV.values[ijk + s_varshift] = Th.entropy(RS.p0_half[k], T[k], qt[k], 0.0, 0.0)

    return



def get_dt_stats(path, path_ref):
    time_prof = read_in_netcdf('t', 'timeseries', path_ref)
    dt_stats = time_prof[1] - time_prof[0]

    k_hours = np.ndarray(shape=(0), dtype=np.int)
    for t in range(time_prof.shape[0]):
        if np.abs(np.mod(time_prof[t], 3600)) < dt_stats:
            k_hours = np.append(k_hours, np.int(t))

    plt.figure()
    plt.plot(time_prof, '-x')
    plt.title('stats time (dt = ' + str(dt_stats) + 's)')
    plt.savefig(os.path.join(path,'figs_stats','time_stats.pdf'))
    plt.close()
    return dt_stats, k_hours, time_prof


def get_zrange(path, path_ref):
    zrange = read_in_netcdf('z', 'reference', path_ref)
    plt.figure()
    plt.plot(zrange, '-x')
    plt.title('z-profile')
    plt.savefig(os.path.join(path,'figs_stats','z_stats.pdf'))
    plt.close()
    return zrange


def plot_profile(var_name, time, path_ref, path):
    global dt_stats, k_hours, zrange

    var_mean = read_in_netcdf(var_name+'_mean', 'profiles', path_ref)
    var_mean2 = read_in_netcdf(var_name + '_mean2', 'profiles', path_ref)
    var_variance = var_mean2 - var_mean**2
    var_max = read_in_netcdf(var_name+'_max', 'profiles', path_ref)

    plt.figure(figsize=(16,6))
    plt.subplot(1, 3, 1)
    for k in k_hours:
        plt.plot(var_mean[k,:], zrange, label='t='+str(time[k]/3600)+'h')
    plt.xlabel('E[' + var_name + ']')
    plt.ylabel('height z [m]')
    plt.legend(loc='best')
    plt.title('E['+var_name+']')
    plt.subplot(1, 3, 2)
    for k in k_hours:
        plt.plot(var_variance[k, :], zrange, label='t='+str(time[k]/3600)+'h')
    plt.xlabel('Var['+var_name+']')
    plt.ylabel('height z [m]')
    plt.legend(loc='best')
    plt.title('Var[' + var_name + ']')
    plt.subplot(1, 3, 3)
    for k in k_hours:
        plt.plot(var_max[k, :], zrange, label='t=' + str(time[k] / 3600) + 'h')
    plt.xlabel('max(' + var_name + '[:,:,k])')
    plt.ylabel('height z [m]')
    plt.legend(loc='best')
    plt.title('max(' + var_name + ')')
    plt.savefig(os.path.join(path,'figs_stats', var_name+'_mean_var.pdf'))
    return


if __name__ == '__main__':
    main()

