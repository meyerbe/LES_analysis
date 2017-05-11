import os, sys
import argparse
import json as simplejson
import numpy as np
import pylab as plt

sys.path.append("..")
from io_read_in_files import read_in_netcdf

label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 8


def main():
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("path")
    parser.add_argument("casename")
    args = parser.parse_args()
    path = args.path
    global case_name
    case_name = args.casename
    path_ref = os.path.join(path, 'Stats.'+case_name+'.nc')

    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    global nx, ny, nz, dz, dt_stats
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dz = nml['grid']['dz']
    dt_stats = nml['stats_io']['frequency']
    print('Grid: dz='+str(dz))
    print('')

    files = os.listdir(os.path.join(path, 'fields'))
    print('All Files: ', files)
    t = np.int(files[0][0:-3])

    day = np.int(24 * 3600)
    hour = np.int(3600)
    # t = files_[0]
    # files_ = [files[0][0:-3]]
    # files_ = [10800, 12000, 13200, 14400]
    # files_ = [18000, 18900, 19800, 20700]

    # files_ = [0, 1*day, 2*day, 3*day, 4*day, 5*day]
    # files_ = [0, 1 * day, 2 * day, 3 * day, 4 * day, 5 * day, 6*day, 7*day, 8*day]
    # files_ = [6*day, 7*day, 8*day, 9*day, 10*day, 11*day, 12*day, 13*day]
    # files_ = files
    # files_ = [0]

    # Bomex_test
    levels = [400, 500, 600, 720, 800, 900, 1000]
    # files_ = ['21600.nc']
    files_ = ['3600.nc']
    files_cum = files_

    #  Bomex 170314_weno7
    # files_ = [0, 1 * hour, 2 * hour, 3 * hour, 4 * hour, 5 * hour, 6 * hour]
    # files_cum = [5*hour, np.int(5.5*hour), 6*hour]

    #  Rico
    # files_ = [0, 2 * hour, 4 * hour, 6 * hour, 8 * hour, 10 * hour, 12 * hour,
    #                 14 * hour, 16 * hour, 18 * hour, 20 * hour, 22 * hour, 24 * hour]
    # files_cum = [23 * hour, np.int(23.5 * hour), 24 * hour]

    #  DYCOMS RF01
    # files_ = ['1800.nc','3600.nc', '7200.nc']
    # files_ = [0, 1 * hour, 2 * hour, 3 * hour, 4 * hour]
    # files_cum = ['10800.nc', '12600.nc', '14400.nc']

    # DYCOMS RF02
    # files_ = [0, 1 * hour, 2 * hour, 3 * hour, 4 * hour, 5 * hour, 6 * hour]
    # files_cum = ['10800.nc', '12600.nc', '14400.nc']
    # files_cum = ['18000.nc', '19800.nc', '21600.nc']

    # ZGILS S6
    # files_ = files
    # files_cum = files[0:44:2]

    # ZGILS S12
    # files_ = files
    # files_cum = ['86400.nc', '172800.nc', '259200.nc', '345600.nc', '432000.nc', '518400.nc', '604800.nc', '691200.nc']
    # files_cum = ['345600.nc', '432000.nc', '518400.nc', '604800.nc', '691200.nc']

    # TRMM_LBA
    # files_ = files
    # files_cum = ['1012600.nc', '1014400.nc', '1016200.nc']


    # files_ = [1317600, 1339200, 1360800, 1382400]      # ZGILS 6

    print('Selected Files: ', files_)

    # path_fields = os.path.join(path, 'fields', str(t)+'.nc')
    # ql = read_in_netcdf('ql', 'fields', path_fields)
    print('path_ref', path_ref)
    zrange = read_in_netcdf('z', 'reference', path_ref)
    print('')
    print('zrange: ', zrange)

    # plot_ql_n_all(files_, zrange, path, prof=False)
    # plot_mean('ql', files_, zrange, path)
    # plot_mean('s', files_, zrange, path)
    # plot_mean('thetali', files_, zrange, path)
    plot_mean_levels('qt', files_, zrange, levels, path)
    plot_mean_levels('ql', files_, zrange, levels, path)

    # plot_mean('qt', files_, zrange, levels, path)
    # plot_mean('ql', files_, zrange, levels, path)
    # plot_mean_var('ql', files_, zrange, path)
    # plot_mean_var('qt', files_, zrange, path)
    # plot_mean_var('s', files_, zrange, path)
    # plot_mean_var('thetali', files_, zrange, path)

    # plot_mean_cumulated('ql', files_cum, zrange, path)
    # plot_mean_cumulated('s', files_cum, zrange, path)

    # plot_max_var('ql', zrange, path)

    return





def plot_mean_cumulated(var_name, files_cum, zrange, path):
    global case_name
    global nz
    mean_all = np.zeros(nz)
    # path_references = os.path.join(path, 'Stats.' + case_name + '.nc')
    # time = read_in_netcdf('t', 'timeseries', path_references)

    print('')
    print('files_cum', files_cum, len(files_cum))

    plt.figure(figsize=(9, 6))
    cm1 = plt.cm.get_cmap('bone')
    cm2 = plt.cm.get_cmap('winter')
    count_color = 0.0
    day = 24 * 3600

    for t in files_cum:
        if str(t)[-1] == 'c':
            path_fields = os.path.join(path, 'fields', str(t))
            it = np.int(t[0:-3])
        else:
            path_fields = os.path.join(path, 'fields', str(t) + '.nc')
            it = np.double(t)
        var_field = read_in_netcdf(var_name, 'fields', path_fields)
        var_mean_field = np.mean(np.mean(var_field, axis=0), axis=0)
        mean_all += var_mean_field


        if it < 3600:
            lab = str(it) + 's'
            print('sss')
        elif it < 3600 * 24:

            lab = str(np.round(it/3600,0)) + 'h'
            print('hhh')
        else:
            d = it / day
            h = np.mod(it, day) / 3600
            s = np.mod(np.mod(it, day), 3600)
            # lab = str(np.round(d,0)) + 'days ' + str(np.round(h,0)) + 'h ' + str(np.round(s,0)) + 's'
            lab = str(np.round(d, 0)) + 'days ' + str(np.round(h, 1)) + 'h '
        plt.plot(var_mean_field[:], zrange, '--', color=cm1(count_color / len(files_cum)), label='t=' + lab + ' (from field)')
        count_color += 1.0

    mean_all /= len(files_cum)
    plt.plot(mean_all[:], zrange, 'k', label='time mean')
    plt.legend()
    plt.xlabel('mean ' + var_name)
    plt.ylabel('height z [m]')
    plt.title('mean ' + var_name + ' (' + case_name + ', n*nx*ny=' + str(nx * ny * len(files_cum)) + ')')
    plt.savefig(
        os.path.join(path, 'figs_stats', var_name + '_mean_fromfield_cum.pdf'))
    # plt.show()
    plt.close()

    return


def plot_mean_var(var_name, files_, zrange, path):
    global case_name
    global nx, ny, nz, dz, dt_stats
    cm1 = plt.cm.get_cmap('bone')
    cm2 = plt.cm.get_cmap('winter')
    count_color = 0.0

    plt.figure(figsize=(14, 7))
    for t in files_:
        if str(t)[-1] == 'c':
            path_fields = os.path.join(path, 'fields', str(t))
            it = np.int(t[0:-3])
        else:
            path_fields = os.path.join(path, 'fields', str(t) + '.nc')
            it = t

        # path_fields = os.path.join(path, 'fields', str(t) + '.nc')
        var = read_in_netcdf(var_name, 'fields', path_fields)
        var_mean_fields = np.mean(np.mean(var, axis=0), axis=0)
        var2_mean = np.mean(np.mean(var*var, axis=0), axis=0)
        var_variance = var2_mean - var_mean_fields
        plt.subplot(1, 2, 1)
        plt.plot(var_mean_fields, zrange, color = cm2(count_color/len(files_)), label='t=' + str(it) + 's')
        plt.subplot(1, 2, 2)
        plt.plot(var_variance, zrange, color = cm2(count_color/len(files_)), label='t=' + str(it) + 's')
        count_color += 1.0

    plt.subplot(1, 2, 1)
    plt.legend()
    plt.xlabel('mean ' + var_name)
    plt.ylabel('height z [m]')
    plt.title('mean '+var_name+' (' + case_name + ', nx*ny=' + str(nx * ny) + ')')
    plt.subplot(1, 2, 2)
    plt.legend()
    plt.xlabel('Var['+var_name+']')
    plt.ylabel('height z [m]')
    plt.title('Var['+var_name+'] (' + case_name + ', nx*ny=' + str(nx * ny) + ')')
    plt.savefig(
        os.path.join(path, 'figs_stats', var_name+'_mean_var_fromfield.png'))
    # plt.show()
    plt.close()
    return


def plot_max_var(var_name, zrange, path):
    global case_name
    global nx, ny, nz, dz, dt_stats
    cm1 = plt.cm.get_cmap('bone')
    cm2 = plt.cm.get_cmap('winter')
    count_color = 0.0

    path_references = os.path.join(path, 'Stats.' + case_name + '.nc')
    var_max = read_in_netcdf(var_name + '_max', 'profiles', path_references)
    time = read_in_netcdf('t', 'timeseries', path_references)

    plt.figure(figsize=(14, 7))

    for it in range(var_max.shape[0]):
        tt = np.int((np.double(it) - time[0]) / np.double(dt_stats))
        plt.plot(var_max[tt, :], zrange, 'k--', label='t=' + str(tt * dt_stats) + 's')

    plt.legend()
    plt.xlabel('max' + var_name)
    plt.ylabel('height z [m]')
    plt.title('max' + var_name + ' (' + case_name + ', nx*ny=' + str(nx * ny) + ')')
    plt.title('max[' + var_name + '] (' + case_name + ', nx*ny=' + str(nx * ny) + ')')
    plt.savefig(
        os.path.join(path, 'figs_stats', var_name + '_max_fromprofile.png'))
    # plt.show()
    plt.close()
    return


def plot_ql_n_all(files_, zrange, path, prof=True):
    global case_name
    global nx, ny, nz, dz, dt_stats
    cm1 = plt.cm.get_cmap('bone')
    cm2 = plt.cm.get_cmap('winter')

    count_color = 0.0
    plt.figure(figsize=(14,7))
    for t in files_:
        print('')
        if str(t)[-1] == 'c':
            path_fields = os.path.join(path, 'fields', str(t))
            it = np.int(t[0:-3])
        else:
            path_fields = os.path.join(path, 'fields', str(t) + '.nc')
            it = t

        ql = read_in_netcdf('ql', 'fields', path_fields)
        zql = []
        nql = []
        nql_all = []
        for z in range(nz):
            n = np.count_nonzero(ql[:, :, z])
            nql_all.append(n)
            if n > 0:
                zql.append(z)
                nql.append(n)

        ql_mean_fields = np.mean(np.mean(ql,axis=0),axis=0)
        ax = plt.subplot(1,2,1)
        plt.plot(nql_all, zrange, color = cm2(count_color/len(files_)), label='t=' + str(it) + 's')
        if np.any(nql_all) > 0.0:
            try:
                ax.set_xscale('log')
            except:
                print('no log-scaling possible (nql_all=' + str(nql_all) + ')')
                pass
        plt.subplot(1, 2, 2)
        plt.plot(ql_mean_fields, zrange, color=cm2(count_color/len(files_)), label='t=' + str(it) + 's')
        if prof:
            path_references = os.path.join(path, 'Stats.' + case_name + '.nc')
            ql_mean = read_in_netcdf('ql_mean', 'profiles', path_references)
            time = read_in_netcdf('t', 'timeseries', path_references)
            tt = np.int((np.double(it)-time[0]) / np.double(dt_stats))
            plt.plot(ql_mean[tt,:], zrange,  'k--', label='t=' + str(tt*dt_stats) + 's')
        count_color += 1.0

    plt.subplot(1, 2, 1)
    plt.legend()
    plt.xlabel('# non-zero ql')
    plt.ylabel('height z [m]')
    plt.title('# non-zero ql (' + case_name + ', nx*ny=' + str(nx * ny) + ')')
    plt.subplot(1, 2, 2)
    plt.legend()
    plt.xlabel('mean ql')
    plt.ylabel('height z [m]')
    plt.title('mean ql (' + case_name + ', nx*ny=' + str(nx * ny) + ')')

    plt.savefig(
        os.path.join(path, 'figs_stats', 'ql_number_fromfield.png'))
    # plt.show()
    plt.close()
    return



def plot_mean(var_name, files_, zrange, levels, path):
    print('plotting mean')
    global case_name
    global nx, ny, nz, dz, dt_stats
    path_references = os.path.join(path, 'Stats.' + case_name + '.nc')
    var_mean = read_in_netcdf(var_name+'_mean', 'profiles', path_references)
    time = read_in_netcdf('t', 'timeseries', path_references)
    # print(time)
    # cm1 = plt.cm.get_cmap('viridis')
    cm1 = plt.cm.get_cmap('bone')
    cm2 = plt.cm.get_cmap('winter')
    count_color = 0.0

    plt.figure(figsize=(9, 6))
    # print('files'), files_
    for t in files_:
        print('t', t)
        if str(t)[-1] == 'c':

            path_fields = os.path.join(path, 'fields', str(t))
            if case_name == 'TRMM_LBA':
                it = np.int(t[3:-3])
            else:
                it = np.int(t[0:-3])
        else:
            path_fields = os.path.join(path, 'fields', str(t) + '.nc')
            if case_name == 'TRMM_LBA':
                it = np.int(t[3:-1])
            else:
                it = np.double(t)
        # path_fields = os.path.join(path, 'fields', str(t) + '.nc')
        var_field = read_in_netcdf(var_name, 'fields', path_fields)
        var_mean_field = np.mean(np.mean(var_field, axis=0), axis=0)

        tt = np.int((np.double(it)-time[0]) / np.double(dt_stats))
        # print('tt', tt, 'dt_stats', dt_stats, 'it', it, 'time[0]', time[0])
        plt.plot(var_mean_field[:], zrange, color=cm1(count_color/len(files_)), label='t=' + str(it) + 's (from field)')
        plt.plot(var_mean[tt, :], zrange, '--', color=cm2(count_color/len(files_)), label='t=' + str(tt * dt_stats) + 's (from Stats)')
        count_color += 1.0
    mini = np.min([np.amin(var_mean[tt,:]),np.amin(var_mean_field)])
    maxi = np.max([np.amax(var_mean[tt,:]),np.amax(var_mean_field)])
    # if case_name == 'Bomex':
    #     # Bomex 170314
    #     plt.plot([mini, maxi], [800, 800], 'k', linewidth=1, label='800m')
    #     plt.plot([mini, maxi], [1560, 1560], 'k', linewidth=1, label='1520m')
    #     # Bomex Kyle
    #     # plt.plot([mini, maxi], [1080, 1080], 'k', linewidth=1)
    #     # plt.plot([mini, maxi], [3600, 3600], 'k', linewidth=1)
    #     # plt.plot([mini, maxi], [3640, 3640], 'k', linewidth=1)
    # elif case_name[0:8] == 'ZGILS_S12':
    #     # ZGILS 12
    #     # plt.plot([mini, maxi], [700, 700], 'k', linewidth=1, label='700m')
    #     # plt.plot([mini, maxi], [800, 800], 'k', linewidth=1, label='800m')
    #     # plt.plot([mini, maxi], [900, 900], 'k', linewidth=1, label='900m')
    # elif case_name[0:7] == 'ZGILS_S6':
    #     # ZGILS 6
    #     plt.plot([mini, maxi], [700, 700], 'k', linewidth=1, label='700m')
    #     plt.plot([mini, maxi], [1000, 1000], 'k', linewidth=1, label='1000m')
    #     plt.plot([mini, maxi], [1500, 1500], 'k', linewidth=1, label='1500m')
    #     plt.plot([mini, maxi], [2500, 2500], 'k', linewidth=1, label='2500m')
    # elif case_name[0:10] == 'DYCOMS_RF01':
    #     # Dycoms RF01
    #     plt.plot([mini, maxi], [600, 600], 'k', linewidth=1, label='600m')
    #     plt.plot([mini, maxi], [830, 830], 'k', linewidth=1, label='830m')
    #     plt.plot([mini, maxi], [900, 900], 'k', linewidth=1, label='900m')
    # elif case_name[0:10] == 'DYCOMS_RF02':
    #     # DYCOMS RF02
    #     plt.plot([mini, maxi], [600, 600], 'k', linewidth=1, label='600m')
    #     plt.plot([mini, maxi], [850, 850], 'k', linewidth=1, label='850m')
    for l in levels:
        plt.plot([mini, maxi], [l, l], 'k', linewidth=1, label=str(l)+'m')
    plt.legend()
    plt.xlabel('mean '+var_name)
    plt.ylabel('height z [m]')
    plt.title('mean '+var_name + ' (' + case_name + ', nx*ny=' + str(nx * ny) + ')')

    plt.savefig(
        os.path.join(path, 'figs_stats', var_name + '_mean_fromfield.png'))
    # plt.show()
    plt.close()
    return


def plot_mean_levels(var_name, files_, zrange, levels, path):
    print('plotting mean levels')
    global case_name
    global nx, ny, nz, dz, dt_stats
    path_references = os.path.join(path, 'Stats.' + case_name + '.nc')
    var_mean = read_in_netcdf(var_name+'_mean', 'profiles', path_references)
    time = read_in_netcdf('t', 'timeseries', path_references)
    # print(time)
    # cm1 = plt.cm.get_cmap('viridis')
    cm1 = plt.cm.get_cmap('bone')
    cm2 = plt.cm.get_cmap('winter')
    count_color = 0.0

    plt.figure(figsize=(9, 6))
    # print('files'), files_
    for t in files_:
        print('t', t)
        if str(t)[-1] == 'c':
            path_fields = os.path.join(path, 'fields', str(t))
            if case_name == 'TRMM_LBA':
                it = np.int(t[3:-3])
            else:
                it = np.int(t[0:-3])
        else:
            path_fields = os.path.join(path, 'fields', str(t) + '.nc')
            if case_name == 'TRMM_LBA':
                it = np.int(t[3:-1])
            else:
                it = np.double(t)
        # path_fields = os.path.join(path, 'fields', str(t) + '.nc')
        var_field = read_in_netcdf(var_name, 'fields', path_fields)
        var_mean_field = np.mean(np.mean(var_field, axis=0), axis=0)

        tt = np.int((np.double(it)-time[0]) / np.double(dt_stats))
        # print('tt', tt, 'dt_stats', dt_stats, 'it', it, 'time[0]', time[0])
        mini = np.min([np.amin(var_mean[tt, :]), np.amin(var_mean_field)])
        maxi = np.max([np.amax(var_mean[tt, :]), np.amax(var_mean_field)])

        if var_name == 'ql':
            location = 4
            ql = np.ndarray(shape=(0))
            z_ = np.ndarray(shape=(0))
            k_ = np.ndarray(shape=(0), dtype=np.int)
            for k in range(zrange.shape[0]):
                if var_mean_field[k] > 0.0:
                    ql = np.append(ql, var_mean_field[k])
                    z_ = np.append(z_, zrange[k])
                    k_ = np.append(k_, k)
            for l in z_:
                if np.mod(l, 50) == 0:
                    plt.plot([mini, maxi], [l, l], color='0.5', linewidth=1.5, label=str(l) + 'm')
                else:
                    plt.plot([mini, maxi], [l, l], color='0.5', linewidth=0.5)
            plt.plot(ql[:], z_, color=cm1(count_color / len(files_)),
                         label='t=' + str(it) + 's (from field)')
            # plt.plot(var_mean[tt, :], zrange, '--', color=cm2(count_color / len(files_)),
            #              label='t=' + str(tt * dt_stats) + 's (from Stats)')
        else:
            location = 3
            for l in zrange:
                if np.mod(l,100) == 0:
                    plt.plot([mini, maxi], [l, l], color='0.5', linewidth=1.0, label=str(l)+'m')
                elif np.mod(l,10*dz) == 0:
                    plt.plot([mini, maxi], [l, l], color='0.2', linewidth=0.2)
            plt.plot(var_mean_field[:], zrange, color=cm1(count_color/len(files_)), label='t=' + str(it) + 's (from field)')
            plt.plot(var_mean[tt, :], zrange, '--', color=cm2(count_color/len(files_)), label='t=' + str(tt * dt_stats) + 's (from Stats)')
        count_color += 1.0


    plt.legend(loc=location, fontsize=6)
    plt.xlabel('mean '+var_name)
    plt.ylabel('height z [m]')
    plt.title('mean '+var_name + ' (' + case_name + ', nx*ny=' + str(nx * ny) + ', dz='+str(dz)+')')

    plt.savefig(
        os.path.join(path, 'figs_stats', var_name + '_mean_fromfield_levels.png'))
    # plt.show()
    plt.close()
    return


if __name__ == "__main__":
    main()