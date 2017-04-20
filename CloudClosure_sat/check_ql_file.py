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
    t = np.int(files[0][0:-3])

    day = np.int(24 * 3600)
    hour = np.int(3600)
    # t = files_[0]
    files_ = [files[0][0:-3]]
    # files_ = [10800, 12000, 13200, 14400]
    # files_ = [18000, 18900, 19800, 20700]

    # files_ = [0, 1*day, 2*day, 3*day, 4*day, 5*day]
    files_ = [0, 1 * day, 2 * day, 3 * day, 4 * day, 5 * day, 6*day, 7*day, 8*day]
    # files_ = [6*day, 7*day, 8*day, 9*day, 10*day, 11*day, 12*day, 13*day]
    # files_ = files
    # files_ = [0]

    #  Bomex 170314_weno7
    files_ = [0, 1 * hour, 2 * hour, 3 * hour, 4 * hour, 5 * hour, 6 * hour]
    files_cum = [5*hour, np.int(5.5*hour), 6*hour]

    #  Rico
    # files_ = [0, 2 * hour, 4 * hour, 6 * hour, 8 * hour, 10 * hour, 12 * hour,
    #                 14 * hour, 16 * hour, 18 * hour, 20 * hour, 22 * hour, 24 * hour]

    #  DYCOMS RF01
    # files_ = ['1800.nc','3600.nc']

    # DYCOMS RF02
    # files_ = [0, 1 * hour, 2 * hour, 3 * hour, 4 * hour, 5 * hour, 6 * hour]

    # files_ = [1317600, 1339200, 1360800, 1382400]      # ZGILS 6
    print('All Files: ', files)
    print('Selected Files: ', files_)

    # path_fields = os.path.join(path, 'fields', str(t)+'.nc')
    # ql = read_in_netcdf('ql', 'fields', path_fields)
    zrange = read_in_netcdf('z', 'reference', path_ref)

    # plot_ql_n_all(files_, zrange, path, prof=True)
    # plot_mean('ql', files_, zrange, path)
    # plot_mean('s', files_, zrange, path)
    # plot_mean('thetali', files_, zrange, path)
    # plot_mean('qt', files_, zrange, path)
    # plot_mean_var('ql', files_, zrange, path)
    # plot_mean_var('qt', files_, zrange, path)
    # plot_mean_var('s', files_, zrange, path)
    # plot_mean_var('thetali', files_, zrange, path)

    plot_mean_cumulated('ql', files_cum, zrange, path)

    return





def plot_mean_cumulated(var_name, files_cum, zrange, path):
    global case_name
    global nz
    mean_all = np.zeros(nz)
    # path_references = os.path.join(path, 'Stats.' + case_name + '.nc')
    # time = read_in_netcdf('t', 'timeseries', path_references)

    print('')
    print(files_cum, len(files_cum))

    plt.figure(figsize=(9, 6))
    cm1 = plt.cm.get_cmap('bone')
    cm2 = plt.cm.get_cmap('winter')
    count_color = 0.0

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

        plt.plot(var_mean_field[:], zrange, '--', color=cm1(count_color / len(files_cum)), label='t=' + str(it) + 's (from field)')
        count_color += 1.0

    mean_all /= len(files_cum)
    plt.plot(mean_all[:], zrange, 'k', label='time mean')
    plt.legend()
    plt.xlabel('mean ' + var_name)
    plt.ylabel('height z [m]')
    plt.title('mean ' + var_name + ' (' + case_name + ', nx*ny=1' + str(nx * ny) + ')')
    plt.savefig(
        os.path.join(path, 'figs_stats', var_name + '_mean_fromfield_cum.png'))
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
    plt.title('mean '+var_name+' (' + case_name + ', nx*ny=1' + str(nx * ny) + ')')
    plt.subplot(1, 2, 2)
    plt.legend()
    plt.xlabel('Var['+var_name+']')
    plt.ylabel('height z [m]')
    plt.title('Var['+var_name+'] (' + case_name + ', nx*ny=1' + str(nx * ny) + ')')
    plt.savefig(
        os.path.join(path, 'figs_stats', var_name+'_mean_var_fromfield.png'))
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
        try:
            ax.set_xscale('log')
        except:
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
    plt.title('# non-zero ql (' + case_name + ', nx*ny=1' + str(nx * ny) + ')')
    plt.subplot(1, 2, 2)
    plt.legend()
    plt.xlabel('mean ql')
    plt.ylabel('height z [m]')
    plt.title('mean ql (' + case_name + ', nx*ny=1' + str(nx * ny) + ')')

    plt.savefig(
        os.path.join(path, 'figs_stats', 'ql_number_fromfield.png'))
    # plt.show()
    plt.close()
    return



def plot_mean(var_name, files_, zrange, path):
    global case_name
    global nx, ny, nz, dz, dt_stats
    path_references = os.path.join(path, 'Stats.' + case_name + '.nc')
    var_mean = read_in_netcdf(var_name+'_mean', 'profiles', path_references)
    time = read_in_netcdf('t', 'timeseries', path_references)
    print(time)
    # cm1 = plt.cm.get_cmap('viridis')
    cm1 = plt.cm.get_cmap('bone')
    cm2 = plt.cm.get_cmap('winter')
    count_color = 0.0

    plt.figure(figsize=(9, 6))
    # print('files'), files_
    for t in files_:
        if str(t)[-1] == 'c':
            path_fields = os.path.join(path, 'fields', str(t))
            it = np.int(t[0:-3])
        else:
            path_fields = os.path.join(path, 'fields', str(t) + '.nc')
            it = np.double(t)
        # path_fields = os.path.join(path, 'fields', str(t) + '.nc')
        var_field = read_in_netcdf(var_name, 'fields', path_fields)
        var_mean_field = np.mean(np.mean(var_field, axis=0), axis=0)

        tt = np.int((np.double(it)-time[0]) / np.double(dt_stats))
        print('tt', tt, 'dt_stats', dt_stats)
        plt.plot(var_mean_field[:], zrange, color=cm1(count_color/len(files_)), label='t=' + str(it) + 's (from field)')
        plt.plot(var_mean[tt, :], zrange, '--', color=cm2(count_color/len(files_)), label='t=' + str(tt * dt_stats) + 's (from Stats)')
        count_color += 1.0
    mini = np.min([np.amin(var_mean[tt,:]),np.amin(var_mean_field)])
    maxi = np.max([np.amax(var_mean[tt,:]),np.amax(var_mean_field)])
    # Bomex 170314
    # plt.plot([mini, maxi], [800, 800], 'k', linewidth=1, label='800m')
    # plt.plot([mini, maxi], [1560, 1560], 'k', linewidth=1, label='1520m')
    # Bomex Kyle
    # plt.plot([mini, maxi], [1080, 1080], 'k', linewidth=1)
    # plt.plot([mini, maxi], [3600, 3600], 'k', linewidth=1)
    # plt.plot([mini, maxi], [3640, 3640], 'k', linewidth=1)
    # ZGILS 12
    # plt.plot([mini, maxi], [700, 700], 'k', linewidth=1, label='700m')
    # plt.plot([mini, maxi], [800, 800], 'k', linewidth=1, label='800m')
    # plt.plot([mini, maxi], [900, 900], 'k', linewidth=1, label='900m')
    # ZGILS 6
    # plt.plot([mini, maxi], [700, 700], 'k', linewidth=1, label='700m')
    # plt.plot([mini, maxi], [1000, 1000], 'k', linewidth=1, label='1000m')
    # plt.plot([mini, maxi], [1500, 1500], 'k', linewidth=1, label='1500m')
    # plt.plot([mini, maxi], [2500, 2500], 'k', linewidth=1, label='2500m')
    # Dycoms RF01
    # plt.plot([mini, maxi], [600, 600], 'k', linewidth=1, label='600m')
    # plt.plot([mini, maxi], [830, 830], 'k', linewidth=1, label='830m')
    # plt.plot([mini, maxi], [900, 900], 'k', linewidth=1, label='900m')
    # DYCOMS RF02
    plt.plot([mini, maxi], [600, 600], 'k', linewidth=1, label='600m')
    plt.plot([mini, maxi], [850, 850], 'k', linewidth=1, label='850m')
    plt.legend()
    plt.xlabel('mean '+var_name)
    plt.ylabel('height z [m]')
    plt.title('mean '+var_name + ' (' + case_name + ', nx*ny=1' + str(nx * ny) + ')')

    plt.savefig(
        os.path.join(path, 'figs_stats', var_name + '_mean_fromfield.png'))
    # plt.show()
    plt.close()
    return



if __name__ == "__main__":
    main()