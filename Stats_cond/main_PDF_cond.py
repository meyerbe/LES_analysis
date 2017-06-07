import os
import argparse
import json as simplejson
import numpy as np
import pylab as plt
import netCDF4 as nc

label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['xtick.direction']='out'
plt.rcParams['ytick.direction']='out'
plt.rcParams['legend.fontsize'] = 'small'
plt.rcParams['figure.titlesize'] = 15
plt.rcParams['lines.linewidth'] = 0.8
plt.rcParams['image.cmap'] = 'viridis'



import PDF_conditional

def main():
    Lx = 5e3
    Ly = 5e3
    dk_range = [0, 2, 4]
    # dk_range = [2, 4]

    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("path")
    parser.add_argument("casename")
    parser.add_argument("--time")
    # parser.add_argument("--Lx")
    parser.add_argument('--dk', nargs='+', type=int)
    # # This is the correct way to handle accepting multiple arguments (or list)
    # # '+' == 1 or more.
    # # '*' == 0 or more.
    # # '?' == 0 or 1.
    # # An int is an explicit number of arguments to accept.
    # parser.add_argument('--nargs', nargs='+')
    # # To make the input integers
    # parser.add_argument('--nargs-int-type', nargs='+', type=int)

    args = parser.parse_args()
    # path = '/Volumes/Data/ClimatePhysics/LES/updrafts_colleen/'
    # case_name = 'Bomex'
    time = 21600
    # if args.path:
    #     path = args.path
    # if args.casename:
    #     case_name = args.casename
    path = args.path
    case_name = args.casename
    if args.time:
        time = np.int(args.time)
    # dk_range: number of layers below and above added to data
    if args.dk:
        dk_range = args.dk

    path_ref = os.path.join(path, 'Stats.' + case_name + '.nc')
    path_fields = os.path.join(path, 'fields')

    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = nml['grid']['dx']
    dy = nml['grid']['dy']
    dz = nml['grid']['dz']

    PDF_cond = PDF_conditional.PDF_conditional()

    files = os.listdir(os.path.join(path, 'fields'))
    print('Found the follwing files: ' + str(files))
    ncomp_range = [1, 2, 3, 4, 5, 6, 7, 8]
    # ncomp_range = [1]
    krange, files = set_zrange(case_name)
    N = len(files)
    print('Use the following files' + str(files) + ', ' +str(N))
    print('')
    print('krange' + ', ' + str(type(krange)) + ', ' + str(type(krange[0])))
    print('zrange: ' + str(krange * dz))
    print('dkrange: ' + str(dk_range) + ', ' + str(type(dk_range)) + ', ' + str(type(dk_range[0])))
    print('Lx, Ly:' + ', ' + str(Lx) + ', ' + str(Ly) + ', nx*dx: ' + str(nx * dx))
    print('dz: ' + str(dz))
    print('ncomp: ' + str(ncomp_range))
    print('')
    print(path)
    print('')

    PDF_cond.initialize(krange, path, case_name)
    n_sample = 1e6

    # for Lx in [1000, 5000, 10000, 20000]:
    for Lx in [5000]:
        Ly = Lx
        for dk in dk_range:
            PDF_cond.predict_pdf(files, path, n_sample, ncomp_range, Lx, Ly, dk, krange, nml)

    return


# ----------------------------------
def set_zrange(case_name):
    if case_name[0:8] == 'ZGILS_S6':
        files = ['1382400.nc']
        # files_ = [1317600, 1339200, 1360800, 1382400]  # ZGILS 6
        # files = files[0:25:2]
        krange = np.asarray([25, 25, 40, 50, 60, 65], dtype=np.int32)
    elif case_name[0:9] == 'ZGILS_S12':
        # files = ['432000.nc']
        files = ['86400.nc']
        # files = ['345600.nc', '432000.nc', '518400.nc', '604800.nc', '691200.nc']
        krange = np.asarray([35,40,45])
    elif case_name == 'DYCOMS_RF01':
        # DYCOMS RF01 large
        krange = np.asarray([135, 140, 145, 150, 155, 160, 166, 180], dtype=np.int32)
        # files = ['3600.nc']
        # DYCOMS RF01
        # krange = np.asarray([140, 150, 160, 166, 180], dtype=np.int32)
        files = ['10800.nc', '14400.nc']
        # files = ['10800.nc']
    elif case_name == 'DYCOMS_RF02':
        # krange = np.asarray([140,150], dtype=np.int32)
        # krange = np.asarray(150, dtype=np.int32)
        krange = np.asarray([120, 140, 150, 160, 170, 200], dtype=np.int32)
        # files = ['18000.nc', '19800.nc', '21600.nc']
        files = ['14400.nc']
    elif case_name == 'Bomex':
        ## Bomex large, kyle
        ## Bomex 170314_weno7 (dz=40)
        # krange = np.asarray([10, 12, 15, 18, 20, 22, 25, 40, 50], dtype=np.int32)
        # files = ['21600.nc']
        ## Bomex (dz=20)
        # krange = np.asarray([50, 60], dtype=np.int32)
        krange = np.asarray([15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90, 100, 125], dtype=np.int32)
        files = ['18000.nc', '19800.nc', '21600.nc']
        # files = ['21600.nc']
        # Bomex test
        # files = ['21600.nc']
        # krange = np.asarray([10, 17, 20, 25, 50])
        # krange = np.asarray([10, 12, 15, 18, 20, 22, 25, 40, 50])
        # krange = np.asarray([20, 50])
        # krange = np.asarray([18,30,38])
    elif case_name == 'TRMM_LBA':
        # TRMM
        # files = ['1012600.nc', '1014400.nc', '1016200.nc']
        files = ['1014400.nc']
        krange = np.asarray([10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60], dtype=np.int32)

    return krange, files





if __name__ == "__main__":
    main()