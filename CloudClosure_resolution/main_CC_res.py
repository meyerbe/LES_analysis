import sys, os
import argparse
import json as simplejson
import numpy as np

import CloudClosure_res, CloudClosure_res_acc

def main():
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("path")
    parser.add_argument("casename")
    args = parser.parse_args()
    path = args.path
    case_name = args.casename

    path_ref = os.path.join(path, 'Stats.'+case_name+'.nc')
    path_fields = os.path.join(path, 'fields')

    nml = simplejson.loads(open(os.path.join(path, case_name+'.in')).read())
    nz = nml['grid']['nz']
    dz = nml['grid']['dz']
    nx = nml['grid']['nx']
    dx = nml['grid']['dx']

    ClCl = CloudClosure_res_acc.CloudClosure()


    files = os.listdir(os.path.join(path, 'fields'))
    print('Found the follwing files: ' + str(files))
    ncomp_range = [1, 2, 3, 4, 5, 6, 7, 8]
    # ncomp_range = [1]
    krange, files = set_zrange(case_name)
    N = len(files)
    print('Use the following files', files, N)
    print('zrange: ' + str(krange * dz))
    print('dz: ' + str(dz))
    print('ncomp: ' + str(ncomp_range))
    print('')
    print(path)
    print('')

    ClCl.initialize(krange, path, case_name)
    n_sample = 1e6
    Lx = 20e3
    Ly = 20e3
    print('Lx, Ly', Lx, Ly, nx*dx)
    print('krange', type(krange), type(krange[0]))

    # dk_range: number of layers below and above added to data
    # for dk_range in [0, 1, 2, 3, 4]:
    for dk_range in [0]:
        ClCl.predict_pdf(files, path, n_sample, ncomp_range, Lx, Ly, dk_range, krange, nml)

    return


# ----------------------------------
def set_zrange(case_name):
    if case_name[0:8] == 'ZGILS_S6':
        # ZGILS 6
        files = ['1382400.nc']
        # files_ = [1317600, 1339200, 1360800, 1382400]  # ZGILS 6
        # files = files[0:25:2]
        krange = np.asarray([25, 25, 40, 50, 60, 65], dtype=np.int32)
    elif case_name[0:9] == 'ZGILS_S12':
        # ZGILS S12
        # files = ['432000.nc']
        files = ['86400.nc']
        # files = ['345600.nc', '432000.nc', '518400.nc', '604800.nc', '691200.nc']
        krange = np.asarray([35,40,45])
    elif case_name == 'DYCOMS_RF01':
        # DYCOMS RF01 large
        krange = np.asarray([140, 150, 160, 166, 180], dtype=np.int32)
        # files = ['3600.nc']
        # DYCOMS RF01
        krange = np.asarray([140, 150, 160, 166, 180], dtype=np.int32)
        # files = ['10800.nc', '12600.nc', '14400.nc']
        files = ['10800.nc']
    elif case_name == 'DYCOMS_RF02':
        # DYCOMS RF02
        # krange = np.asarray([140,150], dtype=np.int32)
        # krange = np.asarray(150, dtype=np.int32)
        krange = np.asarray([120, 140, 150, 160, 170, 200], dtype=np.int32)
        # files = ['18000.nc', '19800.nc', '21600.nc']
        files = ['7200.nc']
        # files = ['3600.nc']
    elif case_name == 'Bomex':
        # Bomex large, kyle
        # krange = np.asarray([27, 91])
        # Bomex 170314_weno7
        ## krange = np.asarray([15, 20, 25, 30, 35, 40, 45])
        krange = np.asarray([10, 12, 15, 18, 20, 22, 25, 40, 50], dtype=np.int32)
        ## files = ['18000.nc', '19800.nc', '21600.nc']
        files = ['21600.nc']
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
        krange = np.asarray([10, 15, 20, 30, 40, 50, 60], dtype=np.int32)

    return krange, files





if __name__ == "__main__":
    main()