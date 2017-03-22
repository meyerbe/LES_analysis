import sys, os
import argparse
import json as simplejson
import numpy as np

import CloudClosure

# from CloudClosure import do_everything_with_pycles
# from CloudClosure import verification_CC

def main():
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("inpath")
    parser.add_argument("outpath")
    parser.add_argument("casename")
    args = parser.parse_args()
    path_in = args.inpath
    path_out = args.outpath
    case_name = args.casename
    path_ref = os.path.join(path_in, 'Stats.'+case_name+'.nc')

    # (0) Namelist File
    # path = '../test_ZGILS6/'
    # path_ref = os.path.join(path, 'Stats.ZGILS_S6_1xCO2_SST_FixSub.nc')
    # case_name = 'ZGILS_S6_1xCO2_SST_FixSub'

    # path = '../test/'
    # path_ref = os.path.join(path, 'Stats.Bomex.nc')
    # do_everything(path, path_ref)

    # path = '../test_bomex_n1024/'
    # path_ref = os.path.join(path, 'Stats.Bomex.nc')
    # case_name = 'Bomex'

    # path = '../test_bomex/'
    # case_name = 'Bomex'
    # path_ref = os.path.join(path, 'Stats.'+case_name+'.nc')


    nml = simplejson.loads(open(os.path.join(path_in, case_name+'.in')).read())
    nz = nml['grid']['nz']
    dz = nml['grid']['dz']

    ClCl = CloudClosure.CloudClosure()
    ClCl.initialize(path_in, path_ref, case_name)
    # ClCl.verification_CC(path, path_ref)

    files = os.listdir(os.path.join(path_in, 'fields'))
    # ncomp_range = [1, 2, 3, 4, 5, 6, 8, 10]
    ncomp_range = [1, 2]

    # ZGILS 6
    files = ['1382400.nc']
    # files_ = [1317600, 1339200, 1360800, 1382400]  # ZGILS 6
    krange = np.asarray([16, 37, 50])
    # ZGILS S12
    # files = ['432000.nc']
    # krange = np.asarray([35,40,45])
    # DYCOMS RF01
    # krange = np.asarray([140,166])
    # files = ['13800.nc', '14400.nc']
    # files = ['13800.nc']
    # DYCOMS RF02
    # krange = np.asarray([120, 170])
    # files = ['18000.nc']
    # Bomex large, kyle
    # krange = np.asarray([27, 91])
    # Bomex 170413
    # krange = np.asarray([20, 39])
    krange = np.asarray([20])
    files = ['21600.nc']
    # Bomex test
    # files = ['14400.nc']
    # krange = np.asarray([10, 17, 20, 25, 50])
    # krange = np.asarray([10, 20, 50])
    N = len(files)
    print('Found the following directories', files, N)
    print('zrange: ' + str(krange * dz) + ', dz: ' + str(dz), 'ncomp: ' + str(ncomp_range))
    print('')
    ClCl.predict_pdf(files, path_in, path_out, path_ref, ncomp_range, krange, nml)

    # for ncomp in [1,2,3,4,5,6,7,8,9,10]:
    # for ncomp in [1, 3, 5, 10, 15]:
    # # for ncomp in [15]:
    # # ncomp = 1
    #     # Bomex
    #     # krange = np.arange(6, 20, 8)
    #     krange = np.asarray([25, 85])
    #     # krange = np.asarray([10, 17, 25, 50])
    #     # krange = np.asarray([20])
    #     print('zrange: '+str(krange*dz)+', dz: '+str(dz), 'ncomp: '+str(ncomp))
    #     print('')
    #     ClCl.predict_pdf(files, path, path_ref, ncomp, krange, nml)





    return



if __name__ == "__main__":
    main()