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
    print('Found the following files: ' + str(files))
    # ncomp_range = [1, 2, 3, 4, 5, 6, 7, 8]
    ncomp_range = [1, 2, 3, 4, 5]

    # ZGILS 6
    # files = ['1382400.nc']
    # files_ = [1317600, 1339200, 1360800, 1382400]  # ZGILS 6
    # files = files[0:25:2]
    # krange = np.asarray([25,25,40,50,60,65])
    # ZGILS S12
    # files = ['432000.nc']
    # files = ['345600.nc', '432000.nc', '518400.nc', '604800.nc', '691200.nc']
    # krange = np.asarray([35,40,45])
    # DYCOMS RF01
    # krange = np.asarray([140,150,160,166,180])
    # files = ['10800.nc', '12600.nc', '14400.nc']
    # files = ['14400.nc']
    # files = ['3600.nc']
    # DYCOMS RF02
    # krange = np.asarray([120, 170])
    krange = np.asarray([120, 140, 160, 170, 200])
    # files = ['18000.nc', '19800.nc', '21600.nc']
    files = ['10800.nc']
    # Bomex large, kyle
    # krange = np.asarray([27, 91])
    # Bomex 170314_weno7
    # krange = np.asarray([15, 20, 25, 30, 35, 40, 45])
    # files = ['18000.nc', '19800.nc', '21600.nc']
    # Bomex test
    # files = ['21600.nc']
    # krange = np.asarray([10, 17, 20, 25, 50])
    # krange = np.asarray([20, 25])
    # krange = np.asarray([18,30,38])
    # TRMM
    # files = ['1012600.nc', '1014400.nc', '1016200.nc']
    # krange = np.asarray([10, 15, 20, 30, 40, 50, 60])


    N = len(files)
    n_sample = np.int(1e6)
    print('Use the following files: ', files, N)
    print('zrange: ' + str(krange * dz) + ', dz: ' + str(dz), 'ncomp: ' + str(ncomp_range))
    print('')
    ClCl.predict_pdf(files, path_in, path_out, path_ref, n_sample, ncomp_range, krange, nml)

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