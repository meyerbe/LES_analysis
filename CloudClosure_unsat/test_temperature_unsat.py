import sys, os
import json as simplejson
import netCDF4 as nc

import numpy as np

from CC_thermodynamic_functions import do_everything
from CC_thermodynamic_functions import do_everything_with_pycles

def main():
    # path = '../test_ZGILS6/'
    # path_ref = os.path.join(path, 'Stats.ZGILS_S6_1xCO2_SST_FixSub.nc')
    # do_everything(path, path_ref)
    # path = '../test/'
    # path_ref = os.path.join(path, 'Stats.Bomex.nc')
    # do_everything(path, path_ref)

    path = '../test_bomex/'
    path_ref = os.path.join(path, 'Stats.Bomex.nc')
    # do_everything(path, path_ref)
    do_everything_with_pycles(path, path_ref)

    return





def read_in_netcdf(var_name, group_name, path):
    print('read in '+ var_name + ': ' +path)
    rootgrp = nc.Dataset(path, 'r')
    grp = rootgrp.groups[group_name]
    var = grp.variables[var_name]
    return

if __name__ == "__main__":
    main()