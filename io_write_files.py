import netCDF4 as nc
import numpy as np
import json as  simplejson
import os

# ____________________
def dump_pickle(data,out_path,file_name):
    data_ = (1.4,42)
    # output = open(os.path.join(out_path,'data.pkl'), 'w')
    output = open(os.path.join(out_path, file_name), 'w')
    pickle.dump(data, output)
    output.close()
    return
def test_pickle(in_path,file_name):
    print('')
    print('------- test pickle ------')
    fullpath_in = os.path.join(in_path,file_name)
    f = open(fullpath_in)
    data = pickle.load(f)
    print(data)
    print('')
    var = data['w']
    print(var)
    print
    ''
    means_ = var['means']
    print(means_)
    print('-------------------------')
    print('')
    retur