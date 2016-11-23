__author__ = 'meyerbe'
from netCDF4 import Dataset


def test():
    rootgrp = Dataset('test.nc','w',format='NETCDF4')

    # create Groups, Sub-Groups
    grp_a = rootgrp.createGroup('test_data')
    print(rootgrp.groups)

    grp_a_a = grp_a.createGroup('test_data_group')
    print(rootgrp.groups)
    print('___')
    print(grp_a.groups)

    print(grp_a.groups.values)

    # create Variables
    # (1) define dimensions
    time = grp_a_a.createDimension('time', size =None)
    print('Dimensions:', grp_a_a.dimensions)

    grp_a_a.createDimension('lat', 73)
    print('Dimensions:', grp_a_a.dimensions)

    print('Time dimension:', 'time'.__len__())
    print('Time', len('time'))
    print('Time', len(time))
    print('Lat dimension:', 'lat'.__len__())
    #print(values)
    #print('Time unl', is_unlimited(grp_a_a.dimensions['time']))
    dim = len(grp_a_a.dimensions['time'])
    print('dd', dim)

    print(grp_a_a.dimensions.values())


    # (2) create Variable
    times = grp_a_a.createVariable('time','f8',('time',))
    latitudes = grp_a_a.createVariable('latitude','f4',('lat',))





    rootgrp.close()

    # ___________
    print('___________')

    rootgrp = Dataset('test.nc','w',format='NETCDF4')

    # create Groups, Sub-Groups
    grp_a = rootgrp.createGroup('test_data')
    print(rootgrp.groups)
    grp_a_a = grp_a.createGroup('test_data_group')
    # create Variables
    # (1) define dimensions
    time = grp_a_a.createDimension('time', size =None)
    print('Dimensions:', grp_a_a.dimensions)
    print('Time', len(time))


def walktree(top):
    values = top.groups.values()
    yield values
    for value in top.groups.values():
        for children in walktree(value):
            yield children


f = Dataset('100/0.nc')
for grp in walktree(f):
    for a in grp:
        print(a)
#print(f.groups)

