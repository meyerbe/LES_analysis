import numpy as np
import matplotlib.cm as cm
import pylab as plt

from netCDF4 import Dataset
import netCDF4 as nc

from h5py import File

import time
import sys

#include "Parameters.pxi"

g = 9.80665   #[m/s^2]





def read_in_netcdf(variable_name, fullpath_in):
    print(fullpath_in, variable_name)
    ncfile = Dataset(fullpath_in,'r')
    data = ncfile.variables[variable_name][:]    
    ncfile.close()
    return data

def read_in_hdf5(variable_name, group_name, fullpath_in):
    f = File(fullpath_in)
    print('read in:', group_name, variable_name, fullpath_in)
    profiles_group = f[group_name]
    variable_dataset = profiles_group[variable_name]
    variable_dataset_shape = variable_dataset.shape
    variable = np.ndarray(shape = variable_dataset_shape)
    for t in range(variable_dataset_shape[0]):
        if group_name == "profiles":
            variable[t,:] = variable_dataset[t, :]
        elif group_name == "correlations":
            variable[t,:] = variable_dataset[t, :]
        elif group_name == "timeseries":
            variable[t] = variable_dataset[t]
        elif group_name == "fields":
            variable[0,:,:,:] = variable_dataset[T[t],:,:,:]
    f.close()
    return variable

def covariance(data1,data2,profile1,profile2,tt):
    nx = data1.shape[0]
    ny = data1.shape[1]
    nz = data1.shape[2]
    if nz != data2.shape[2] or nx != data2.shape[0]:
        print('data size error')
        sys.exit()
    print('nx', nx, 'nz', nz)

    cov = np.ndarray(shape = (nz,))
    mean1 = profile1[tt,:]
    mean2 = profile2[tt,:]
    #mean1 = np.zeros(shape = (nz,))
    #mean2 = np.zeros(shape = (nz,))
    prod = np.zeros(shape = (nz,))
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                prod[k] += data1[i,j,k]*data2[i,j,k]
                #mean1 += data1[i,j,k]
                #mean2 += data2[i,j,k]
    prod /= nx*ny
    #mean1 /= nx*ny
    #mean2 /= nx*ny
    for k in range(nz):
        cov[k] = prod[k] - mean1[k]*mean2[k]
    return cov

def variance(data,profile,tt):
    nx = data.shape[0]
    ny = data.shape[1]
    nz = data.shape[2]

    var = np.ndarray(shape = (nz,))
    mean = profile[tt,:]
    prod = np.zeros(shape = (nz,))
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                prod[k] += data[i,j,k]*data[i,j,k]
    prod /= nx*ny
    for k in range(nz):
        var[k] = prod[k] - mean[k]*mean[k]
    return var

def compare(data,profile):
    nx = data.shape[0]
    ny = data.shape[1]
    nz = data.shape[2]

    mean = np.zeros(shape = (nz,))
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                mean[k] += data[i,j,k]
    mean /= nx*ny
    print('profile',profile.shape,'mean',mean.shape)
    print('mean max',np.amax(np.abs(mean)))
    plt.figure(1,figsize=(10,10))
    plt.subplot(1,2,1)
    plt.plot(profile[tt,:],'-x',label='profile')
    plt.plot(mean,label='mean')
    plt.legend()
    plt.subplot(1,2,2)
    plt.plot(mean,label='mean')
    plt.legend()
    plt.savefig('comp.png')
    return

def gradient(data,dz):
    nz = data.shape[-1]
    ddz = 1/(2*dz)
    data_grad = np.zeros(shape = (nz,))
    for k in range(1,nz-1):
        data_grad[k] = ddz*(data[k+1]-data[k-1])

    return data_grad


def gradient3d(data,dz):
    nz = data.shape[-1]
    ddz = 1/(2*dz)
    data_grad = np.zeros(shape = data.shape)
    for k in range(1,nz-1):
        data_grad[:,:,k] = ddz*(data[:,:,k+1]-data[:,:,k-1])
    return data_grad


def stats(data1,data2,profile1,profile2,tt):
    # data1 = w
    # data2 \in {u,v}
    nx = data1.shape[0]
    ny = data1.shape[1]
    nz = data1.shape[2]
    if nz != data2.shape[2] or nx != data2.shape[0]:
        print('data size error')
        sys.exit()
    
    mean1 = profile1[tt,:]
    mean2 = profile2[tt,:]

    cov = np.ndarray(shape = (nz,))
    skew = np.ndarray(shape = (nz,))
    
    prod = np.zeros(shape = (nz,))
    var = np.zeros(shape = (nz,))
    triple = np.zeros(shape = (nz,))
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                prod[k] += data1[i,j,k]*data2[i,j,k]
                var[k] += data2[i,j,k]*data2[i,j,k]
                triple[k] += data1[i,j,k]*data2[i,j,k]*data2[i,j,k]
    prod /= nx*ny
    var /= nx*ny
    triple /= nx*ny
    #for k in range(nz):
    #cov[k] = prod[k] - mean1[k]*mean2[k]
    #skew[k] = triple[k] - 2*prod[k]*mean2[k] - mean1[k]*var[k] + 2*mean1[k]*mean2[k]*mean2[k]
    cov = prod - mean1*mean2
    skew = triple - 2*prod*mean2 - mean1*var + 2*mean1*mean2*mean2
    
    return var, cov, skew


def average(data,tt):
    nx = data.shape[0]
    ny = data.shape[1]
    nz = data.shape[2]
    
    mean = np.zeros(nz)
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                mean[k] += data[i,j,k]
    mean /= nx*ny
    return mean


# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------

T = '7200'
field = '7200.nc'
tt = 2
dt_profile = 3600
dt_field = 3600
data_full = 'dcbl_full.statistics.hdf5'
data_ql = 'dcbl_ql.statistics.hdf5'

path_full = 'dcbl/150918_2nd/full_SGSfull/'
path_ql = 'dcbl/150918_2nd/ql_SGSfull/'
path_out = 'dcbl/150918_2nd/'
print(path_full)
print(path_ql)

time_full = read_in_hdf5('time', "timeseries", path_full + data_full)
time_ql = read_in_hdf5('time', "timeseries", path_ql + data_ql)
print('time full',time_full[1:5]/3600, 't_max: ', time_full[-1]/3600, np.shape(time_full))
print('time ql',time_ql[1:5]/3600, 't_max: ', time_ql[-1]/3600, np.shape(time_ql))

## READ IN FIELD DATA
tic = time.time()
th_l_full = read_in_netcdf('theta_l', path_full + 'fields/' + field)
th_l_ql = read_in_netcdf('theta_l', path_ql + 'fields/' + field)
#u_full = read_in_netcdf('u', path_full + 'fields/' + field)
#u_ql = read_in_netcdf('u', path_ql + 'fields/' + field)
#v_full = read_in_netcdf('v', path_full + 'fields/' + field)
#v_ql = read_in_netcdf('v', path_ql + 'fields/' + field)
w_full = read_in_netcdf('w', path_full + 'fields/' + field)
w_ql = read_in_netcdf('w', path_ql + 'fields/' + field)
toc = time.time()
print(toc-tic)

## READ IN PROFILES
th_mean_full = read_in_hdf5('potential_temperature', "profiles", path_full + data_full)
th_mean_ql = read_in_hdf5('potential_temperature', "profiles", path_ql + data_ql)
#u_mean_full = read_in_hdf5('u', "profiles", path_full + data_full)
#u_mean_ql = read_in_hdf5('u', "profiles", path_ql + data_ql)
#v_mean_full = read_in_hdf5('v', "profiles", path_full + data_full)
#v_mean_ql = read_in_hdf5('v', "profiles", path_ql + data_ql)
w_mean_full = read_in_hdf5('w', "profiles", path_full + data_full)
w_mean_ql = read_in_hdf5('w', "profiles", path_ql + data_ql)

nx = w_full.shape[0]
ny = w_full.shape[1]
nz = w_full.shape[2]
##
nz = nz/2

dz = 25

## TEST
#compare(u_full,u_mean_full)
#compare(th_l_full,th_mean_full)


# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------

var_th_full, wth_full, wthth_full = stats(w_full,th_l_full,w_mean_full,th_mean_full,tt)
var_th_ql, wth_ql, wthth_ql = stats(w_ql,th_l_ql,w_mean_ql,th_mean_ql,tt)

a = read_in_hdf5('potential_temperature_resolved_variance',"profiles",path_full + data_full)
b = read_in_hdf5('potential_temperature_resolved_variance',"profiles",path_ql + data_ql)
var_th_full[:] = a[tt,:]
var_th_ql[:] = b[tt,:]

## TH_VARIANCE & TH_VARIANCE RATE
R_full = np.ndarray(shape = var_th_full.shape)
R_ql = np.ndarray(shape = var_th_ql.shape)
#th_var_full = read_in_hdf5('potential_temperature_resolved_variance',"profiles",path_full + data_full)
ddt = 1/(dt_profile)
for t in range(1,var_th_full.shape[0]):
    R_full[:] = ddt*(var_th_full[t-1]-var_th_full[t])
    R_ql[:] = ddt*(var_th_ql[t-1]-var_th_ql[t])



## BUOYANCY PRODUCTION
th_grad_full = gradient(th_mean_full[tt,:],dz)
th_grad_ql = gradient(th_mean_ql[tt,:],dz)

B_full = np.zeros(shape = th_grad_full.shape)
B_ql = np.zeros(shape = th_grad_ql.shape)
B_full[:] = - 2*th_grad_full[:]*wth_full[:]
B_ql[:] = -2*th_grad_ql[:]*wth_ql[:]



## TURBULENT TRANSPORT
T_full = np.zeros((dz,))
T_ql = np.zeros((dz,))
T_full = - gradient(wthth_full,dz)
T_ql = - gradient(wthth_ql,dz)


#sys.exit()

## DISSIPATION
#nu_full = read_in_netcdf('eddy_viscosity',path_full + 'fields/' + field)
#nu_ql = read_in_netcdf('eddy_viscosity',path_ql + 'fields/' + field)
nu_profile_full = read_in_hdf5('eddy_diffusivity',"profiles",path_full + data_full)
nu_profile_ql = read_in_hdf5('eddy_diffusivity',"profiles",path_ql + data_ql)

grad_th_full = gradient3d(th_l_full,dz)
grad_th_ql = gradient3d(th_l_ql,dz)

grad_th_full2 = np.multiply(grad_th_full,grad_th_full)
grad_th_ql2 = np.multiply(grad_th_ql,grad_th_ql)

diss_th_full = average(grad_th_full2,tt)
diss_th_ql = average(grad_th_ql2,tt)

epsilon_full = np.zeros(nz)
epsilon_full = - 2 * nu_profile_full[tt,:] * diss_th_full[:]
epsilon_ql = np.zeros(nz)
epsilon_ql = - 2 * nu_profile_ql[tt,:] * diss_th_ql[:]


#zvector = dz*np.linspace(0,nz-1,nz)
#plt.figure(1,figsize=(25,10))
#plt.subplot(1,3,1)
#plt.plot(epsilon_full,zvector,linewidth=3,label=r'$\epsilon$ full')
#plt.plot(nu_profile_full[tt,:]*diss_u_full,zvector,linewidth=2,label=r'$\epsilon$ u')
#plt.plot(nu_profile_full[tt,:]*diss_v_full,zvector,linewidth=2,label=r'$\epsilon$ v')
#plt.plot(nu_profile_full[tt,:]*diss_w_full,zvector,linewidth=2,label=r'$\epsilon$ w')
#plt.title('dissipation: full LES')
#plt.legend()
#plt.subplot(1,3,2)
#plt.plot(epsilon_ql,zvector,linewidth=3,label=r'$\epsilon$ ql')
#plt.plot(nu_profile_ql[tt,:]*diss_u_ql,zvector,linewidth=2,label=r'$\epsilon$ u')
#plt.plot(nu_profile_ql[tt,:]*diss_v_ql,zvector,linewidth=2,label=r'$\epsilon$ v')
#plt.plot(nu_profile_ql[tt,:]*diss_w_ql,zvector,linewidth=2,label=r'$\epsilon$ w')
#plt.legend()
#plt.title('dissipation: QL LES')
#plt.subplot(1,3,3)
#plt.plot(nu_profile_ql[tt,:],zvector,'x-',linewidth=2,label='nu ql')
#plt.plot(nu_profile_full[tt,:],zvector,'-',linewidth=2,label='nu full')
#plt.legend()
#plt.title('viscosities')
##plt.plot(diss_u_full,zvector,linewidth=2,label='diss u full')
##plt.plot(diss_v_full,zvector,linewidth=2,label='diss v full')
##plt.plot(diss_w_full,zvector,linewidth=2,label='diss w full')
##plt.plot(diss_u_ql,zvector,'--',linewidth=2,label='diss u ql')
##plt.plot(diss_v_ql,zvector,'--',linewidth=2,label='diss v ql')
##plt.plot(diss_w_ql,zvector,'--',linewidth=2,label='diss w ql')
##plt.legend()
#
#plt.savefig('TKE_dissipation.png')
#plt.close()



# ------------------------------------------------------------------------------------------

zvector = dz*np.linspace(0,nz-1,nz)
plt.figure(1,figsize=(20,10))
plt.subplot(1,2,1)
plt.plot(R_full[0:nz],zvector,linewidth=2,label='R full')
plt.plot(B_full[0:nz],zvector,linewidth=2,label='B full')
plt.plot(T_full[0:nz],zvector,linewidth=2,label='T full')
plt.plot(epsilon_full[0:nz],zvector,linewidth=2,label='D full')
plt.title(r'$\overline{(\theta-\overline{\theta})^2}$: full LES (' + path_full + ', ' + field + ')')
plt.legend()
plt.subplot(1,2,2)
plt.plot(R_ql[0:nz],zvector,linewidth=2,label='R ql')
plt.plot(B_ql[0:nz],zvector,linewidth=2,label='B ql')
plt.plot(T_ql[0:nz],zvector,linewidth=2,label='T ql')
plt.plot(epsilon_ql[0:nz],zvector,linewidth=2,label='D ql')
plt.title(r'$\overline{(\theta-\overline{\theta})^2}$: QL LES (' + path_ql + ', ' + field + ')')
plt.legend()
plt.savefig('Th_Budget_' + T + '.png')
plt.savefig(path_out + 'Th_Budget_' + T + '.png')
plt.close()


sys.exit()


