''' offline calculation of TKE from LES output (old LES; hdf5) '''

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




def stats_new(data1,data2,tt):
    nx = data1.shape[0]
    ny = data1.shape[1]
    nz = data1.shape[2]
    if nz != data2.shape[2] or nx != data2.shape[0]:
        print('data size error')
        sys.exit()

    


def stats(data1,data2,profile1,profile2,tt):
    # data1 = w
    # data2 \in {u,v}
    nx = data1.shape[0]
    ny = data1.shape[1]
    nz = data1.shape[2]
    if nz != data2.shape[2] or nx != data2.shape[0]:
        print('data size error')
        sys.exit()
    
    #mean1 = profile1[tt,:]
    #mean2 = profile2[tt,:]

    cov = np.ndarray(shape = (nz,))
    var = np.ndarray(shape = (nz,))
    skew = np.ndarray(shape = (nz,))

    mean1 = np.zeros(shape = (nz,))
    mean2 = np.zeros(shape = (nz,))
    co_prod = np.zeros(shape = (nz,))
    prod = np.zeros(shape = (nz,))
    triple = np.zeros(shape = (nz,))
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                mean1[k] += data1[i,j,k]
                mean2[k] += data2[i,j,k]
                co_prod[k] += data1[i,j,k]*data2[i,j,k]
                prod[k] += data2[i,j,k]*data2[i,j,k]
                triple[k] += data1[i,j,k]*data2[i,j,k]*data2[i,j,k]
    mean1 /= nx*ny
    mean2 /= nx*ny
    co_prod /= nx*ny
    prod /= nx*ny
    triple /= nx*ny
    #for k in range(nz):
    #cov[k] = prod[k] - mean1[k]*mean2[k]
    #skew[k] = triple[k] - 2*prod[k]*mean2[k] - mean1[k]*var[k] + 2*mean1[k]*mean2[k]*mean2[k]
    cov = co_prod - mean1*mean2
    var = prod - mean1*mean1
    skew = triple - 2*co_prod*mean2 - mean1*prod + 2*mean1*mean2*mean2
    
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
print('time full',time_full[0:5]/3600, 't_max: ', time_full[-1]/3600, np.shape(time_full))
print('time ql',time_ql[0:5]/3600, 't_max: ', time_ql[-1]/3600, np.shape(time_ql))

## READ IN FIELD DATA
tic = time.time()
th_l_full = read_in_netcdf('theta_l', path_full + 'fields/' + field)
th_l_ql = read_in_netcdf('theta_l', path_ql + 'fields/' + field)
u_full = read_in_netcdf('u', path_full + 'fields/' + field)
u_ql = read_in_netcdf('u', path_ql + 'fields/' + field)
v_full = read_in_netcdf('v', path_full + 'fields/' + field)
v_ql = read_in_netcdf('v', path_ql + 'fields/' + field)
w_full = read_in_netcdf('w', path_full + 'fields/' + field)
w_ql = read_in_netcdf('w', path_ql + 'fields/' + field)
#p_full = read_in_netcdf('p',path_full + 'fields/' + field)
#p_ql = read_in_netcdf('p',path_ql + 'fields/' + field)
toc = time.time()
print(toc-tic)

## READ IN PROFILES
th_mean_full = read_in_hdf5('potential_temperature', "profiles", path_full + data_full)
th_mean_ql = read_in_hdf5('potential_temperature', "profiles", path_ql + data_ql)
u_mean_full = read_in_hdf5('u', "profiles", path_full + data_full)
u_mean_ql = read_in_hdf5('u', "profiles", path_ql + data_ql)
v_mean_full = read_in_hdf5('v', "profiles", path_full + data_full)
v_mean_ql = read_in_hdf5('v', "profiles", path_ql + data_ql)
w_mean_full = read_in_hdf5('w', "profiles", path_full + data_full)
w_mean_ql = read_in_hdf5('w', "profiles", path_ql + data_ql)
p_mean_full = read_in_hdf5('pressure',"profiles", path_full + data_full)
p_mean_ql = read_in_hdf5('pressure',"profiles", path_ql + data_ql)

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
## TKE TRANSPORT
print('TKE Transport')
var_u_full, wu_full, wuu_full = stats(w_full,u_full,w_mean_full,u_mean_full,tt)
var_v_full, wv_full, wvv_full = stats(w_full,v_full,w_mean_full,v_mean_full,tt)
var_w_full, ww_full, www_full = stats(w_full,w_full,w_mean_full,w_mean_full,tt)
var_u_ql, wu_ql, wuu_ql = stats(w_ql,u_ql,w_mean_ql,u_mean_ql,tt)
var_v_ql, wv_ql, wvv_ql = stats(w_ql,v_ql,w_mean_ql,v_mean_ql,tt)
var_w_ql, ww_ql, www_ql = stats(w_ql,w_ql,w_mean_ql,w_mean_ql,tt)

trans_u_full = gradient(wuu_full,dz)
trans_v_full = gradient(wvv_full,dz)
trans_w_full = gradient(www_full,dz)
trans_u_ql = gradient(wuu_ql,dz)
trans_v_ql = gradient(wvv_ql,dz)
trans_w_ql = gradient(www_ql,dz)

T_full = np.zeros(shape = trans_u_full.shape)
T_ql = np.zeros(shape = trans_u_ql.shape)
T_full = - (trans_u_full + trans_v_full + trans_w_full)
T_ql = - (trans_u_ql + trans_v_ql + trans_w_ql)



## TKE & TKE RATE
print('TKE Rate')
#u2_full = variance(u_full,u_mean_full,tt)
#v2_full = variance(v_full,v_mean_full,tt)
#w2_full = variance(w_full,w_mean_full,tt)
u2_full = var_u_full
v2_full = var_v_full
w2_full = var_w_full
u_var_full = read_in_hdf5('u_resolved_variance',"profiles",path_full + data_full)
v_var_full = read_in_hdf5('v_resolved_variance',"profiles",path_full + data_full)
w_var_full = read_in_hdf5('w_resolved_variance',"profiles",path_full + data_full)
e_full = np.ndarray(shape = (nz,))
e_full_ = np.ndarray(shape = u_var_full.shape)
R_full = np.ndarray(shape = u_var_full.shape)

#u2_ql = variance(u_ql,u_mean_ql,tt)
#v2_ql = variance(v_ql,v_mean_ql,tt)
#w2_ql = variance(w_ql,w_mean_ql,tt)
u2_ql = var_u_ql
v2_ql = var_v_ql
w2_ql = var_w_ql
u_var_ql = read_in_hdf5('u_resolved_variance',"profiles",path_ql + data_ql)
v_var_ql = read_in_hdf5('v_resolved_variance',"profiles",path_ql + data_ql)
w_var_ql = read_in_hdf5('w_resolved_variance',"profiles",path_ql + data_ql)
e_ql = np.ndarray(shape = (nz,))
e_ql_ = np.ndarray(shape = u_var_ql.shape)
R_ql = np.ndarray(shape = u_var_ql.shape)

e_full = 0.5*(u2_full + v2_full + w2_full)
e_ql = 0.5*(u2_ql + v2_ql + w2_ql)

#ddt = 1/(2*dt_profile)
ddt = 1/dt_profile
for t in range(0,u_var_full.shape[0]):
    e_full_[t,:] = 0.5*(u_var_full[t,:] + v_var_full[t,:] + w_var_full[t,:])
    e_ql_[t,:] = 0.5*(u_var_ql[t,:] + v_var_ql[t,:] + w_var_ql[t,:])
    if t > 0:
        R_full[t,:] = ddt*(e_full_[t-1,:] - e_full_[t,:])
        #R_ql[t,:] = ddt*(e_ql_[t-1,:]-e_ql_[t+1,:])
        R_ql[t,:] = ddt*(e_ql_[t-1,:]-e_ql_[t,:])
#e_full_ = 0.5*(u_var_full[tt,:] + v_var_full[tt,:] + w_var_full[tt,:])

plt.figure()
plt.subplot(1,3,1)
plt.plot(u2_full,'x-',label='u2')
plt.plot(u_var_full[tt,0:nz],label='u res_var')
plt.plot(v2_full,'x-',label='v2')
plt.plot(v_var_full[tt,0:nz],label='v res_var')
plt.plot(w2_full,'x-',label='w2')
plt.plot(w_var_full[tt,0:nz],label='w res_var')
plt.legend()
plt.subplot(1,3,2)
plt.plot(e_full[0:nz],label = 'e, full')
plt.plot(e_full_[tt,0:nz],label = 'e res var, full')
plt.plot(e_ql[0:nz],label = 'e, ql')
plt.plot(e_ql_[tt,0:nz],label = 'e res var, ql')
plt.legend()
plt.subplot(1,3,3)
plt.plot(R_full[tt,0:nz],label = 'rate tke, full')
plt.plot(R_ql[tt,0:nz],label = 'rate tke, ql')
plt.legend()
plt.savefig('TKE_comp_var_' + T + '.png')
plt.savefig(path_out + 'TKE_comp_var_' + T + '.png')
plt.close()




## BUOYANCY PRODUCTION
print('Buoyancy Production')
heat_flux_full = covariance(th_l_full,w_full,th_mean_full,w_mean_full,tt)
heat_flux_ql = covariance(th_l_ql,w_ql,th_mean_ql,w_mean_ql,tt)

B_full = np.zeros(shape = heat_flux_full.shape)
B_ql = np.zeros(shape = heat_flux_ql.shape)
B_full[:] = g*1/th_mean_full[tt,:]*heat_flux_full[:]
B_ql[:] = g*1/th_mean_ql[1,:]*heat_flux_ql[:]




## SHEAR PRODUCTION
print('Shear Production')
#wu_full = covariance(u_full,w_full,u_mean_full,w_mean_full,tt)
#wv_full = covariance(v_full,w_full,v_mean_full,w_mean_full,tt)
#wu_ql = covariance(u_ql,w_ql,u_mean_ql,w_mean_ql,tt)
#wv_ql = covariance(v_ql,w_ql,v_grad_ql,w_mean_ql,tt)

u_grad_full = gradient(u_mean_full[tt,:],dz)
u_grad_ql = gradient(u_mean_ql[tt,:],dz)
v_grad_full = gradient(v_mean_full[tt,:],dz)
v_grad_ql = gradient(v_mean_ql[tt,:],dz)

S_full = np.zeros(shape = wu_full.shape)
S_ql = np.zeros(shape = wu_ql.shape)
S_full[:] = -(wu_full[:]*u_grad_full[:] + wv_full*v_grad_full[:])
S_ql[:] = -(wu_ql[:]*u_grad_ql[:] + wv_ql*v_grad_ql[:])



#sys.exit()

## DISSIPATION
print('Dissipation')
#nu_full = read_in_netcdf('eddy_viscosity',path_full + 'fields/' + field)
#nu_ql = read_in_netcdf('eddy_viscosity',path_ql + 'fields/' + field)
nu_profile_full = read_in_hdf5('eddy_viscosity',"profiles",path_full + data_full)
nu_profile_ql = read_in_hdf5('eddy_viscosity',"profiles",path_ql + data_ql)

grad_u_full = gradient3d(u_full,dz)
grad_v_full = gradient3d(v_full,dz)
grad_w_full = gradient3d(w_full,dz)

grad_u_full2 = np.multiply(grad_u_full,grad_u_full)
grad_v_full2 = np.multiply(grad_v_full,grad_v_full)
grad_w_full2 = np.multiply(grad_w_full,grad_w_full)

diss_u_full = average(grad_u_full2,tt)
diss_v_full = average(grad_v_full2,tt)
diss_w_full = average(grad_w_full2,tt)


grad_u_ql = gradient3d(u_ql,dz)
grad_v_ql = gradient3d(v_ql,dz)
grad_w_ql = gradient3d(w_ql,dz)

grad_u_ql2 = np.multiply(grad_u_ql,grad_u_ql)
grad_v_ql2 = np.multiply(grad_v_ql,grad_v_ql)
grad_w_ql2 = np.multiply(grad_w_ql,grad_w_ql)

diss_u_ql = average(grad_u_ql2,tt)
diss_v_ql = average(grad_v_ql2,tt)
diss_w_ql = average(grad_w_ql2,tt)

epsilon_full = np.zeros(nz)
epsilon_full = - nu_profile_full[tt,:] * (diss_u_full + diss_v_full + diss_w_full)
epsilon_ql = np.zeros(nz)
epsilon_ql = - nu_profile_ql[tt,:] * (diss_u_ql + diss_v_ql + diss_w_ql)

zvector = dz*np.linspace(0,nz-1,nz)
plt.figure(1,figsize=(25,10))
plt.subplot(1,3,1)
plt.plot(epsilon_full[0:nz],zvector,linewidth=3,label=r'$\epsilon$ full')
plt.plot(nu_profile_full[tt,0:nz]*diss_u_full[0:nz],zvector,linewidth=2,label=r'$\epsilon$ u')
plt.plot(nu_profile_full[tt,0:nz]*diss_v_full[0:nz],zvector,linewidth=2,label=r'$\epsilon$ v')
plt.plot(nu_profile_full[tt,0:nz]*diss_w_full[0:nz],zvector,linewidth=2,label=r'$\epsilon$ w')
plt.title('dissipation: full LES')
plt.legend()
plt.subplot(1,3,2)
plt.plot(epsilon_ql[0:nz],zvector,linewidth=3,label=r'$\epsilon$ ql')
plt.plot(nu_profile_ql[tt,0:nz]*diss_u_ql[0:nz],zvector,linewidth=2,label=r'$\epsilon$ u')
plt.plot(nu_profile_ql[tt,0:nz]*diss_v_ql[0:nz],zvector,linewidth=2,label=r'$\epsilon$ v')
plt.plot(nu_profile_ql[tt,0:nz]*diss_w_ql[0:nz],zvector,linewidth=2,label=r'$\epsilon$ w')
plt.legend()
plt.title('dissipation: QL LES')
plt.subplot(1,3,3)
plt.plot(nu_profile_ql[tt,0:nz],zvector,'-',linewidth=2,label='nu ql')
plt.plot(nu_profile_full[tt,0:nz],zvector,'-',linewidth=2,label='nu full')
plt.legend()
plt.title('viscosities')
#plt.plot(diss_u_full,zvector,linewidth=2,label='diss u full')
#plt.plot(diss_v_full,zvector,linewidth=2,label='diss v full')
#plt.plot(diss_w_full,zvector,linewidth=2,label='diss w full')
#plt.plot(diss_u_ql,zvector,'--',linewidth=2,label='diss u ql')
#plt.plot(diss_v_ql,zvector,'--',linewidth=2,label='diss v ql')
#plt.plot(diss_w_ql,zvector,'--',linewidth=2,label='diss w ql')
#plt.legend()

plt.savefig('TKE_dissipation_' + T + '.png')
plt.savefig(path_out + 'TKE_dissipation_' + T + '.png')
plt.close()



# ------------------------------------------------------------------------------------------
print('!!!',nz,dz)
#zvector = dz*np.linspace(0,nz-1,nz)
plt.figure(1,figsize=(20,10))
plt.subplot(1,2,1)
plt.plot(S_full[0:nz],zvector,linewidth=2,label='S full')
plt.plot(B_full[0:nz],zvector,linewidth=2,label='B full')
#plt.plot(P_full[0:nz],zvector,linewidth=2,label='P full')
plt.plot(T_full[0:nz],zvector,linewidth=2,label='T full')
plt.plot(R_full[tt,0:nz],zvector,linewidth=2,label='R full')
plt.plot(epsilon_full[0:nz],zvector,linewidth=2,label='D full')
plt.title('TKE: full LES (' + path_full + ', ' + field + ')')
plt.legend()
plt.subplot(1,2,2)
plt.plot(S_ql[0:nz],zvector,linewidth=2,label='S ql')
plt.plot(B_ql[0:nz],zvector,linewidth=2,label='B ql')
#plt.plot(P_ql[0:nz],zvector,linewidth=2,label='P full')
plt.plot(T_ql[0:nz],zvector,linewidth=2,label='T ql')
plt.plot(R_ql[tt,0:nz],zvector,linewidth=2,label='R ql')
plt.plot(epsilon_ql[0:nz],zvector,linewidth=2,label='D ql')
plt.title('TKE: QL LES (' + path_ql + ', ' + field + ')')
plt.legend()
plt.savefig('TKE_Budget_' + T + '.png')
plt.savefig(path_out + 'TKE_Budget_' + T + '.png')
plt.close()


sys.exit()
## PRESSURE CORRELATION
print('Pressure Correlation')
residual_full = np.ndarray(nz)
residual_full[:] = - (R_full[tt,:] - B_full[:] + S_full[:] + T_full[:])

density_full = read_in_hdf5('density',"profiles",path_full + data_full)
#density_ql = read_in_hdf5('density',"profiles",path_ql + data_ql)

print('pressure',p_full.shape,w_full.shape,p_mean_full.shape,w_mean_full.shape)

pw_full = covariance(p_full,w_full,p_mean_full,w_mean_full,tt)
#pw_ql = covariance(p_ql,w_ql,p_mean_ql,w_mean_ql,tt)

grad_full = gradient(pw_full,dz)
#grad_ql = gradient(pw_ql,dz)

P_full = np.zeros(shape = heat_flux_full.shape)
#P_ql = np.zeros(shape = heat_flux_ql.shape)
P_full[:] = -(1/density[:]*grad_full[:])


zvector = dz*np.linspace(0,nz-1,nz)
plt.figure(1,figsize=(20,10))
plt.subplot(1,2,1)
plt.plot(e_full,zvector)
plt.plot(u2_full)
plt.title('TKE')
plt.subplot(1,2,2)
plt.plot(S_full[0:nz],zvector,label='S full')
plt.plot(B_full[0:nz],zvector,label='B full')
#plt.plot(P_full[0:nz],zvector,label='P full')
plt.plot(T_full[0:nz],zvector,label='T full')
plt.plot(R_full[0:nz],zvector,label='T full')
plt.legend()
plt.savefig('pressure_' + T + '.png')
plt.savefig(path_out + 'pressure_' + T + '.png')
plt.close()

