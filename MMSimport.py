#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 16:28:43 2021

@author: young
"""


import pyspedas
import pytplot
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interpn
import happi
time_range = ['2017-07-06/08:34:09.5','2017-07-06/08:34:13.5'] # Sun2019
# time_range = ['2015-10-16/13:04:31','2015-10-16/13:04:37'] # Zhao2016_1
# time_range = ['2015-12-14/00:58:55','2015-12-14/00:59:03'] # Zhao2016_2


fgm = pyspedas.mms.fgm(trange=time_range,
                       data_rate='brst',
                       probe = ['1','2','3','4'],
                       time_clip = True,
                       # notplot = True,
                       )
edp = pyspedas.mms.edp(trange=time_range,
                        data_rate='brst',
                        probe = ['1','2','3','4'],
                        time_clip = True,
                        # notplot = True,
                        )
def get_fpi_data(datatype):
    return pyspedas.mms.fpi(trange=time_range, 
                            datatype=datatype,
                            data_rate='brst',
                            probe = ['1','2','3','4'],
                            time_clip = True,
                            # notplot = True,
                            )
fpi_vi = get_fpi_data('dis-moms')
fpi_ve = get_fpi_data('des-moms')
# fpi_i = get_fpi_data('dis-dist')
# fpi_e = get_fpi_data('des-dist')
    
# positions = ['mms1_fgm_r_gse_brst_l2',
#              'mms2_fgm_r_gse_brst_l2',
#              'mms3_fgm_r_gse_brst_l2',
#              'mms4_fgm_r_gse_brst_l2']
# fields = ['mms1_fgm_b_gse_brst_l2',
#           'mms2_fgm_b_gse_brst_l2',
#           'mms3_fgm_b_gse_brst_l2',
#           'mms4_fgm_b_gse_brst_l2']

vfields = ['mms1_des_bulkv_gse_brst',
           'mms2_des_bulkv_gse_brst',
           'mms3_des_bulkv_gse_brst',
           'mms4_des_bulkv_gse_brst']
Pfields = ['mms1_des_temptensor_gse_brst',
           'mms2_des_temptensor_gse_brst',
           'mms3_des_temptensor_gse_brst',
           'mms4_des_temptensor_gse_brst']
#### convert to gsm
[pyspedas.cotrans(vfield,vfield.replace('gse','gsm'),
                  coord_in='gse',coord_out='gsm') for vfield in vfields]

[pyspedas.cotrans(vfield,vfield.replace('gse','gsm'),
                  coord_in='gse',coord_out='gsm') for vfield in Pfields]
# def get_tensor_and_store_vector(string):
#     for i in range(4):
#         data = pytplot.get_data('mms' + str(i+1) + string)
#         for j in range(3):
#             pytplot.store_data('mms' + str(i+1) + string + '_vec_' + str(j+1),
#                                 data = {'x':data[0],'y':data[1][...,j]})
# get_tensor_and_store_vector('_des_prestensor_gse_brst')
#%%Rearrange data into LMN coords
def average_spacecraft(string):
    return np.sum(np.array([pytplot.data_quants['mms'+probe+string]
                   for probe in ['1']]),0)
# J = pyspedas.mms.curlometer(fields = fields,
#                             positions = positions)
# J = np.array(pytplot.data_quants['curlB'])

B = average_spacecraft('_fgm_b_gsm_brst_l2')
# B_gsm = average_spacecraft('_fgm_b_gsm_brst_l2')
n_e = average_spacecraft('_des_numberdensity_brst')
T_e = average_spacecraft('_des_temptensor_gse_brst')
# T_par = average_spacecraft('_des_temppara_brst')
# T_perp = average_spacecraft('_des_tempperp_brst')
# T_i_par = average_spacecraft('_dis_temppara_brst')
# T_i_perp = average_spacecraft('_dis_tempperp_brst')
v_e = average_spacecraft('_des_bulkv_gsm_brst')
v_e -= np.array([811,-24,-61])
# v_i = average_spacecraft('_dis_bulkv_gse_brst')
# f_i = average_spacecraft('_dis_dist_brst')
# f_e = average_spacecraft('_des_dist_brst')
# P = average_spacecraft('_des_prestensor_gsm_brst')
E = average_spacecraft('_edp_dce_gse_brst_l2')
time = np.array(pytplot.data_quants['mms1_fgm_b_gsm_brst_l2'].coords['time'])
# time_i = np.array(pytplot.data_quants['mms2_dis_bulkv_gsm_brst'].coords['time'])
time_e = np.array(pytplot.data_quants['mms1_des_bulkv_gsm_brst'].coords['time'])
time_E = np.array(pytplot.data_quants['mms2_edp_dce_gse_brst_l2'].coords['time'])

##### define rotation matrix function (Rodriguez rotation formula)
def rot_mat(theta,vec):
    W = np.array([[0,-vec[2],vec[1]],
                  [vec[2],0,-vec[0]],
                  [-vec[1],vec[0],0]])
    return np.identity(3) + np.sin(theta) * W + 2 * np.sin(theta/2)**2 * W @ W
    # c = np.cos(theta)
    # s = np.sin(theta)
    
#### Determine GSM rotation angle
# B_tot = (B[:,1]**2+B[:,2]**2)**.5
# unit_B = [B[:,1],B[:,2]]/B_tot
# unit_B_gsm = [B_gsm[:,1],B_gsm[:,2]]/B_tot
                  


##### Determine LMN coords ######
# variance_mat = np.zeros((3,3))
# for i in range(3):
#     for j in range(3):
#         variance_mat[i,j] = np.mean(B[:,i]*B[:,j])-np.mean(B[:,i])*np.mean(B[:,j])
# trans_mat = np.linalg.eig(variance_mat)[1]
# trans_mat[:,0] = -trans_mat[:,0]
# trans_mat[:,1] = -trans_mat[:,1]

# Transform matrix of Sun2019
trans_mat = np.array([[-0.96,0.291,-0.0042],
                      [-0.29,-0.95,0.12],
                      [0.03,0.116,0.99]])
trans_mat = trans_mat[:,[2,0,1]]


# Zhao2016_1
# trans_mat = np.array([[-0.259,0.895,-0.362],
#                       [-0.146,0.335,0.931],
#                       [0.955,0.294,0.044]]).transpose()
# trans_mat = trans_mat[:,[2,0,1]]

# Zhao2016_2
# trans_mat = np.array([[-0.555,-0.695,-0.458],
#                       [-0.562,0.719,-0.410],
#                       [0.614,0.030,-0.789]]).transpose()
# trans_mat = trans_mat[:,[2,1,0]]


trans_mat = rot_mat(np.pi-np.pi/8,trans_mat[:,2]) @ trans_mat


# #### do LMN transforms of variables ####
def LMN_transform(var):
    for t in range(np.size(var,0)):
        var[t] = trans_mat.transpose() @ var[t]
    return var
for var in [v_e]:
    var = LMN_transform(var)
for t in range(np.size(B,0)):
    B[t,:-1] = trans_mat.transpose() @ B[t,:-1]
for t in range(np.size(B,0)):
    E[t] = trans_mat.transpose() @ E[t]
T_e = trans_mat.transpose() @ T_e @ trans_mat
# P = trans_mat.transpose() @ P @ trans_mat
# pytplot.store_data('temp',data={'x':time_e,
#                                 'y':T_e})
# pytplot.store_data('B',data={'x':time,'y':B})
# pytplot.store_data('v_e',data={'x':time_e,
#                                 'y':v_e})
# pytplot.store_data('J',data={'x':time,
#                               'y':J})
# B_interp = interp1d(time,B.transpose())(time_e).transpose()
# E_interp = interp1d(time_E,E.transpose())(time_e).transpose()
#%% Plot
fig,ax = plt.subplots(5,2,sharex='col',figsize=[7,6])
ax[0,0].plot(time,B)
ax[1,0].plot(time_e, -n_e[:,np.newaxis] *v_e)
# ax[1,0].plot(time_E[:-1],E)
# ax[1,0].plot(time_e, -n_e * (v_e[:,1]+800),c='C1',alpha=0.3)
# ax[1,0].plot(time_i,v_i)
ax[2,0].plot(time_e,n_e)
[ax[3,0].plot(time_e,T) for T in [T_e[:,0,0],T_e[:,1,1],T_e[:,2,2]]]
[ax[4,0].plot(time_e,T) for T in [T_e[:,0,1],T_e[:,1,2],T_e[:,0,2]]]
ax[4,0].set_ylim([-100,100])
ax[4,0].text(0.9,-0.23,'[sec]',transform=ax[4,0].transAxes)
[ax.set_ylabel(label) for ax,label in zip(ax[:,0],
                                          ['[nT]',
                                           '[q$_e$ km/s cm$^{-3}$]',
                                           '[cm$^{-3}$]',
                                           '[eV]',
                                           '[eV]'])]

S = happi.Open("../PIC/BennettHeatingGuide11")
[Bx,By,Bz,Bt] = [S.Field(0,B) for B in ["Bx","-By","-Bz",
                                        "(Bx**2 + By**2 + Bz**2)**.5"]]
[Jx,Jy,Jz] = [S.Field(0,J) for J in ["Jx_eons",
                                     "-Jy_eons",
                                     "-Jz_eons"]]
[Pxx,Pyy,Pzz,
 Pxy,Pyz,Pzx] = [S.ParticleBinning(P) for P in ['#18','#19','#20',
                                                '-#21','#22','-#23']]
Rho = S.Field(0,"-Rho_eons")
tsteps = Bx.getTimesteps()
Lx = S.namelist.Lx
Ly = S.namelist.Ly
Nx = S.namelist.Nx
Ny = S.namelist.Ny
dx = S.namelist.dx
dy = S.namelist.dy
B0 = S.namelist.B0
n0 = S.namelist.n0
mi = S.namelist.mass_ratio
vA = B0/(n0*mi)**.5

# Spacecraft trajectory
def y(x):
    return -x + Ly

x_range = np.linspace(Lx/2+.7,Lx/2-.7,1000)
y_range = y(x_range)
l_range = np.linspace(-((x_range[0]-Lx/2)**2 + (y_range[0]-Ly/2)**2)**.5,
                      ((x_range[-1]-Lx/2)**2 + (y_range[-1]-Ly/2)**2)**.5,
                      1000)

def interp(func,t):
    return interpn((np.linspace(0,Lx,Nx+1),np.linspace(0,Ly,Ny+1)),
                   func.getData(tsteps[t])[0],(x_range,y_range))
t = 99
[ax.yaxis.tick_right() for ax in ax[:,1]]
[ax.yaxis.set_label_position("right") for ax in ax[:,1]]
[ax[0,1].plot(l_range,interp(B,t)/B0) for B in [Bx,By,Bz,Bt]]
[ax[0,1].text(0.34+0.08*i,0.6,string,transform=ax[0,1].transAxes,c=color) 
 for i,[string,color] in enumerate(zip(['$B_x$','$B_y$','$B_z$','$B_t$'],
                                       ['C0','C1','C2','C3']))]
[ax[1,1].plot(l_range,interp(J,t)/n0/vA) for J in [Jx,Jy,Jz]]
[ax[1,1].text(0.1+0.08*i,0.2,string,transform=ax[1,1].transAxes,c=color) 
 for i,[string,color] in enumerate(zip(['$J_x$','$J_y$','$J_z$'],
                                       ['C0','C1','C2']))]
ax[2,1].plot(l_range,interp(Rho,t)/n0)
ax[2,1].text(0.5,0.5,'$n_e$',transform=ax[2,1].transAxes,c='C0')
[ax[3,1].plot(l_range,(interp(P,t)/interp(Rho,t))/vA**2) for P in [Pxx,Pyy,Pzz]]
[ax[3,1].text(0.7+0.08*i,0.1,string,transform=ax[3,1].transAxes,c=color) 
 for i,[string,color] in enumerate(zip(['$T_{xx}$','$T_{yy}$','$T_{zz}$'],
                                       ['C0','C1','C2']))]
[ax[4,1].plot(l_range,(interp(P,t)/interp(Rho,t))/vA**2) for P in [Pxy,Pyz,Pzx]]
ax[4,1].set_ylim([-7,7])
[ax[4,1].text(0.7+0.08*i,0.1,string,transform=ax[4,1].transAxes,c=color) 
 for i,[string,color] in enumerate(zip(['$T_{xy}$','$T_{yz}$','$T_{zx}$'],
                                       ['C0','C1','C2']))]
ax[4,1].set_xlabel('$l/d_i$')
[ax.set_ylabel(label) for ax,label in zip(ax[:,1],
                                          ['[$B_0$]',
                                           '[$n_0 q_e v_A$]',
                                           '[$n_0$]',
                                           '[$m_e v_A^2$]',
                                           '[$m_e v_A^2$]'])]
[ax.text(0.05,0.75,label,transform=ax.transAxes) 
 for ax,label in zip(ax.transpose().flat,['(a)','(b)','(c)','(d)','(e)',
                                          '(f)','(g)','(h)','(i)','(j)'])]
fig.tight_layout()
# plt.savefig('../manuscript/mms_pic_comparison.pdf',dpi=300,bbox_inches='tight')
