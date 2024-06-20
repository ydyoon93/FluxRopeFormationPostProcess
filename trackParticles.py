#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 15:02:26 2023

@author: young
"""

import happi
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.collections import LineCollection

S = happi.Open('PIC/BennettHeatingGuide10')
Lx = S.namelist.Lx
Ly = S.namelist.Ly
Nx = S.namelist.Nx
dx = S.namelist.dx
dy = S.namelist.dy
dt = S.namelist.dt
mi = S.namelist.mass_ratio
l = S.namelist.l*2
Bg = S.namelist.Bg
B0 = S.namelist.B0
omega_ratio = S.namelist.omega_ratio
vA = omega_ratio / mi ** .5
v_ref = l * omega_ratio * mi**.5
ions = S.TrackParticles(species = "test_ions",axes=['x','y','px','py','pz']).getData()
eons = S.TrackParticles(species = "test_eons",axes=['x','y','px','py','pz'],
                        # select='all(t<1,',
                        ).getData()

def add_r(Diag):
    Diag['r'] = ((Diag['x'] - Lx/2) ** 2 + (Diag['y'] - Ly/2) ** 2) ** .5
def add_phi(Diag):
    Diag['phi'] = np.arctan2((Diag['y'] - Ly/2),(Diag['x'] - Lx/2))
add_r(ions), add_r(eons)
add_phi(ions), add_phi(eons)

def add_pr(Diag):
    Diag['pr'] = Diag['px'] * np.cos(Diag['phi']) + Diag['py'] * np.sin(Diag['phi'])
def add_pphi(Diag):
    Diag['pphi'] = (- Diag['px'] * np.sin(Diag['phi']) + Diag['py'] * np.cos(Diag['phi']))
add_pr(ions), add_pr(eons)
add_pphi(ions), add_pphi(eons)
ions['mass'] = mi
eons['mass'] = 1
ions['q'] = 1
eons['q'] = -1
r_range = np.linspace(0.0001,1,Nx*10+1)
def gamma(pr,pphi,pz):
    return()
# @np.vectorize
def vphi(r,pphi,bg):
    return pphi/r*l + bg*r * B0
# @np.vectorize
def vz(r,pz):
    return pz + np.log(1+r**2)*B0*l/2
@np.vectorize
def vr(r,pz,pphi,bg,H):
    return (2 * H - vphi(r,pphi,bg)**2 - vz(r,pz)**2)**.5
#%% plot
t_fraction = 0.3 #fraction of the total time to plot
def phase_space_plot(ax,track,p_direction,par_num):
    # set up a list of (x,y) points
    t   = track["times"][:int(len(track["times"])*t_fraction)]*dt
    x   = track["r"][:len(t),par_num]/l
    g   = np.sqrt(1+track['px'][:len(t),par_num]**2 
                  + track['py'][:len(t),par_num]**2 
                  + track['pz'][:len(t),par_num]**2)
    y   = track[p_direction][:len(t),par_num]/track['mass']/v_ref*track['q']
    # x   = track["x"][:len(t),par_num]
    # y   = track["y"][:len(t),par_num]
    norm = plt.Normalize(t[0],t[-1])
    points = np.array([x,y]).transpose().reshape(-1,1,2)
    
    # set up a list of segments
    segs = np.concatenate([points[:-1],points[1:]],axis=1)
    
    # make the collection of segments
    lc = LineCollection(segs, cmap=plt.get_cmap('viridis'),norm=norm,alpha=0.2)
    lc.set_array(t) # color the segments by our parameter
    
    # plot the collection
    ax.add_collection(lc) # add the collection to the plot
    ax.scatter(x[0],y[0],c='C0')
    return lc

fig, ax = plt.subplots(1,3,sharey='row',sharex=1,figsize=[5,1.7])
for row_num,n in enumerate([40]):
    for i,p in enumerate(['pr','pphi','pz']):
        phase_space_plot(ax[i],eons,p,n)
    pz = eons['pz'][0,n] + np.log(1+eons['r'][0,n]**2/l**2)*B0*l/2
    pphi = eons['pphi'][0,n]*eons['r'][0,n]/l**2 - Bg * eons['r'][0,n]**2/2/l**2
    H  = 0.5 * (eons['px'][0,n]**2+eons['py'][0,n]**2+eons['pz'][0,n]**2)
    ax[0].plot(r_range,vr(r_range,-pz,-pphi,-Bg/B0,H)/v_ref,c='r')
    ax[0].plot(r_range,-vr(r_range,-pz,-pphi,-Bg/B0,H)/v_ref,c='r')
    ax[1].plot(r_range,vphi(r_range,-pphi,-Bg/B0)/v_ref,c='r')
    ax[2].plot(r_range,vz(r_range,-pz)/v_ref,c='r')
    ax[0].set_ylabel('$v/\lambda \omega_{ce}$')
ax[0].set_xlim([0,0.2])
ax[0].set_ylim([-0.05,0.05])
[ax[i].set_xlabel('$r/\lambda$') for i in range(3)]
# [ax.set_xlabel('$x/d_i$') for ax in ax.flatten()]
[ax.text(0.08,0.85,string,transform=ax.transAxes) for ax,string in zip(ax.flatten(),['(a) $(r,v_r)$','(b) $(r,v_\phi)$','(c) $(r,v_z)$',
                                                                                   '(d) $(r,v_r)$','(e) $(r,v_\phi)$','(f) $(r,v_z)$'])]
cax = fig.add_axes((0.13,0.92,0.77,0.05))
cax.set_yticks([])
cax.tick_params(axis='x',direction='in')

plt.colorbar(cm.ScalarMappable(colors.Normalize(0,eons["times"][int(len(eons["times"])*t_fraction)]*dt),
                                cmap='viridis'),
              cax=cax,
              orientation = 'horizontal')
cax.xaxis.tick_top()
# plt.savefig('eon_trajectory.pdf',bbox_inches='tight',dpi=300)
