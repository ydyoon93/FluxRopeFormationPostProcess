#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 12:33:46 2023

@author: young
"""

import numpy as np
import matplotlib.pyplot as plt
import happi

S = happi.Open("PIC/BennettHeatingGuide10")
[Bx,By,Bz,Bt] = [S.Field(0,B) for B in ["Bx","By","Bz",
                                        "(Bx**2 + By**2 + Bz**2)**.5"]]
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
d = {}
for i,direction in enumerate(['x','y','z']):
    d['x_v{0}_i'.format(direction)] = S.ParticleBinning(i)
    d['x_v{0}_e'.format(direction)] = S.ParticleBinning(i+12)

tsteps = d['x_vx_i'].getTimesteps()
times = d['x_vx_i'].getTimes()
v_i = d['x_vx_i'].getAxis('vx')
x_i = d['x_vx_i'].getAxis('x')
v_e = d['x_vx_e'].getAxis('vx')/v_ref
x_e = d['x_vx_e'].getAxis('x')/l - Lx/l/2


dd = {}
for keys in d:
    dd[keys] = [(d[keys].getData(tsteps[t])[0]-
                 d[keys].getData(tsteps[0])[0]).transpose()
                for t in range(tsteps.size)]
#%% Plot 
fig,ax = plt.subplots(1,3,figsize=[7,2])
for axes,diag in zip(ax.flatten(),['x_vx_e','x_vy_e','x_vz_e']):
    # p = axes.imshow(dd[diag][-1],extent=[x_e[0],x_e[-1],v_e[0],v_e[-1]],
    #               aspect='auto',origin='upper',cmap='PiYG',
    #               vmax=35,vmin=-35,
    #               )
    p = axes.contourf(x_e,v_e,np.flipud(dd[diag][-1]),
                  # aspect='auto',
                  # origin='upper',
                  cmap='PuOr_r',
                  levels = np.linspace(-35,35,15,endpoint=1)
                  )
    axes.contour(x_e,v_e,np.flipud(d[diag].getData(0)[0].transpose()),
                 colors='k',
                 levels=np.linspace(5,60,6),
                 linewidths=0.1,
                 )
plt.colorbar(p,ax=ax)

[ax.set_xlim([0,1]) for ax in ax.flatten()]
[ax.set_ylim([-0.045,0.045]) for ax in ax.flatten()]
[ax.set_xlabel(r'$r/\lambda$') for ax in ax.flatten()]
[ax.text(0.1,0.85,title,transform=ax.transAxes) for ax,title in zip(ax.flatten(),
                                            ['(a) $\Delta f_e(r,v_r)$',
                                             '(b) $\Delta f_e(r,v_\phi)$',
                                             '(c) $\Delta f_e(r,v_z)$'])]
ax[1].set_yticks([])
ax[2].set_yticks([])
ax[0].set_ylabel('$v/\lambda\omega_{ce}$')
# plt.savefig('manuscript/phase_space_e.pdf',bbox_inches='tight',dpi=300)
# plt.savefig('phase_space_e.pdf',bbox_inches='tight',dpi=300)