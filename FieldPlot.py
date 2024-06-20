#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 12:46:21 2023

@author: young
"""
import happi
import numpy as np
import matplotlib.pyplot as plt

S = happi.Open('PIC/BennettHeatingGuide10')
Lx = S.namelist.Lx
Ly = S.namelist.Ly
Nx = S.namelist.Nx
Ny = S.namelist.Ny
dx = S.namelist.dx
dy = S.namelist.dy
B0 = S.namelist.B0
n0 = S.namelist.n0
x_range = np.linspace(0,Lx,Nx+1)
y_range = np.linspace(0,Ly,Ny+1)
tsteps = S.Field(0,'Bx').getTimesteps()
fields = {}
for field in ['Bx','By','Bz',
              'Jx_ions+Jx_eons','Jy_ions+Jy_eons','Jz_ions+Jz_eons',
              '-Rho_eons']:
    fields[field] = S.Field(0,field)
    
P   = S.ParticleBinning('#6+#7+#8+#18+#19+#20')
Azi = -np.cumsum(fields['By'].getData(tsteps[0])[0],axis=0)*dx
f   = np.cumsum(fields['Bx'].getData(tsteps[0])[0][0],axis=0)*dy
f   = np.tile(f,(int(Nx+1),1))
Azi += f
Azf =-np.cumsum(fields['By'].getData(tsteps[-1])[0],axis=0)*dx
f   = np.cumsum(fields['Bx'].getData(tsteps[-1])[0][0],axis=0)*dy
f   = np.tile(f,(int(Nx+1),1))
Azf += f
#%% Plot
fig,ax = plt.subplots(4,3,figsize=[6,6],sharex=True,sharey=True)
p = ax[0,0].imshow(fields['Bz'].getData(tsteps[0])[0]/B0,
               extent = [x_range[0],x_range[-1],y_range[0],y_range[-1]],
               vmax=0.7,vmin=0.1)
plt.colorbar(p)
ax[0,0].contour(x_range,y_range,Azi,colors='w',linewidths=0.5,levels=20)
arrow_xrange = np.linspace(x_range - Lx/2 - 0.7, x_range - Lx/2 + 0.7,100)
plt.arrow(Lx/2 + 0.7, Lx/2 + 0.7, -1.4, -1.4, color='r',head_width=.1)

p = ax[1,0].imshow(fields['Bz'].getData(tsteps[-1])[0]/B0,
               extent = [x_range[0],x_range[-1],y_range[0],y_range[-1]],
               vmax=0.7,vmin=0.1)
plt.colorbar(p)
ax[1,0].contour(x_range,y_range,Azf,colors='w',linewidths=0.5,levels=20)

p = ax[0,1].imshow(fields['Jz_ions+Jz_eons'].getData(tsteps[0])[0]/n0,
               extent = [x_range[0],x_range[-1],y_range[0],y_range[-1]],
               vmax=.4,
               cmap = 'cividis')
plt.colorbar(p)
# ax[1,0].contour(x_range,y_range,fields['Bz'],colors='w')
p = ax[1,1].imshow(fields['Jz_ions+Jz_eons'].getData(tsteps[-1])[0]/n0,
               extent = [x_range[0],x_range[-1],y_range[0],y_range[-1]],
               vmax=.4,
               cmap = 'cividis')
plt.colorbar(p)
ax[1,1].contour(x_range,y_range,fields['Bz'].getData(tsteps[-1])[0],
                colors='w',
                linewidths=0.5,
                levels=np.linspace(4,7,7)
                )
extent = 27
xlim = slice(int(Nx/2)-extent,int(Nx/2)+extent)
Range = np.arange(int(Nx/2)-extent,int(Nx/2)+extent,1)*dx
ax[1,1].contour(Range,Range,fields['Bz'].getData(tsteps[-1])[0][xlim,xlim],
                colors='k',
                linewidths=0.7,
                levels=np.linspace(4,7,7)
                )

p = ax[0,2].imshow(P.getData(tsteps[0])[0]/n0/3,
               extent = [x_range[0],x_range[-1],y_range[0],y_range[-1]],
               vmax=1.5,
               cmap = 'hot')
plt.colorbar(p)
p = ax[1,2].imshow(P.getData(tsteps[-1])[0]/n0/3,
                extent = [x_range[0],x_range[-1],y_range[0],y_range[-1]],
                vmax=1.5,
               cmap = 'hot')
plt.colorbar(p)

for axes,binNum in zip([ax[2,0],ax[2,1],ax[2,2]],
                       [18,19,20]):
    p = axes.imshow(S.ParticleBinning(binNum).getData(tsteps[-1])[0]/n0,
                    extent = [x_range[0],x_range[-1],y_range[0],y_range[-1]],
                    vmax=1,
                    cmap = 'hot')
    plt.colorbar(p)
for axes,binNum in zip([ax[3,0],ax[3,1],ax[3,2]],
                       [21,22,23]):
    p = axes.imshow(S.ParticleBinning(binNum).getData(tsteps[-1])[0]/n0,
                    extent = [x_range[0],x_range[-1],y_range[0],y_range[-1]],
                    vmin=-0.07,
                    vmax=0.07,
                    cmap = 'coolwarm')
    plt.colorbar(p)

[ax[3,i].set_xlabel('$x/d_i$') for i in range(3)]
# [[ax[j,i].set_xticks([]) for i in range(2)]for j in range(2)]
[ax[i,0].set_ylabel('$y/d_i$') for i in range(4)]
# [ax[i,1].set_yticks([]) for i in range(3)]

for ax,string in zip(ax.flatten(),
                     ['(a) $\mathbf{B}_i$','(b) $\mathbf{J}_i$','(c) $P_{i}$',
                      '(d) $\mathbf{B}_f$','(e) $\mathbf{J}_f$','(f) $P_{f}$',
                      '(g) $P_{exx}$','(h) $P_{eyy}$','(i) $P_{ezz}$',
                      '(j) $P_{exy}$','(k) $P_{eyz}$','(l) $P_{ezx}$']):
    ax.set_xlim([4,6])
    ax.set_ylim([4,6])
    ax.text(0.1,0.8,string,transform=ax.transAxes,c='w')

# plt.savefig('fields.pdf',bbox_inches='tight',dpi=300)
# plt.savefig('manuscript/fields.pdf',bbox_inches='tight',dpi=300)
