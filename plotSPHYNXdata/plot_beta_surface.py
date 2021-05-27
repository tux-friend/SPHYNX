import numpy  as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
import matplotlib
import os
import sys
from scipy.io import FortranFile
from matplotlib.gridspec import GridSpec
from mpl_toolkits import mplot3d
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition
from pathlib import Path

SMALL_SIZE = 9
MEDIUM_SIZE = 11
BIGGER_SIZE = 12

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title

n=100000   #Number of particles
betasurf = np.loadtxt('betasurface')
betaline = np.loadtxt('betaline')
x=betasurf[:,0]
y=betasurf[:,1]
dmu = betasurf[:,2]
x2=betaline[:,0]
y2=betaline[:,1]
minb=betaline[:,2]

zero=0*dmu
fig =plt.figure(0,figsize=(14,9))
ax=plt.axes(projection='3d')
ax.scatter(x,y,dmu,marker='o', s=(72./fig.dpi)**2,alpha=0.2,label=r'Table values of LS220 for $T=0.5$ MeV')
ax.scatter(x2,y2,minb,marker='o', s=(72./fig.dpi)**2,c='red',label=r'line with minimal $\mu_p+\mu_e-\mu_n\geqslant 0$')
ax.set_xlabel(r'log($\rho$ [gcm$^{-3}$])')
ax.set_ylabel(r'Electron fraction $Y_e$')
ax.set_zlabel(r'$\mu_p+\mu_e-\mu_n$')
ax.legend()
plt.savefig('betasurface.png',dpi=300,bbox_inches='tight')
plt.show()

plt.figure(1)
fig,ax1=plt.subplots()
ls220 = np.loadtxt('LS220_dmutot')
r=ls220[:,0]
ye=ls220[:,1]
dmutot=ls220[:,2]
plt.grid(linewidth=0.2)
plt.xlim([0,1.5e6])
ax1.plot(r,ye,'-',c='tab:blue',linewidth=1,label=r'LS220, $T=0.5$ MeV')
ax1.set_xlabel(r'Radius $r$ [cm]')
ax1.set_ylabel(r'Electron fraction $Y_e$',color='tab:blue')
ax1.tick_params('y',colors='tab:blue')

ax2=ax1.twinx()
ax2.plot(r,dmutot,'-',c='tab:red',linewidth=1,alpha=0.6,label=r'LS220, $T=0.5$ MeV')
ax2.set_ylabel(r'$\mu_p+\mu_e-\mu_n$',color='tab:red')
ax2.tick_params('y',colors='tab:red')
fig.tight_layout()
fig.legend(loc="upper left", bbox_to_anchor=(0,1), bbox_transform=ax1.transAxes)
plt.savefig('beta_ye_dmutot.png',dpi=300,bbox_inches='tight')
plt.show()

