import numpy  as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
import matplotlib
import os
import sys
from scipy.io import FortranFile
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition
from pathlib import Path
SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title

LS220T05 = np.loadtxt('LS220_T05_Newton')
LS220TOV = np.loadtxt('LS220_T05_TOV')
LS220Kaplan = np.loadtxt('LS220_Kaplan.dat')

profiles = [LS220T05,LS220TOV,LS220Kaplan]
n = [r'LS220 $T=0.5$ MeV (newTOV, Newtonian)',r'LS220 $T=0.5$ MeV (newTOV, TOV)', r'LS220 $T=0.5$ MeV (Kaplan, TOV)']

fig, sub = plt.subplots()
i=0
for prof,n in zip(profiles,n):
    r = prof[:,0]
    if i==2:
        line='dotted'
        rho = prof[:,2]
        sub.plot(r,rho,label=n,linestyle=line,color='black')
    else:
        line='solid'
        rho = prof[:,3]
        sub.plot(r,rho,label=n,linestyle=line)
    i=i+1

sub.grid()
sub.legend(loc='best')
plt.xlim([0,1.45e6])
plt.ylim([0,3.2e34])
sub.set_xlabel(r"Radius $r$ [cm]")
sub.set_ylabel(r"Pressure $p$ [dyne/cm$^{2}$]")
plt.savefig('check_newtov_P.png',dpi=300,bbox_inches='tight')
plt.show()

fig, sub = plt.subplots()
i=0
n = [r'LS220 $T=0.5$ MeV (newTOV, Newtonian)',r'LS220 $T=0.5$ MeV (newTOV, TOV)', r'LS220 $T=0.5$ MeV (Kaplan, TOV)']
for prof,n in zip(profiles,n):
    r = prof[:,0]
    if i==2:
        line='dotted'
        rho = prof[:,1]
        sub.plot(r,rho,label=n,linestyle=line,color='black')
    else:
        line='solid'
        rho = prof[:,1]
        sub.plot(r,rho,label=n,linestyle=line)
    i=i+1

sub.grid()
sub.legend(loc='best')
plt.xlim([0,1.45e6])
plt.ylim([0,0.5e15])
sub.set_xlabel(r"Radius $r$ [cm]")
sub.set_ylabel(r"Pressure $p$ [dyne/cm$^{2}$]")
plt.savefig('check_newtov_rho.png',dpi=300,bbox_inches='tight')
plt.show()

