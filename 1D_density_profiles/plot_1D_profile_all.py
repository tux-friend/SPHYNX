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
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title

LS220T05 = np.loadtxt('LS220_T05_Newton')
LS220S1 = np.loadtxt('LS220_S1_Newton')
P05 = np.loadtxt('1D_density_profile_0.5.txt')
P07 = np.loadtxt('1D_density_profile_0.7.txt')
MPA1 = np.loadtxt('1D_density_profile_MPA1.txt')

profiles = [P05,P07,MPA1,LS220T05,LS220S1]
n = ['P-05','P-07','PWP-MPA1',r'LS220, $T=0.5$ MeV ($T$-slice)', r'LS220, $S=1.0\,k_b$/baryon ($S$-slice)']

fig, sub = plt.subplots()

for prof,n in zip(profiles,n):
    r = prof[:,0]
    rho = prof[:,1]
    print(n,np.max(rho))
    sub.plot(r,rho,label=n)

sub.grid(linewidth=0.2)
sub.legend(loc='best')
sub.set_xlim([0, 1.45e6])
sub.set_ylim([0, 0.75e15])
sub.set_xlabel(r"Radius $r$ [cm]")
sub.set_ylabel(r"Density $\rho$ [gcm$^{-3}$]")
plt.savefig('1D_profiles_all.png',dpi=300,bbox_inches='tight')
plt.show()


