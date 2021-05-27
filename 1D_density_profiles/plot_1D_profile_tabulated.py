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

profiles = [LS220T05,LS220S1]
n = [r'$T=0.5$ MeV (temperature slice)', r'$S=1.0\,k_b$/baryon (entropy slice)']

fig, sub = plt.subplots()

for prof,n in zip(profiles,n):
    r = prof[:,0]
    rho = prof[:,1]
    sub.plot(r,rho,label=n)

sub.grid()
sub.legend(loc='best')
sub.set_xlim([0, 1.5e6])
sub.set_ylim([0, 0.5e15])
sub.set_xlabel(r"Radius $r$ [cm]")
sub.set_ylabel(r"Density $\rho$ [gcm$^{-3}$]")
plt.title('1D density profile for LS220')
plt.savefig('density_profiles_tabEoS.png',dpi=300,bbox_inches='tight')
plt.show()


fig, sub = plt.subplots()
n = [r'$T=0.5$ MeV (temperature slice)', r'$S=1.0\,k_b$/baryon (entropy slice)']
for prof,n in zip(profiles,n):
    r = prof[:19200,0]
    T = prof[:19200,4]
    sub.plot(r,T,label=n)

sub.grid()
sub.legend(loc='best')
sub.set_xlim([0, 1.5e6])
sub.set_ylim([0, 14.0])
sub.set_xlabel(r"Radius $r$ [cm]")
sub.set_ylabel(r"Temperature $T$ [MeV]")
plt.title('1D temperature profile for LS220')
plt.savefig('temperature_profiles_tabEoS.png',dpi=300,bbox_inches='tight')
plt.show()

fig, sub = plt.subplots()
n = [r'$T=0.5$ MeV (temperature slice)', r'$S=1.0\,k_b$/baryon (entropy slice)']
for prof,n in zip(profiles,n):
    r = prof[:,0]
    Ye = prof[:,5]
    sub.plot(r,Ye,label=n)

sub.grid()
sub.legend(loc='best')
sub.set_xlim([0, 1.5e6])
sub.set_ylim([0, 0.55])
sub.set_xlabel(r"Radius $r$ [cm]")
sub.set_ylabel(r"Electron fraction $Y_e$")
plt.title('1D electron fraction profile for LS220')
plt.savefig('ye_profiles_tabEoS.png',dpi=300,bbox_inches='tight')
plt.show()
