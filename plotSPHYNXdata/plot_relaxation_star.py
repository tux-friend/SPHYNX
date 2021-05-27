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

n=100000   #Number of particles
iterations =[1,100,10000] #Choose iterations to compare state of relaxation
a, b = 1,3 #rows (a) and columns (b) of subplots

fig, axes= plt.subplots(a,b,figsize=(b*5,a*5))
plt.set_cmap('plasma')
i=0
for it in iterations:
    ax = axes.flatten()[i]
    star=np.loadtxt("s1.{:06d}".format(it)) #load iteration file
    r1 = star[:,1]
    r2 = star[:,2]
    x=r1-np.sum(r1)/n #get x coordinate of particles 
    y=r2-np.sum(r2)/n #get y coordinate of particles
    im = ax.scatter(x,y,marker='o', lw=0, s=(72./(1.2*fig.dpi))**2)
    axes.flatten()[i].set_aspect("equal")
    if it==1:
        a=' iteration'
    else:
        a=' iterations'
    ax.set_title(str(it)+a,loc='center',fontsize=12)
    ax.set(xlabel='x [cm]',ylabel='y [cm]')
    i += 1

plt.savefig('relaxation_star.png',dpi=300,bbox_inches='tight')
plt.show()



