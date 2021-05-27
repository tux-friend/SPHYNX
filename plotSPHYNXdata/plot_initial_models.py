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
n05 = np.loadtxt('n05_initial')
n07 = np.loadtxt('n07_initial')
mpa1 = np.loadtxt('mpa1_initial')
tab = np.loadtxt('TAB_initial')

names = ['P-05','P-07','PWP-MPA1']
initial_models = [n05,n07,mpa1]

fig, axes= plt.subplots(1,3,figsize=(12,3))
plt.set_cmap('magma')

i = 0
densmin=20
densmax=0
for model in initial_models:
    density=np.log10(model[:,6])
    if np.min(density) < densmin:
        densmin=np.min(density)
    else:
        densmin=densmin
    if np.max(density) > densmax:
        densmax=np.max(density)
    else:
        densmax=densmax
normalizer=Normalize(densmin,densmax)
im=cm.ScalarMappable(norm=normalizer)

for model, name in zip(initial_models,names):
    ax = axes.flatten()[i]
    r1 = model[:,1]
    r2 = model[:,2]
    r3 = model[:,3]
    h=model[:,4]
    mass=model[:,25]
    radius=model[:,10]
    print(str(name),'M:',mass[0]*100000,'R:',np.max(radius),'Rho_c:',np.max(model[:,6]))
    density=np.log10(model[:,6])
    x=r1-np.sum(r1)/n
    y=r2-np.sum(r2)/n
    z=r3-np.sum(r2)/n
    x,y,z,density=x.flatten(),y.flatten(),z.flatten(),density.flatten()
    use_points=(abs(z)/h<=2)
    x,y,z,density=x[use_points],y[use_points],z[use_points],density[use_points]
    print(len(z))
    #idx = density.argsort() #sort by density -> high density points are plotted last
    #x, y, density = x[idx], y[idx], density[idx]
    im = ax.scatter(x,y,c=density,marker='o', lw=0, s=(72./(0.5*fig.dpi))**2,norm=normalizer)
    ax.set_aspect("equal")
    ax.axis('off')
    ax.set_title(name,loc='center',fontsize=12)
    i += 1

fig.subplots_adjust(right=0.8)
cbar = fig.colorbar(im, ax=axes.ravel().tolist(),ticks=np.linspace(densmin,densmax,num=5),format="%.2f",shrink=1.0)
cbar.set_label(r'log($\rho$ [gcm$^{-3}$])',rotation=90)
plt.savefig('initial_models.png',dpi=300,bbox_inches='tight')
plt.show()



