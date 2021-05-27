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

iterations =[1,4600,9200,13800,18400,23000,27600,34400,37000] #Choose iterations to compare state of relaxation
a, b = 3,3 #rows (a) and columns (b) of subplots
lim =[-7.5e6,7.5e6] #x- and y-axes limits

fig, axes= plt.subplots(a,b,figsize=(5*b,5*a))
plt.set_cmap('plasma')


i = 0
densmin=20
densmax=0
for it in iterations:
    star=np.loadtxt("s1.{:06d}".format(it))
    density=np.log10(star[:,6])
    if np.min(density) < densmin:
        densmin=np.min(density)
    else:
        densmin=densmin
    if np.max(density) > densmax:
        densmax=np.max(density)
    else:
        densmax=densmax
    star=np.loadtxt("s2.{:06d}".format(it))
    density=np.log10(star[:,6])
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

G=6.674*10**-8
star=np.loadtxt("s1.000001")
M=np.sum(star[:,25])
x = star[:,1]-np.sum(star[:,1])/len(star[:,1])
y = star[:,2]-np.sum(star[:,2])/len(star[:,2])
z = star[:,3]-np.sum(star[:,3])/len(star[:,3])
R=np.max((x**2+y**2+z**2)**0.5)
tD = (R**3/(G*M))**0.5
print(M,R,tD)
estabil = np.loadtxt(str(Path(os.getcwd()).parent)+'/estabil.d')

for it in iterations:
    ax = axes.flatten()[i]
    star=np.loadtxt("s1.{:06d}".format(it)) #load and plot star 1
    x = star[:,1]
    y = star[:,2]
    density=np.log10(star[:,6])
    idx = density.argsort() #sort by density -> high density points are plotted last
    x, y, density = x[idx], y[idx], density[idx]
    
    star=np.loadtxt("s2.{:06d}".format(it)) #load and plot star 2
    x2 = star[:,1]
    y2 = star[:,2]
    density2=np.log10(star[:,6])
    idx = density2.argsort() #sort by density -> high density points are plotted last
    x2, y2, density2 = x2[idx], y2[idx], density2[idx]
    
    x =np.concatenate([x,x2])
    y =np.concatenate([y,y2]) 
    density = np.concatenate([density,density2])
    im = ax.scatter(x,y,c=density,marker='o', lw=0, s=(72./(1.2*fig.dpi))**2,norm=normalizer)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect('equal')
    t=estabil[it-1,0]/tD
    ax.set_title(r't = {:.2f} t$_D$'.format(t),loc='right',fontsize=10)
    ax.set(xlabel='x [cm]',ylabel='y [cm]')
    i += 1

#fig.subplots_adjust(right=0.8)
cbar = fig.colorbar(im, ax=axes.ravel().tolist(),ticks=np.linspace(densmin,densmax,num=5),format="%.2f",shrink=0.6)
cbar.set_label('hallo',rotation=90)
plt.savefig('initial_models.png',dpi=300,bbox_inches='tight')
plt.show()



