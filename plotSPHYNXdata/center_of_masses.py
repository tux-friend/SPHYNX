from scipy.io import FortranFile
import os
import matplotlib
import numpy  as np
import matplotlib.pyplot as plt
SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title


n=100000   #Number of particles
G=6.674*1e-8
cm=np.zeros((1,7))

for i in [1]+list(range(100,5000,100)):
    ns1 = np.loadtxt("s1.{:06d}".format(i))
    ns2 = np.loadtxt("s2.{:06d}".format(i))
    c1=ns1[:,1]
    c2=ns1[:,2]
    c3=ns1[:,3]
    xc=np.sum(c1)/n
    yc=np.sum(c2)/n
    zc=np.sum(c3)/n
    c1=ns2[:,1]
    c2=ns2[:,2]
    c3=ns2[:,3]
    xc2=np.sum(c1)/n
    yc2=np.sum(c2)/n
    zc2=np.sum(c3)/n
    cm=np.append(cm,[[int(i),xc,yc,zc,xc2,yc2,zc2]],axis=0)

theta = np.linspace(-np.pi, np.pi, 200)    
r=6.753132875*1e6
x1=r*np.sin(theta)
x2=r*np.cos(theta)  
cm=np.delete(cm,(0),axis=0)    
np.savetxt('center_of_masses.txt',cm)   

plt.figure()
plt.grid(linewidth=0.2)
plt.ylim([-1e7, 1e7])
plt.xlim([-1e7, 1e7])
plt.plot(x1,x2,linewidth=0.3)
plt.plot(cm[:,1],cm[:,2],'.',label=r'Center of mass (NS$_1$)')
plt.plot(cm[:,4],cm[:,5],'.',label=r'Center of mass (NS$_2$)')
plt.legend()
plt.xlabel(r'$x$ [cm]')
plt.ylabel(r'$y$ [cm]')
plt.gca().set_aspect('equal',adjustable='box')
plt.savefig('plot_CM.png',dpi=300,bbox_inches='tight')
plt.close()

