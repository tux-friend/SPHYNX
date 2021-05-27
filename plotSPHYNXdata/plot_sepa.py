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
c=29979245800.0
R=1286311.1203265896
M=2.7837612956656893e33
tD=np.sqrt(R**3/(G*M))
Ainit=4*R
t0=5/256*c**5/(G**3)*(Ainit**4)/(2*M**3)



sepa = np.loadtxt('separation')
time = sepa[:,0]/tD
Ao=sepa[:,1]/Ainit
cm1x =sepa[:,2]
cm1y =sepa[:,3]
cm2x =sepa[:,4]
cm2y =sepa[:,5]

theoA=(1-time*tD/t0)**0.25

plt.figure(0)
plt.grid(linewidth=0.2)
plt.ylim([0.2,1.05])
plt.xlim([0, 155])
plt.plot(time,Ao,label='Simulation')
plt.plot(time,theoA,label='Point-mass approx.')
plt.hlines(0.61,0,155,linestyles=(0, (5, 5)),linewidth=0.7,colors='black')
plt.legend()
plt.ylabel(r'Orbital separation $A_o/A_{o,init}$')
plt.xlabel(r'Time $t/t_{D}$')
plt.savefig('plot_sepa.png',dpi=300,bbox_inches='tight')
plt.show()
plt.close()

plt.figure(1)
plt.grid(linewidth=0.2)
plt.plot(cm1x,cm1y,'-',label=r'Center of mass (NS$_1$)',linewidth=0.4)
plt.plot(cm2x,cm2y,'-',label=r'Center of mass (NS$_2$)',linewidth=0.4)
plt.legend()
plt.xlabel(r'$x$ [cm]')
plt.ylabel(r'$y$ [cm]')
plt.gca().set_aspect('equal',adjustable='box')
plt.savefig('plot_cm_merge.png',dpi=300,bbox_inches='tight')
plt.show()


