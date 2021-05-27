from scipy.io import FortranFile
import os
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

model=10000

data = np.loadtxt('s1.{0:06d}'.format(model))
c1=data[:,1]
c2=data[:,2]
c3=data[:,3]
xc=np.sum(c1)/len(c1)
yc=np.sum(c2)/len(c2)
zc=np.sum(c3)/len(c3)
r = ((c1-xc)**2+(c2-yc)**2+(c3-zc)**2)**0.5
density=data[:,6]
smooth=data[:,4]
neighbors=data[:,27]
c16=data[:,15]
c17=data[:,16]
c18=data[:,17]
c19=data[:,18]
c20=data[:,19]
c21=data[:,20]
G=6.674*1e-8
press = ((c1-xc)*c19+(c2-yc)*c20+(c3-zc)*c21)/r
grav = G*((c1-xc)*c16+(c2-yc)*c17+(c3-zc)*c18)/r
 
#Density profile plot
plt.figure(0)
plt.plot(r,density,',',label=str(model)+' iterations')
plt.ylim([0, 0.8*1e15])
plt.xlim([0, 1.4*1e6])
plt.grid()
plt.title('Density profile for '+str(model)+' iterations')
plt.xlabel(r"Radius $r$ [cm]")
plt.ylabel(r"Density $\rho$ [g/cm$^3$]")
plt.legend(loc='best',markerscale=12,labelcolor='linecolor')
plt.savefig('density_profile_'+str(model)+'_iterations.png',dpi=300,bbox_inches='tight')

#Smoothing length plot
plt.figure(1)
plt.plot(r,smooth,',',label=str(model)+' iterations')
plt.title(r'Smoothing length $h$ after '+str(model)+' iterations')
plt.grid(linewidth=0.5)
plt.xlim([0,1.4e6])
plt.xlabel(r"Radius $r$ [cm]")
plt.ylabel(r"Smoothing length $h$")
plt.legend(loc='best',markerscale=12,labelcolor='linecolor')
plt.savefig('h_'+str(model)+'_iterations.png',dpi=300,bbox_inches='tight')
    


#gravity and pressure gradient plot
plt.figure(2)
plt.plot(r,press,',',label='pressure gradient')
plt.plot(r,grav,',',label='gravity gradient')
plt.title('Gradients for model with '+str(model)+' iterations')
plt.xlabel(r"Radius $r$ [cm]")
plt.ylabel(r"Gradients [?]")
plt.legend(loc='best',markerscale=12,labelcolor='linecolor')
plt.savefig('pressure_gravity_gradient_model_'+str(model)+'.png',dpi=300,bbox_inches='tight')



