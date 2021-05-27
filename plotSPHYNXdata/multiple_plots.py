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

os.mkdir('plots')

n=100000   #Number of particles
G=6.674*1e-8
fortran=FortranFile('neutron_star','r')

nsarray=[]
for i in range(n):
    nsarray.append(fortran.read_reals(dtype=float))
np.savetxt('neutron_star.txt', nsarray)
nsarray = np.loadtxt('neutron_star.txt')

c1 = nsarray[:, 0]
c2 = nsarray[:, 1]
c3 = nsarray[:, 2]
c4 = nsarray[:, 3]
c5 = nsarray[:, 4]
x = ((c1*c1+c2*c2+c3*c3)**0.5)
y = nsarray[:, 5]

for i in range(1,100):
    data = np.loadtxt("s1.0{0:03d}00".format(i))
    c1=data[:,1]
    c2=data[:,2]
    c3=data[:,3]
    xc=np.sum(c1)/n
    yc=np.sum(c2)/n
    zc=np.sum(c3)/n
    r = ((c1-xc)**2+(c2-yc)**2+(c3-zc)**2)**0.5
    density=data[:,6]
    c16=data[:,15]
    c17=data[:,16]
    c18=data[:,17]
    c19=data[:,18]
    c20=data[:,19]
    c21=data[:,20]
    press = ((c1-xc)*c19+(c2-yc)*c20+(c3-zc)*c21)/r
    grav = G*((c1-xc)*c16+(c2-yc)*c17+(c3-zc)*c18)/r
    plt.figure()
    plt.grid()
    plt.plot(r,density,',',label=str(i*100)+' iterations')
    plt.plot(x,y,',',label='initial model')
    plt.ylim([0, 8*1e14])
    plt.xlim([0, 1.4*1e6])
    plt.title('Density profile for '+str(i*100)+' iterations')
    plt.xlabel(r"Radius $r$ [cm]")
    plt.ylabel(r"Density $\rho$ [g/cm$^3$]")
    plt.legend(loc='best',markerscale=12,labelcolor='linecolor')
    plt.savefig('plots/density_profile_'+str(i*100)+'_iterations.png',dpi=300)
    plt.close()
    plt.figure()
    plt.plot(r,press,',',label='pressure gradient')
    plt.plot(r,grav,',',label='gravity gradient')
    plt.title('Gradients for model with '+str(i*100)+' iterations')
    plt.xlabel(r"Radius $r$ [cm]")
    plt.ylabel(r"Gradients [?]")
    plt.legend(loc='best',markerscale=12,labelcolor='linecolor')
    plt.savefig('plots/pressure_gravity_gradient_model_'+str(i*100)+'.png',dpi=300)
    plt.close()
   






