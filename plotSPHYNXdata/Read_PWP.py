import numpy  as np
import matplotlib.pyplot as plt

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

c2=29979245800.0**2

Gamma = [1.35692, 3.446,3.572,2.887]
K0=3.99874e-8*c2
rho1=14.7
rho2=15.0
p1=34.495

K1=10**p1/((10**rho1)**Gamma[1])
K2=K1*(10**rho1)**Gamma[1]/((10**rho1)**Gamma[2])
K3=K2*(10**rho2)**Gamma[2]/((10**rho2)**Gamma[3])
rho0=np.log10((K1/K0)**(1/(Gamma[0]-Gamma[1])))
p0=np.log10(K0*(10**rho0)**Gamma[0])

x0=np.linspace(13.5,1.003*rho0,100)
x1=np.linspace(0.997*rho0,rho1,100)
x2=np.linspace(rho1,rho2,100)
x3=np.linspace(rho2,15.5,100)

y0=K0*(10**x0)**Gamma[0]
y1=K1*(10**x1)**Gamma[1]
y2=K2*(10**x2)**Gamma[2]
y3=K3*(10**x3)**Gamma[3]

p2=np.log10(K2*(10**rho2)**Gamma[2])

print(rho0)

plt.figure(0)
plt.grid(linewidth=0.2)
plt.plot(x0,np.log10(y0),label=r'$\Gamma_0=1.35692$ (SLy)')
plt.plot(x1,np.log10(y1),label=r'$\Gamma_1=3.446$')
plt.plot(x2,np.log10(y2),label=r'$\Gamma_2=3.572$')
plt.plot(x3,np.log10(y3),label=r'$\Gamma_3=2.887$')
plt.vlines(rho1,30,38,linestyles=(0, (5, 5)),linewidth=0.7,colors='black')
plt.vlines(rho2,30,38,linestyles=(0, (5, 5)),linewidth=0.7,colors='black')
plt.vlines(rho0,30,p0,linestyles=(0, (5, 5)),linewidth=0.7,colors='black')
plt.hlines(p1,13.9,rho1,linestyles=(0, (5, 5)),linewidth=0.7,colors='black')
plt.text(14.05, p1+0.1, r'$p_1$')
plt.text(14.65, 32.2, r'$\rho_1$')
plt.text(14.95, 32.2, r'$\rho_2$')
plt.text(rho0-0.07, 32.2, r'$\rho_{int}$')
plt.xlim([14.0,15.3])
plt.ylim([32.0,37.0])
plt.legend(loc='best')
plt.xlabel(r"log($\rho$ [gcm$^{-3}$])")
plt.ylabel(r"log($p$ [dyne/cm$^{2}$])")
plt.title('Piecewise polytropic EoS (MPA1)')
plt.savefig('Read_PWP.png',dpi=300,bbox_inches='tight')
plt.show()



