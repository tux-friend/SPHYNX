import numpy as np
from scipy.integrate import ode, odeint, quad
from scipy.interpolate import interp1d
from scipy import optimize
from matplotlib import pyplot as plt

# constants [in cgs]
G = 6.674*1.e-8	
c = 2.99792458*1e+10
c2 = c*c	
M_sun = 1.989*1e+33
R = 1.3733e+6 #

#define PWP EoS
Gamma = np.array([1.35692e0, 3.446e0, 3.572e0, 2.887e0])
Gammashort = np.array([3.446e0, 3.572e0, 2.887e0])
p1 = 10**34.495
K0=3.59389e+13 #K0 of fixed crust EoS
K1=p1/((10**14.7)**3.446) #calculate K1 from p1
rhoth0 = (K1/K0)**(1/(Gamma[0]-Gamma[1])) #intersecting density of fixed crust EoS and piece with Gamma1 
rho_th=np.array([rhoth0, 10**14.7, 10**15.0]) #dividing densities 
rhoc=0.525819*1e+15
print(rhoth0)
K1=K0*rho_th[0]**(Gamma[0]-Gamma[1])
K2=K1*rho_th[1]**(Gamma[1]-Gamma[2])
K3=K2*rho_th[2]**(Gamma[2]-Gamma[3])
K=np.array([K0,K1,K2,K3])
Kshort=np.array([K1,K2,K3])
P_th=Kshort*rho_th**Gammashort

#rho_th =np.array([2.6278e+12, 10**14.7,10**15.0])
#K=np.array([3.99874e-07,10**34.495/(rho_th[1]**Gamma[1]),10**34.495/(rho_th[1]**Gamma[2]),10**34.495/(rho_th[1]**Gamma[2])*rho_th[2]**Gamma[2]/(rho_th[2]**Gamma[3])])
#P_th=np.array([K[0]*(rho_th[0])**Gamma[0],10**34.495,K[2]*rho_th[2]**Gamma[2]])
a = np.array([1.1428e-02, 1.2818e-02, -4.5072e-02])

#Internal energy
e = np.array([8.4377e-04, 1.8981e-03])

def eos(rho):
    if rho < rho_th[0]: 
        return K[0]*rho**Gamma[0]
    elif rho_th[0] <= rho < rho_th[1]:
        return K[1]*rho**Gamma[1]
    elif rho_th[1] <= rho < rho_th[2]: 
        return K[2]*rho**Gamma[2]
    else:
        return K[3]*rho**Gamma[3]
      
def inveos(P):
    if P < P_th[0]:
        return (P/K[0])**(1.0/Gamma[0])
    elif P_th[0] <= P < P_th[1]:
        return (P/K[1])**(1.0/Gamma[1])
    elif P_th[1] <= P < P_th[2]:
        return (P/K[2])**(1.0/Gamma[2])
    else:
        return (P/K[3])**(1.0/Gamma[3])

#Newtonian limit TOV --> equations of hydrostatic equilibrium
def tov(y,r):
    P,m=y[0],y[1]
    rho = inveos(P)
    dPdr = -G*m*rho/(r**2) 
    dmdr = 4.0*np.pi*r**2*rho
    return np.array([dPdr,dmdr])

def tovsolve(rhoc):
    r=np.linspace(1.0,R,25003)
    m=np.zeros_like(r)
    P=np.zeros_like(r)
    rho=np.zeros_like(r)
    m[0]=4.0/3.0*np.pi*r[0]**3*rhoc
    P[0]=eos(rhoc)
    y0=np.array([P[0],m[0]])
    y=odeint(tov,y0,r)
    P[:]=y[:,0]
    m[:]=y[:,1]
    for i in range(len(P)):
        rho[i]=inveos(P[i])
    return r,P,rho,m
    
r,P,rho,m=tovsolve(rhoc)
np.savetxt('1D_density_profile_MPA1.txt',np.c_[r,rho,m/M_sun])

plt.plot(r,rho,label='MPA1')
plt.title(r'1D density profile for MPA1')
plt.legend(loc='best')
plt.xlim([0, 1.4e6])
plt.ylim([0, 0.6e15])
plt.xlabel(r"Radius $r$ [cm]")
plt.ylabel(r"Density $\rho$ [gcm$^{-3}$]")
plt.grid()
plt.savefig('densprofile_MPA1.png',dpi=300)

