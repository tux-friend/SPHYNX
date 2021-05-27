import numpy as np
from scipy.integrate import ode, odeint, quad
from scipy.interpolate import interp1d
from scipy import optimize
from matplotlib import pyplot as plt

# constants [in cgs]
G = 6.674*1.e-8	
c = 2.99792458*1e+10
c2 = c*c	
M_sun = 1.989*1e33

#Input
M = 1.4*M_sun #target mass
R = 1.3*1e6 #target radius

n=0.5 #Polytropic index n (Gamma = 1 + 1/n)
rb = 2.753 #right boundary of xi, change accordingly to xi_1

rho_average = 3*M/(4*np.pi*R**3) #average density


#Solve Lane-Emden equation with ODEint routine
def LaneEmden(y, xi, n):
    theta, phi = y 
    dydxi = [-phi/(xi**2),(xi**2)*np.sign(theta)*abs(theta)**n]
    return dydxi

y0=[1.,0.] #boundary conditions

xi = np.linspace(1e-10, rb, 25003)

sol = odeint(LaneEmden,y0,xi,(n,))

#Interpolate the data with cubic spline to find xi_1 and (xi_1)^2*|theta'(xi_1)|
f = interp1d(xi,sol[:,0], kind='cubic')
xi_1 = optimize.root_scalar(f, bracket=[1.0,rb],method='brentq')
g = interp1d(xi,sol[:,1],kind='cubic')

#Plot solution of Lane-Emden equation
plt.figure(0)
plt.plot(xi,sol[:,0])
plt.show()
plt.grid()
print(xi_1) 

#Scaling -- resubstitute xi and theta to R and rho
rho_c=rho_average*xi_1.root**3/(3*g(xi_1.root))
a=R/xi_1.root
K= 4*np.pi*G*a**2/((n+1)*rho_c**(1/n-1))
rho = rho_c*abs(sol[:,0])**n
r = xi[:]*a
P = K*rho**(1+1/n)

#calculate culmulative mass [in M_sun]
cul_m=[4./3.*np.pi*(r[1]-r[0])**3*rho_c/M_sun]
for i in range(1,len(xi)):
   dm = quad(lambda x: 4*np.pi*a**3*rho_c*f(x)**n*x**2,0,xi[i])
   cul_m=np.append(cul_m,dm[0]/M_sun)
   i = i+1

#export data as text file
data = np.c_[r, rho, cul_m]
polytrope = [K,n]
np.savetxt('1D_density_profile_'+str(n)+'.txt', data)
file = open('polytrope_parameters.txt', 'a')
file.write('n: '+str(n)+'\n') 
file.write('K: '+str(K)+'\n') 
file.write('R [cm]: '+str(R)+'\n') 
file.write('M [M_sun]: '+str(M/M_sun)+'\n') 
file.write('rho_c [gcm^-3]: '+str(rho_c)+'\n')

