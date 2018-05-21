import numpy as np
import matplotlib.pyplot as plt

#compression ration
cr = 5
#size (mm)
w = 25.4*12
l = 25.4*12
h = np.array([25.4,25.4,10,10*cr])

finalVf = 0.2
fill = 0.25
rho = np.array([0.28,7.8,1.0,0.2]) *10e3 #denity Mg/m3
k = np.array([0.08,20,20,1] )  #Thermal condictivity
shc = np.array([680, 470, 0.71*finalVf + 1.6*(1-finalVf)*fill ,0.71*finalVf + 1.6*(1-finalVf)*fill ]) # J/(kgK)
D = k/(rho * shc) *10e6  #thermal conductivity/(density * specific heat capacity)

nx , ny = 50 , 50 #element numbers in x dirextion
nz = np.array([10,10,5,5]) # elemnet mnumbers in Z dierction, insulator, steel, compositie note symmetric

dx = w/nx
dy = l/ny
dz = h/nz

dt = 0.1

#initial conditions
T0 = 300 #temp in k

u0 = T0*np.ones(nx,ny,(nz[0]+nz[1]+nz[2]))


def do_timestep(u0,u):
    #use numpy stuff to make this work
    u[1:-1,1:-1,1:-1] = u0[1:-1,1:-1,1:-1]+dt*(
        (u0[2:,1:-1]-2*u0[1:-1,1:-1]+u0[:-2,1:-1])/dx2
        +(u0[2:,1:-1]-2*u0[1:-1,1:-1]+u0[:-2,1:-1])/dy2
        +(u0[2:,1:-1]-2*u0[1:-1,1:-1]+u0[:-2,1:-1])/dz2
    )
    u0 = u.copy()
    return u0,u
