__author__ = 'john_000'
import numpy as np
import matplotlib.pyplot as plt
import time,sys


nx=41
dx=2./(nx-1)
nt=20
nu=0.3
sigma=.2
dt=sigma*dx**2/nu



u=np.ones(nx)
u[.5/dx : 1/dx+1]=2
print u
un=[]

for n in range(nt):
    un.append(u)
    for i in range(1,nx-1):
        u[i]=un[n][i]+nu*dt/dx**2*(un[n][i+1]-2*un[n][i]+un[n][i-1])
    print u
    plt.plot(np.linspace(0,2,nx),u)
    plt.show()
        #raw_input("hit key for next")

nx=10
ny=10
nz=10
nt=10

u=np.ones((nt,nx,ny,nz))
#calculate K's
k=np.ones((nx,ny,nz))
stuff=0



#this is the important bit, here is where we timestep the whole thing calculating the new values from
#the previous timestep.
t=1
for x in range(nx):
    for y in range(ny):
        for z in range(nz):
            u[t+1,x,y,z]=(u[t,x-1,y,z]*k[x-1,y,z]+
                          u[t,x+1,y,z]*k[x+1,y,z]+
                          u[t,x,y-1,z]*k[x,y-1,z]+
                          u[t,x,y+1,z]*k[x,y+1,z]+
                          u[t,x,y,z-1]*k[x,y,z-1]+
                          u[t,x,y,z+1]*k[x,y,z+1]+
                          u(t,x,y,z)(stuff))




