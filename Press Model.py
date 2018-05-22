import numpy as np
import matplotlib.pyplot as plt

#compression ration
cr = 5
#size (mm)
w = 25.4*12
l = 25.4*12
h = np.array([25.4,25.4,10,10*cr])
dz = 1.0


finalVf = 0.2
fill = 0.25
rho = np.array([0.28,7.8,1.0,0.2]) *10e3 #denity Mg/m3
k = np.array([0.08,20,20,1] )  #Thermal condictivity
shc = np.array([680, 470, 0.71*finalVf + 1.6*(1-finalVf)*fill ,0.71*finalVf + 1.6*(1-finalVf)*fill ]) # J/(kgK)
D = k/(rho * shc) *10e6  #thermal conductivity/(density * specific heat capacity)

nx , ny = 50 , 50 #element numbers in x dirextion
nz = np.array([25,25,10,25]) # elemnet mnumbers in Z dierction, insulator, steel, compositie note symmetric

dx = w/nx
dy = l/ny
dx2 = dx*dx
dy2 = dy*dy
dz2 = dz*dz
dt = 0.1

#initial conditions
T0 = 300 #temp in k
SC = 1e9
#define coefficient matrix
a = np.zeros((nx,ny,(nz[0]+nz[1]+nz[2])),dtype=np.float16)
#diffusivities
a[:,:,0:nz[0]]     = D[0]
a[:,:,nz[0]:(nz[0]+nz[1])] = D[1]
a[:,:,(nz[0]+nz[1]):(nz[0]+nz[1]+nz[2])] = D[2]


u0 = T0*np.ones((nx,ny,(nz[0]+nz[1]+nz[2])))
u =  np.empty_like(u0)

def do_timestep(u0,u):
    #sort out the middle volumes
    u[1:-1,1:-1,1:-1] = (u0[1:-1,1:-1,1:-1]+a[1:-1,1:-1,1:-1]*dt*(
        (u0[2:,1:-1,1:-1]-2*u0[1:-1,1:-1,1:-1]+u0[:-2,1:-1,1:-1])/dx2
        +(u0[1:-1,2:,1:-1]-2*u0[1:-1,1:-1,1:-1]+u0[1:-1,:-2,1:-1])/dy2
        +(u0[1:-1,1:-1,2:]-2*u0[1:-1,1:-1,1:-1]+u0[1:-1,1:-1,:-2])/dz2))
    #now for the edges
    #outsides TODO deal with convection losses
    #todo force bottom to constant temp?
    u[0,1:-1,1:-1] = (u0[0,1:-1,1:-1]+a[0,1:-1,1:-1]*dt*(
        (u0[1,1:-1,1:-1]-2*u0[0,1:-1,1:-1])/dx2
        +(u0[0,2:,1:-1]-2*u0[0,1:-1,1:-1]+u0[0,:-2,1:-1])/dy2
        +(u0[0,1:-1,2:]-2*u0[0,1:-1,1:-1]+u0[0,1:-1,:-2])/dz2))
    u[1:-1,0,1:-1] = (u0[1:-1,0,1:-1]+a[1:-1,0,1:-1]*dt*(
        (u0[2:,0,1:-1]-2*u0[1:-1,0,1:-1]+u0[:-2,0,1:-1])/dx2
        +(u0[1:-1,1,1:-1]-2*u0[1:-1,0,1:-1])/dy2
        +(u0[1:-1,0,2:]-2*u0[1:-1,0,1:-1]+u0[1:-1,0,:-2])/dz2))
    #u[1:-1,1:-1,0] = (u0[1:-1,1:-1,0]+dt*(
    #    (u0[2:,1:-1,0]-2*u0[:-2,1:-1,0]+u0[:-2,1:-1,0])/dx2
    #    +(u0[1:-1,2:,0]-2*u0[1:-1,1:-1,0]+u0[1:-1,:-2,0])/dy2
    #    +(u0[1:-1,1:-1,1]-2*u0[1:-1,1:-1,0])/dz2))

    #egdes
    u[0,0,1:-1] = u0[0,0,1:-1]+ a[0,0,1:-1]*dt*(
        (u0[1,0,1:-1]-2*u0[0,0,1:-1])/dx2+
        (u0[0,1,1:-1]-2*u0[0,0,1:-1])/dy2+
        (u0[0,0,2:]-2*u0[0,0,1:-1]+u0[0,0,:-2])/dz2
    )
    #u[0,1:-1,0] = u0[0,1:-1,0]+ a[0,1:-1,0]*dt*(
    #    (u0[1,1:-1,0]-2*u0[0,1:-1,0])/dx2+
   #     (u0[0,2:,0]-2*u0[0,1:-1,0]+u0[0,:-2,0])/dy2+
   #     (u0[0,1:-1,1]-2*u0[0,1:-1,0])/dz2
    #)
    #u[1:-1,0,0] = u0[1:-1,0,0]+ a[1:-1,0,0]*dt*(
    #    (u0[2:,1,0,]-2*u0[1:-1,0,0]+u0[:-2,0,0])/dx2+
    #    (u0[1:-1,1,0]-2*u0[1:-1,0,0])/dy2+
    #    (u0[1:-1,0,1]-2*u0[1:-1,0,0])/dz2
    #)
    #u[0,0,0] = u0[0,0,0]+a[0,0,0]*dt*(
    #    (u0[1,0,0]-2*u0[0,0,0])/dx2+
    #    (u0[0,1,0]-2*u0[0,0,0])/dy2+
    #    (u0[0,0,1]-2*u0[0,0,0])/dz2
    #)

    #symmetrics boundaries at middle
    u[-1,1:-1,1:-1] = u0[-1,1:-1,1:-1]+a[-1,1:-1,1:-1]*dt*(
        (u0[-2,1:-1,1:-1]-u0[-1,1:-1,1:-1])/dx2+
        (u0[-1,2:,1:-1]-2*u0[-1,1:-1,1:-1]+u0[-1,:-2,1:-1])/dy2+
        (u0[-1,1:-1,2:]-2*u0[-1,1:-1,1:-1]+u0[-1,1:-1,:-2])/dz2
    )
    u[1:-1,-1,1:-1] = u0[1:-1,-1,1:-1]+a[1:-1,-1,1:-1]*dt*(
        (u0[2:,-1,1:-1]-2*u0[1:-1,-1,1:-1]+u0[:-2,-1,1:-1])/dx2+
        (u0[1:-1,-2,1:-1]-u0[1:-1,-1,1:-1])/dy2+
        (u0[1:-1,-1,2:]-2*u0[1:-1,-1,1:-1]+u0[1:-1,-1,:-2])/dz2
    )
    u[1:-1,1:-1,-1] = u0[1:-1,1:-1,-1]+a[1:-1,1:-1,-1]*dt*(
        (u0[2:,1:-1,-1]-2*u0[1:-1,1:-1,-1]+u0[:-2,1:-1,-1])/dx2+
        (u0[1:-1,2:,1]-2*u0[1:-1,1:-1,-1]+u0[1:-1,:-2,-1])/dy2+
        (u0[1:-1,1:-1,-2]-u0[1:-1,1:-1,-1])/dz2
    )
    u[1:-1,-1,-1] = u0[1:-1,-1,-1]+a[1:-1,-1,-1]*dt*(
        (u0[2:,-1,-1]-2*u0[1:-1,-1,1]+u0[:-2,-1,-1])/dx2+
        (u0[1:-1,-2,-1]-u0[1:-1,-1,-1])/dy2+
        (u0[1:-1,-1,-2]-u0[1:-1,-1,-1])/dz2
    )
    u[-1,1:-1,-1] = u0[-1,1:-1,-1]+a[-1,1:-1,-1]*dt*(
        (u0[-2,1:-1,-1]-u0[-1,1:-1,1])/dx2+
        (u0[-1,2:,-1]-2*u0[1:-1,-1,-1]-u0[-1,:-2,-1])/dy2+
        (u0[-1,1:-1,-2]-u0[-1,1:-1,-1])/dz2
    )
    u[-1,-1,1:-1] = u0[-1,-1,1:-1]+a[-1,-1,1:-1]*dt*(
        (u0[-2,-1,1:-1]-u0[-1,-1,1:-1])/dx2+
        (u0[-1,-2,1:-1]-u0[-1,-1,1:-1])/dy2+
        (u0[-1,-1,2:]-2*u0[-1,-1,1:-1]+u0[-1,-1,:-2])/dz2
    )
    u[-1,-1,-1] = u0[-1,-1,-1]+a[-1,-1,-1]*dt*(
        (u0[-2,-1,-1]-u0[-1,-1,-1])/dx2+
        (u0[-1,-2,-1]-u0[-1,-1,-1])/dy2+
        (u0[-1,-1,-2]-u0[-1,-1,-1])/dz2
    )

    #insulated boundary, constant temp
    u[:,:,0]=300




    #add some heat
    u[10:20,:,30:35] += SC*dt/(rho[1]*shc[1])
    u[30:40,:,30:35] += SC*dt/(rho[1]*shc[1])
    u0 = u.copy()
    return u0,u



# Number of timesteps
nsteps = 101
# Output 4 figures at these timesteps
mfig = [0, 10, 50, 100]
fignum = 0
fig = plt.figure()
for m in range(nsteps):
    u0, u = do_timestep(u0, u)
    if m in mfig:
        fignum += 1
        print(m, fignum)
        ax = fig.add_subplot(220 + fignum)
        im = ax.imshow(u[:,25,:].copy(), cmap=plt.get_cmap('hot'), vmin=300,vmax=500)
        ax.set_axis_off()
        ax.set_title('{:.1f} ms'.format(m*dt*1000))
fig.subplots_adjust(right=0.85)
cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
cbar_ax.set_xlabel('$T$ / K', labelpad=20)
fig.colorbar(im, cax=cbar_ax)
plt.show()




