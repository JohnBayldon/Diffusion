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
rho = np.array([0.28,7.8,0.2,1.0]) *1e3 #denity kg/m3
k = np.array([0.08,20,1,20] )  #Thermal condictivity
shc = np.array([680, 470, 710*finalVf + 1600*(1-finalVf)*fill ,710*finalVf + 1600*(1-finalVf)*fill ]) # J/(kgK)
D = k/(rho * shc) *1e6  #thermal conductivity/(density * specific heat capacity) (mm2)

nx , ny = 50 , 50 #element numbers in x dirextion
nz = np.array([25,25,10,25]) # elemnet mnumbers in Z dierction, insulator, steel, compositie note symmetric
bz = np.array([nz[0],nz[0]+nz[1],nz[0]+nz[1]+nz[2]])

dx = w/nx
dy = l/ny
dx2 = dx*dx
dy2 = dy*dy
dz2 = dz*dz
print('stability dt(0.001)<{}'.format(1e-6*rho*shc*(dx2+dy2+dz2)/(2*k)))
dt = .02
print(dx2,dy2,dz2)
#initial conditions
T0 = 300 #temp in k
#power = 5kW
Power = 5000
#heater elements


SC = 1e8

#define coefficient matrix
a = np.zeros((nx,ny,bz[2]),dtype=np.float16)
#diffusivities
a[:,:,0:nz[0]]     = D[0]
a[:,:,nz[0]:(nz[0]+nz[1])] = D[1]
a[:,:,(nz[0]+nz[1]):(nz[0]+nz[1]+nz[2])] = D[2]
D01 = 2*D[0]*D[1]/(D[0]+D[1])
D12 = 2*D[1]*D[2]/(D[1]+D[2])

u0 = T0*np.ones((nx,ny,(nz[0]+nz[1]+nz[2])))
u =  np.empty_like(u0)
heater = np.zeros_like(u0)
heater[10:20,:,30:40]=1.0
heater[30:40,:,30:40]=1.0
num_heater_vols = heater.sum()
temp_rise_per_time_step = Power*dt/(num_heater_vols*(dx*dy*dz*rho[1]*shc[1])/1e9)
print('Temp rise in heaters = {}'.format(temp_rise_per_time_step))
heater *= temp_rise_per_time_step
total_energy = heater.sum()


def do_timestep(u0,u):
    #sort out the middle volumes diffusion in x and y directions
    u[1:-1,1:-1,1:] = (u0[1:-1,1:-1,1:]+a[1:-1,1:-1,1:]*dt*(
        (u0[2:,1:-1,1:]-2*u0[1:-1,1:-1,1:]+u0[:-2,1:-1,1:])/dx2
        +(u0[1:-1,2:,1:]-2*u0[1:-1,1:-1,1:]+u0[1:-1,:-2,1:])/dy2))

    #Now for the z direction in the middle no change insulator
    u[:,:,1:(bz[0]-1)] = u[:,:,1:(bz[0]-1)]+D[0]*dt*(
        (u0[:,:,2:bz[0]]-2*u0[:,:,1:(bz[0]-1)]+u0[:,:,:(bz[0]-2)])/dz2
    )

    #boundary from steel to insulation
    u[:,:,bz[0]-1] = u0[:,:,bz[0]-1]+dt*(
        D01 * (u0[:, :, bz[0]] - u0[:, :, bz[0] - 1])
        -D[0] * (u0[:, :, bz[0] - 1]-u0[:, :, bz[0] - 2])
    )
    u[:, :, bz[0]] = u0[:, :, bz[0]] + dt * (
        +D[1] * (u0[:, :, bz[0] + 1] - u0[:, :, bz[0]])
        -D01 * (u0[:, :, bz[0]] -u0[:, :, bz[0] - 1])
    )

    #steel
    u[:, :, (bz[0]+1):(bz[1] - 1)] = u[:, :, (bz[0]+1):(bz[1]-1)] + D[1] * dt * (
        (u0[:, :, (bz[0]+2):(bz[1])]
         - 2 * u0[:, :, (bz[0]+1):(bz[1] - 1)]
         + u0[:, :, (bz[0]):(bz[1] - 2)]) / dz2
    )
    #steel to composite Boundary
    u[:,:,bz[1]-1]= u0[:,:,bz[1]-1]+dt*(
        +D12 * (u0[:, :, bz[1]] - u0[:, :, bz[1] - 1]
        -D[1]*(u0[:,:,bz[1]-1]-u0[:,:,bz[1]-2]))
    )
    u[:,:,bz[1]]= u0[:,:,bz[1]]+dt*(
        D[2]*(u0[:,:,bz[1]+1]-u0[:,:,bz[1]])
        -D12*(u0[:,:,bz[1]]-u0[:,:,bz[1]-1])
    )



    #composite:
    u[1:-1, 1:-1, (bz[1]+1):-1] = u[1:-1, 1:-1, (bz[1]+1):-1] + D[2] * dt * (
        (u0[1:-1, 1:-1, (bz[1]+2):]
         - 2 * u0[1:-1, 1:-1, (bz[1]+1):- 1]
         + u0[1:-1, 1:-1, (bz[1]):-2]) / dz2
    )

    #central symmetry Note X and Y flows dealt with in main step
    u[:,:,-1]= u0[:, :, -1] + D[2] * dt * (
        (-u0[:, :, -1]
         + u0[:, :, -2]) / dz2)
    #x boundary y flows
    u[-1,1:-1,:] = u0[-1,1:-1,:]+a[-1,1:-1,:]*dt*(
        (u0[-1,2:,:]-2*u0[-1,1:-1,:]+u0[-1,:-2,:])/dy2
    )
    #x boundary x flows
    u[-1,:,:] = u0[-1,:,:]+a[-1,:,:]*dt*(
        (u0[-2,:,:]-u0[-1,:,:])/dx2
    )
    #y boundary x flows
    u[1:-1,-1,:] = u0[1:-1,-1,:]+a[1:-1,-1,:]*dt*(
        (u0[2:,-1,:]-2*u0[1:-1,-1,:]+u0[:-2,-1,:])/dx2
    )
    #y boundary y flows
    u[:,-1,:] = u0[:,-1,:]+a[:,-1,:]*dt*(
        (u0[:,-2,:]-u0[:,-1,:])/dy2
    )


    #now for the edges
    #outsides TODO deal with convection losses
    #todo force bottom to constant temp?
    u[0,1:-1,1:-1] = (u0[0,1:-1,1:-1])
    u[1:-1,0,1:-1] = (u0[1:-1,0,1:-1])
    #u[1:-1,1:-1,0] = (u0[1:-1,1:-1,0]+dt*(
    #    (u0[2:,1:-1,0]-2*u0[:-2,1:-1,0]+u0[:-2,1:-1,0])/dx2
    #    +(u0[1:-1,2:,0]-2*u0[1:-1,1:-1,0]+u0[1:-1,:-2,0])/dy2
    #    +(u0[1:-1,1:-1,1]-2*u0[1:-1,1:-1,0])/dz2))

    #egdes
    u[0,0,1:-1] = u0[0,0,1:-1]
    u[0,1:-1,0] = u0[0,1:-1,0]
    u[1:-1,0,0] = u0[1:-1,0,0]
    u[0,0,0] = u0[0,0,0]


    #edges
    u[1:-1,-1,-1] = u0[1:-1,-1,-1]+a[1:-1,-1,-1]*dt*(
        (u0[2:,-1,-1]-2*u0[1:-1,-1,-1]+u0[:-2,-1,-1])/dx2+
        (u0[1:-1,-2,-1]-u0[1:-1,-1,-1])/dy2
    )
    u[-1,1:-1,-1] = u0[-1,1:-1,-1]+a[-1,1:-1,-1]*dt*(
        (u0[-2,1:-1,-1]-u0[-1,1:-1,-1])/dx2+
        (u0[-1,2:,-1]-2*u0[1:-1,-1,-1]+u0[-1,:-2,-1])/dy2
    )
    u[-1,-1,1:-1] = u0[-1,-1,1:-1]+a[-1,-1,1:-1]*dt*(
        (u0[-2,-1,1:-1]-u0[-1,-1,1:-1])/dx2+
        (u0[-1,-2,1:-1]-u0[-1,-1,1:-1])/dy2
    )
    #center
    u[-1,-1,-1] = u0[-1,-1,-1]+a[-1,-1,-1]*dt*(
        (u0[-2,-1,-1]-u0[-1,-1,-1])/dx2+
        (u0[-1,-2,-1]-u0[-1,-1,-1])/dy2+
        (u0[-1,-1,-2]-u0[-1,-1,-1])/dz2
    )

    #combination edges
    u[0,-1,1:-1]=u0[0,-1,1:-1]
    u[-1,0,1:-1]=u0[-1,0,1:-1]
    u[0,1:-1,-1]=u0[0,1:-1,-1]
    u[-1,1:-1,0]=u0[-1,1:-1,0]
    u[1:-1,-1,0]=u0[1:-1,-1,0]
    u[1:-1,0,-1]=u0[1:-1,0,-1]
    u[-1,-1,0]=u0[-1,-1,0]
    u[-1,0,-1]=u0[-1,0,-1]
    u[0,-1,-1]=u0[0,-1,-1]
    u[-1,0,0]=u0[-1,0,0]
    u[0,-1,0]=u0[0,-1,0]
    u[0,0,-1]=u0[0,0,-1]



    #insulated boundary, constant temp
    u[:,:,0]=300



    if u.min()<300:
        print('rats')
    #add some heat
    if u[1:-1,1:-1,bz[0]:bz[1]].max()<1000:
        print(u[1:-1,1:-1,bz[0]:bz[1]].max())
        u += heater
    #u[50:60,1:-1,30:35] += SC*dt/(rho[1]*shc[1])
    #[70:80,1:-1,30:35] += SC*dt/(rho[1]*shc[1])
    u0 = u.copy()
    return u0,u



# Number of timesteps
nsteps = 50001
# Output 4 figures at these timesteps
mfig = [10000, 20000, 40000, 50000]
fignum = 0
fig = plt.figure()
for m in range(nsteps):
    print(m)
    u0, u = do_timestep(u0, u)
    if m in mfig:
        fignum += 1
        print(m, fignum)
        ax = fig.add_subplot(220 + fignum)
        im = ax.imshow(u[:,45,:].copy(), cmap=plt.get_cmap('hot'), vmin=300,vmax=550)
        ax.set_axis_off()
        ax.set_title('{:.1f} min'.format(m*dt/60))
fig.subplots_adjust(right=0.85)
cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
cbar_ax.set_xlabel('$T$ / K', labelpad=20)
fig.colorbar(im, cax=cbar_ax)
plt.show()




