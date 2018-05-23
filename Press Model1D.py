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
rho = np.array([0.28,7.8,1.0,1.0]) *10e3 #denity Mg/m3
k = np.array([0.08,20,20,20] )  #Thermal condictivity
shc = np.array([680, 470, 710*finalVf + 1600*(1-finalVf)*fill,710*finalVf + 1600*(1-finalVf)*fill ]) # J/(kgK)
D = k/(rho * shc) *1e6  #thermal conductivity/(density * specific heat capacity)
print('diffusivities {}'.format(D))
nx , ny = 1 , 1 #element numbers in x dirextion
nz = np.array([25,25,10,25]) # elemnet mnumbers in Z dierction, insulator, steel, compositie note symmetric
bz = np.array([nz[0],nz[0]+nz[1],nz[0]+nz[1]+nz[2]])

dx = w/nx
dy = l/ny
dx2 = dx*dx
dy2 = dy*dy
dz2 = dz*dz
print('stability dt(0.001)<{}'.format(1e-6*rho*shc*(dx2+dy2+dz2)/(2*k)))
dt = .1
print(dx2,dy2,dz2)
#initial conditions
T0 = 300 #temp in k
SC = 1e9
#define coefficient matrix
a = np.zeros((nx,ny,bz[2]),dtype=np.float16)
#diffusivities
a[:,:,0:nz[0]]     = D[0]
a[:,:,nz[0]:(nz[0]+nz[1])] = D[1]
a[:,:,(nz[0]+nz[1]):(nz[0]+nz[1]+nz[2])] = D[2]
D01 = 2*D[0]*D[1]/(D[0]+D[1])
D12 = 2*D[1]*D[2]/(D[1]+D[2])
print(D01,D12)
u0 = T0*np.ones((nx,ny,(nz[0]+nz[1]+nz[2])))
u =  np.empty_like(u0)
u0[0,0,45:48] = 350
print(u0[0,0,45:55])
def do_timestep(u0,u):
    #sort out the middle volumes diffusion in x and y directions


    #Now for the z direction in the middle no change insulator
    u[0,0,1:(bz[0]-1)] = u0[0,0,1:(bz[0]-1)]+D[0]*dt*(
        (u0[0,0,2:bz[0]]-2*u0[0,0,1:(bz[0]-1)]+u0[0,0,:(bz[0]-2)])/dz2
    )
    u[0, 0, bz[0] - 1] = u0[0, 0, bz[0] - 1] + dt * (
        D01 * (u0[0, 0, bz[0]] - u0[0, 0, bz[0] - 1]
        -D[0] * (u0[0, 0, bz[0] - 1]-u0[0, 0, bz[0] - 2]))
    )
    u[0, 0, bz[0]] = u0[0, 0, bz[0]] + dt * (
        +D[1] * (u0[0, 0, bz[0] + 1] - u0[0, 0, bz[0]])
        -D01 * (u0[0, 0, bz[0]] -u0[0, 0, bz[0] - 1])
    )

    #steel
    u[0, 0, (bz[0]+1):(bz[1] - 1)] = u0[0, 0, (bz[0]+1):(bz[1]-1)] + D[1] * dt * (
        (u0[0, 0, (bz[0]+2):(bz[1])]
         - 2 * u0[0, 0, (bz[0]+1):(bz[1] - 1)]
         + u0[0, 0, (bz[0]):(bz[1] - 2)]) / dz2
    )
    #composite Boundary
    u[0,0,bz[1]-1]= u0[0,0,bz[1]-1]+dt*(
        +D12 * (u0[0, 0, bz[1]] - u0[0, 0, bz[1] - 1]
        -D[1]*(u0[0,0,bz[1]-1]-u0[0,0,bz[1]-2]))
    )
    u[0,0,bz[1]]= u0[0,0,bz[1]]+dt*(
        D[2]*(u0[0,0,bz[1]+1]-u0[0,0,bz[1]])
        -D12*(u0[0,0,bz[1]]-u0[0,0,bz[1]-1])
    )

    #composite:
    u[0, 0, (bz[1]+1):-1] = u0[0, 0, (bz[1]+1):-1] + D[2] * dt * (
        (u0[0, 0, (bz[1]+2):]
         - 2 * u0[0, 0, (bz[1]+1):- 1]
         + u0[0, 0, (bz[1]):-2]) / dz2
    )

    #bound = u0[15, 25, (bz[1] - 5):(bz[1] + 5)]

    #other boundary

    #locked boundary
    u[0,0,bz[2]-1]= u0[0, 0, bz[2]-1] + D[2] * dt * (
        (-u0[0, 0, bz[2]-1]
         + u0[0, 0, bz[2]-2]) / dz2)





    #insulated boundary, constant temp
    u[0,0,0]=300



    if u.min()<300:
        print('rats')
    #add some heat
    if u[0,0,40:50].max()<400:
        print('heat on')
        u[0,0,30:45] += SC*dt/(rho[1]*shc[1])
    u0 = u.copy()
    return u0,u



# Number of timesteps
nsteps = 5001
# Output 4 figures at these timesteps
mfig = [300, 500, 800, 1000,2000,3000,4000,5000]
fignum = 0
fig = plt.figure()
sub = fig.add_subplot()
for m in range(nsteps):
    print(m)
    u0, u = do_timestep(u0, u)
    print(u0[0,0,45:60])
    energies = [sum(u0[0,0,:bz[0]])*D[0],sum(u0[0,0,bz[0]:bz[1]])*D[1],sum(u0[0,0,bz[1]:])*D[1]]
    print(energies, sum(energies))
    if m in mfig:
        plt.plot(u0[0,0,:])
        plt.show()




