'''
Create a 1D explicit heat model fro fabric in press
model assumes heated steel plates at top and bottom
uses a symmetrical set up with nod N being symetic.
'''
import numpy as np

L = .05               #length of system in m
n = 10              #numnber if control vols
T_init = 25         # initial temperature in C
T_rate = 0        #temperature ramp rate in C/min
T_set   = 100       #setpoint temperature in C
t = 0               #start time
dt = 60              #time step
k = 800            #diffusion coeff carbon
rho = 0.05*1.75 * 10e3    #kg/m3
shc = 710            #J/kgK

x = np.array([i/n for i in range(n)],dtype = np.float32)  # co-ords of control points in model units
dx = np.diff(x)
T = np.ones_like(x)*T_init
#T[0] = 100
k_list = np.ones_like(x)*k  #references to points
k_cv =2* k_list[0:-1]*k_list[1:]/(k_list[0:-1]+k_list[1:])  #harmonic mean
a_list = k_cv/dx
a_po = np.ones_like(x)*rho*shc
a_po[0]= a_po[0]/2
a_po[-1] = a_po[-1]/2
s_c = np.ones_like(x)*0.1
s_p = np.zeros_like(x)
binder_content = np.ones_like(x)
dt_max = rho*shc*dx[0]**2/(2*k)
pass


def update(dt):
    global T,t, binder_content
    t += dt
    T_old = np.copy(T)
    #calculate new platen temp
    T[0] = T_old[0]+T_rate*dt/60
    #--------------------------source terms-----------------------------------------
    rate = np.minimum(np.maximum(T-220,0)*.001,1)
    s_c = rate*binder_content*400
    binder_content = binder_content - s_c/400




    #----------------------------------------------------------------------------------
    a_po_dt = a_po/dt

    T[1:-1]= (a_list[0:-1]*T[0:-2] + a_list[1:]*T[2:]+
              (a_po_dt[1:-1]-a_list[0:-1]-a_list[1:])*T_old[1:-1])/a_po_dt[1:-1] + s_c[1:-1]
    #symmetrical boundary
    T[-1] = (a_list[-2]*T[-2]+(a_po_dt[-1]-a_list[-2])*T_old[-1])/a_po_dt[-1]+s_c[-1]
    #populate tmeperature matrix

    pass
    #solve for temps




if __name__ == "__main__":
    T_init = 25  # initial temperature in C
    T_rate = 5  # temperature ramp rate in C/min
    T_set = 350  # setpoint temperature in C
    dt = dt_max*.9
    with open("outfile.txt",'w') as f:
        while(T[0]<T_set):
            update(dt=dt)

            t_str = ''.join(f"{int(i)}," for i in T)
            b_str = ''.join(f"{int(i*100)}," for i in binder_content)
            print(f"temperatures at time {t/60}min , {t_str} {b_str}")
            f.write(f"{t},{t_str}\n")
        T_rate = 0
        while(t<(200*60)):
            update(dt = dt)

            t_str = ''.join(f"{int(i)}," for i in T)
            b_str = ''.join(f"{int(i*100)}," for i in binder_content)
            print(f"temperatures at time {t/60}min , {t_str},{b_str}")
            f.write(f"{t},{t_str}\n")