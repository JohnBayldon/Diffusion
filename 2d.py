
"""
A program which uses an explicit finite difference
scheme to solve the diffusion equation with fixed
boundary values and a given initial value for the
density.

Two steps of the solution are stored: the current
solution, u, and the previous step, ui. At each time
step, u is calculated from ui. u is moved to ui at the
end of each time-step to move forward in time.

"""
import scipy as sp
import time

import matplotlib
matplotlib.use('TkAgg') # Change this as desired.
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2TkAgg
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
from pylab import *
import Tkinter as Tk






#set up display stuff
root=Tk.Tk()
root.wm_title("Diffusion")
f=Figure(figsize=(5,4),dpi=100)
sub=f.add_subplot(111)

t = arange(0.0, 3.0, 0.01)
s = sin(2*pi*t)

sub.plot(t, s)



# a tk.DrawingArea
canvas = FigureCanvasTkAgg(f, master=root)
canvas.show()
canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

toolbar = NavigationToolbar2TkAgg(canvas, root)
toolbar.update()
canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

def on_key_event(event):
    updatefig()
    key_press_handler(event, canvas, toolbar)

canvas.mpl_connect('key_press_event', on_key_event)


def _quit():
    root.quit()     # stops mainloop
    root.destroy()  # this is necessary on Windows to prevent
                    # Fatal Python Error: PyEval_RestoreThread: NULL tstate

button = Tk.Button(master=root, text='Quit', command=_quit)
button.pack(side=Tk.BOTTOM)





# Declare some variables:
dx=0.01        # Interval size in x-direction.
dy=0.01        # Interval size in y-direction.
a=0.5          # Diffusion constant.
timesteps=500  # Number of time-steps to evolve system.
nx = int(1/dx)
ny = int(1/dy)


dx2=dx**2 # To save CPU cycles, we'll compute Delta x^2
dy2=dy**2 # and Delta y^2 only once and store them.
# For stability, this is the largest interval possible
# for the size of the time-step:
dt = dx2*dy2/( 2*a*(dx2+dy2) )
# Start u and ui off as zero matrices:
ui = sp.zeros([nx,ny])
u = sp.zeros([nx,ny])

# Now, set the initial conditions (ui).
for i in range(nx):
    for j in range(ny):
      if ( ( (i*dx-0.5)**2+(j*dy-0.5)**2 <= 0.1)
         & ((i*dx-0.5)**2+(j*dy-0.5)**2>=.05) ):
            ui[i,j] = 1

def evolve_ts(u, ui):
   global nx, ny
   u[1:-1, 1:-1] = ui[1:-1, 1:-1] + a*dt*( (ui[2:, 1:-1] - 2*ui[1:-1, 1:-1] + ui[:-2, 1:-1])/dx2 + (ui[1:-1, 2:] - 2*ui[1:-1, 1:-1] + ui[1:-1, :-2])/dy2 )



def updatefig(*args):
    global u, ui, m
    im.set_array(ui)
    manager.canvas.draw()
    # Uncomment the next two lines to save images as png
    # filename='diffusion_ts'+str(m)+'.png'
    # fig.savefig(filename)
    u[1:-1, 1:-1] = ui[1:-1, 1:-1] + a*dt*(
      (ui[2:, 1:-1] - 2*ui[1:-1, 1:-1] + ui[:-2, 1:-1])/dx2
      + (ui[1:-1, 2:] - 2*ui[1:-1, 1:-1] + ui[1:-1, :-2])/dy2 )
    ui = sp.copy(u)
    m+=1
    print "Computing and rendering u for m =", m
    if m >= timesteps:
        return False
    return True



fig = plt.figure(1)
img = subplot(111)
im = img.imshow( ui, cmap=cm.hot, interpolation='nearest', origin='lower')
manager = get_current_fig_manager()
m=1
fig.colorbar( im ) # Show the colorbar along the side

while updatefig():
    show()

Tk.mainloop()






