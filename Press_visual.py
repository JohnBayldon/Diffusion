import tkinter as tk
import scipy as sp
import matplotlib
matplotlib.use('TkAgg') # Change this as desired.
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2Tk
from matplotlib.backend_bases import key_press_handler
from pylab import plt,subplot,get_current_fig_manager,cm,show
import numpy as np




class MainApp(tk.Frame):
    def __init__(self,parent,*args,**kwargs):
        tk.Frame.__init__(self,parent,*args,**kwargs)
        self.parent=parent
        self.running=False
        self.time=0
        self.initWidgets()
        self.SetInitialConditions()
        self.im = self.img.imshow(self.u0[:,45,:], cmap=cm.hot, interpolation='nearest', origin='lower')
        self.fig.colorbar( self.im ) # Show the colorbar along the side
        self.TempPlot,=self.TempGraph.plot(self.u0[25,45,:])
        self.TempGraph.set_ylim([0,250])
        self.TempGraph.set_ylabel('Temperature, T, ($^\circ$ C)')
        self.TempGraph.set_xlabel('X index, j')
        show()






    def initWidgets(self):
        self.fig = plt.figure(1)
        self.img = subplot(111)
        self.manager=get_current_fig_manager()
        self.img = subplot(2,1,2)
        self.TempGraph=subplot(2,1,1)






        row=0
        self.grid()
        self.lblPower=tk.Label(self,text="Power")
        self.lblPower.grid(row=row,column=0)
        self.sclPower=tk.Scale(self,from_=0,to_=100,orient=tk.HORIZONTAL)
        self.sclPower.grid(row=row,column=1,columnspan=3)

        row=row+1
        self.lblTime=tk.Label(self,text="Time={0}".format(self.time))
        self.lblTime.grid(row=row,column=0)

        #lastrow
        row=row+1
        self.btnOne=tk.Button(master=self,text="Run")
        self.btnOne["command"]=self.Run
        self.btnOne.grid(row=row,column=0)
        self.btnTwo=tk.Button(master=self,text="Soak")
        self.btnTwo["command"]=self.Soak
        self.btnTwo.grid(row=row,column=2)

        self.QUIT=tk.Button(master=self,text="QUIT")
        self.QUIT["command"]=self.quit
        self.QUIT.grid(row=row,column=3)

    def Run(self):
        self.running=True
        self.btnOne["command"]=self.Stop
        self.btnOne["text"]="Stop"
        self.btnTwo["state"]=tk.DISABLED

        while self.running:
            self.q=self.sclPower.get()/1000.0
            self.updateLoop()
            self.update()

    def Soak(self):
        self.running=True
        self.btnTwo["command"]=self.Stop
        self.btnTwo["text"]="Stop"
        self.btnOne["state"]=tk.DISABLED
        self.q=0.0
        while self.running:
            self.updateLoop()


    def Stop(self):
        self.running=False
        self.btnOne["command"]=self.Run
        self.btnOne["text"]="Run"
        self.btnOne["state"]=tk.NORMAL
        self.btnTwo["command"]=self.Soak
        self.btnTwo["text"]="Soak"
        self.btnTwo["state"]=tk.NORMAL

    def quit(self):
        plt.close()
        root.quit()
        root.destroy()

    def SetInitialConditions(self):
        cr = 5
        # size (mm)
        w = 25.4 * 12
        l = 25.4 * 12
        h = np.array([25.4, 25.4, 10, 10 * cr])
        dz = 1.0

        finalVf = 0.2
        fill = 0.25
        rho = np.array([0.28, 7.8, 0.2, 1.0]) * 1e3  # denity kg/m3
        k = np.array([0.08, 20, 1, 20])  # Thermal condictivity
        shc = np.array([680, 470, 710 * finalVf + 1600 * (1 - finalVf) * fill,
                        710 * finalVf + 1600 * (1 - finalVf) * fill])  # J/(kgK)
        self.D = k / (rho * shc) * 1e6  # thermal conductivity/(density * specific heat capacity) (mm2)

        nx, ny = 50, 50  # element numbers in x dirextion
        nz = np.array([25, 25, 10, 25])  # elemnet mnumbers in Z dierction, insulator, steel, compositie note symmetric
        self.bz = np.array([nz[0], nz[0] + nz[1], nz[0] + nz[1] + nz[2]])

        dx = w / nx
        dy = l / ny
        self.dx2 = dx * dx
        self.dy2 = dy * dy
        self.dz2 = dz * dz
        print('stability dt(0.001)<{}'.format(1e-6 * rho * shc * (self.dx2 + self.dy2 + self.dz2) / (2 * k)))
        self.dt = .02
        print(self.dx2, self.dy2, self.dz2)
        # initial conditions
        T0 = 25  # temp in k
        # power = 5kW
        Power = 5000/4
        # heater elements


        SC = 1e8

        # define coefficient matrix
        self.a = np.zeros((nx, ny, self.bz[2]), dtype=np.float16)
        # diffusivities
        self.a[:, :, 0:nz[0]] = self.D[0]
        self.a[:, :, nz[0]:(nz[0] + nz[1])] = self.D[1]
        self.a[:, :, (nz[0] + nz[1]):(nz[0] + nz[1] + nz[2])] = self.D[2]
        self.D01 = 2 * self.D[0] * self.D[1] / (self.D[0] + self.D[1])
        self.D12 = 2 * self.D[1] * self.D[2] / (self.D[1] + self.D[2])

        self.u0 = T0 * np.ones((nx, ny, (nz[0] + nz[1] + nz[2])))
        self.u = np.empty_like(self.u0)
        self.heater = np.zeros_like(self.u0)
        self.heater[5:50, :, 30:35] = 1.0






        num_heater_vols = self.heater.sum()
        temp_rise_per_time_step = Power * self.dt / (num_heater_vols * (dx * dy * dz * rho[1] * shc[1]) / 1e9)
        print('Temp rise in heaters = {}'.format(temp_rise_per_time_step))
        self.heater *= temp_rise_per_time_step





    def updateLoop(self):
        for n in range(10):
            self.UpdateFIgure()
            self.time+=self.dt

            self.update()
        self.lblTime['text']="Time={0:f}s".format(self.time)
        self.im.set_array(self.u[:,45,:])
        self.im.set_clim([0,250])
        self.TempPlot.set_ydata(self.u[25,45,:])
           # =self.TempGraph.plot(self.T[self.newVals,20,:])
        self.manager.canvas.draw()


    def UpdateFIgure(self):
        # sort out the middle volumes diffusion in x and y directions
        self.u[1:-1, 1:-1, 1:] = (self.u0[1:-1, 1:-1, 1:] + self.a[1:-1, 1:-1, 1:] * self.dt * (
            (self.u0[2:, 1:-1, 1:] - 2 * self.u0[1:-1, 1:-1, 1:] + self.u0[:-2, 1:-1, 1:]) / self.dx2
            + (self.u0[1:-1, 2:, 1:] - 2 * self.u0[1:-1, 1:-1, 1:] + self.u0[1:-1, :-2, 1:]) / self.dy2))

        # Now for the z direction in the middle no change insulator
        self.u[:, :, 1:(self.bz[0] - 1)] = self.u0[:, :, 1:(self.bz[0] - 1)] + self.D[0] * self.dt * (
            (self.u0[:, :, 2:self.bz[0]] - 2 * self.u0[:, :, 1:(self.bz[0] - 1)] + self.u0[:, :, :(self.bz[0] - 2)]) / self.dz2
        )

        # boundary from steel to insulation
        self.u[:, :, self.bz[0] - 1] = self.u0[:, :, self.bz[0] - 1] + self.dt * (
            self.D01 * (self.u0[:, :, self.bz[0]] - self.u0[:, :, self.bz[0] - 1])
            - self.D[0] * (self.u0[:, :, self.bz[0] - 1] - self.u0[:, :, self.bz[0] - 2])
        )
        self.u[:, :, self.bz[0]] = self.u0[:, :, self.bz[0]] + self.dt * (
            +self.D[1] * (self.u0[:, :, self.bz[0] + 1] - self.u0[:, :, self.bz[0]])
            - self.D01 * (self.u0[:, :, self.bz[0]] - self.u0[:, :, self.bz[0] - 1])
        )

        # steel
        self.u[:, :, (self.bz[0] + 1):(self.bz[1] - 1)] = self.u[:, :, (self.bz[0] + 1):(self.bz[1] - 1)] + self.D[1] * self.dt * (
            (self.u0[:, :, (self.bz[0] + 2):(self.bz[1])]
             - 2 * self.u0[:, :, (self.bz[0] + 1):(self.bz[1] - 1)]
             + self.u0[:, :, (self.bz[0]):(self.bz[1] - 2)]) / self.dz2
        )
        # steel to composite Boundary
        self.u[:, :, self.bz[1] - 1] = self.u0[:, :, self.bz[1] - 1] + self.dt * (
            +self.D12 * (self.u0[:, :, self.bz[1]] - self.u0[:, :, self.bz[1] - 1]
                    - self.D[1] * (self.u0[:, :, self.bz[1] - 1] - self.u0[:, :, self.bz[1] - 2]))
        )
        self.u[:, :, self.bz[1]] = self.u0[:, :, self.bz[1]] + self.dt * (
            self.D[2] * (self.u0[:, :, self.bz[1] + 1] - self.u0[:, :, self.bz[1]])
            - self.D12 * (self.u0[:, :, self.bz[1]] - self.u0[:, :, self.bz[1] - 1])
        )

        # composite:
        self.u[1:-1, 1:-1, (self.bz[1] + 1):-1] = self.u[1:-1, 1:-1, (self.bz[1] + 1):-1] + self.D[2] * self.dt * (
            (self.u0[1:-1, 1:-1, (self.bz[1] + 2):]
             - 2 * self.u0[1:-1, 1:-1, (self.bz[1] + 1):- 1]
             + self.u0[1:-1, 1:-1, (self.bz[1]):-2]) / self.dz2
        )

        # central symmetry Note X and Y flows dealt with in main step
        self.u[:, :, -1] = self.u0[:, :, -1] + self.D[2] * self.dt * (
            (-self.u0[:, :, -1]
             + self.u0[:, :, -2]) / self.dz2)
        # x boundary y flows
        self.u[-1, 1:-1, :] = self.u0[-1, 1:-1, :] + self.a[-1, 1:-1, :] * self.dt * (
            (self.u0[-1, 2:, :] - 2 * self.u0[-1, 1:-1, :] + self.u0[-1, :-2, :]) / self.dy2
        )
        # x boundary x flows
        self.u[-1, :, :] = self.u0[-1, :, :] + self.a[-1, :, :] * self.dt * (
            (self.u0[-2, :, :] - self.u0[-1, :, :]) / self.dx2
        )
        # y boundary x flows
        self.u[1:-1, -1, :] = self.u0[1:-1, -1, :] + self.a[1:-1, -1, :] * self.dt * (
            (self.u0[2:, -1, :] - 2 * self.u0[1:-1, -1, :] + self.u0[:-2, -1, :]) / self.dx2
        )
        # y boundary y flows
        self.u[:, -1, :] = self.u0[:, -1, :] + self.a[:, -1, :] * self.dt * (
            (self.u0[:, -2, :] - self.u0[:, -1, :]) / self.dy2
        )

        # now for the edges
        # outsides TODO deal with convection losses
        # todo force bottom to constant temp?
        self.u[0, 1:-1, 1:-1] = (self.u0[0, 1:-1, 1:-1])
        self.u[1:-1, 0, 1:-1] = (self.u0[1:-1, 0, 1:-1])
        # u[1:-1,1:-1,0] = (self.u0[1:-1,1:-1,0]+self.dt*(
        #    (u0[2:,1:-1,0]-2*u0[:-2,1:-1,0]+u0[:-2,1:-1,0])/self.dx2
        #    +(u0[1:-1,2:,0]-2*u0[1:-1,1:-1,0]+u0[1:-1,:-2,0])/self.dy2
        #    +(u0[1:-1,1:-1,1]-2*u0[1:-1,1:-1,0])/self.dz2))

        # egdes
        self.u[0, 0, 1:-1] = self.u0[0, 0, 1:-1]
        self.u[0, 1:-1, 0] = self.u0[0, 1:-1, 0]
        self.u[1:-1, 0, 0] = self.u0[1:-1, 0, 0]
        self.u[0, 0, 0] = self.u0[0, 0, 0]

        # edges
        self.u[1:-1, -1, -1] = self.u0[1:-1, -1, -1] + self.a[1:-1, -1, -1] * self.dt * (
            (self.u0[2:, -1, -1] - 2 * self.u0[1:-1, -1, -1] + self.u0[:-2, -1, -1]) / self.dx2 +
            (self.u0[1:-1, -2, -1] - self.u0[1:-1, -1, -1]) / self.dy2
        )
        self.u[-1, 1:-1, -1] = self.u0[-1, 1:-1, -1] + self.a[-1, 1:-1, -1] * self.dt * (
            (self.u0[-2, 1:-1, -1] - self.u0[-1, 1:-1, -1]) / self.dx2 +
            (self.u0[-1, 2:, -1] - 2 * self.u0[1:-1, -1, -1] + self.u0[-1, :-2, -1]) / self.dy2
        )
        self.u[-1, -1, 1:-1] = self.u0[-1, -1, 1:-1] + self.a[-1, -1, 1:-1] * self.dt * (
            (self.u0[-2, -1, 1:-1] - self.u0[-1, -1, 1:-1]) / self.dx2 +
            (self.u0[-1, -2, 1:-1] - self.u0[-1, -1, 1:-1]) / self.dy2
        )
        # center
        self.u[-1, -1, -1] = self.u0[-1, -1, -1] + self.a[-1, -1, -1] * self.dt * (
            (self.u0[-2, -1, -1] - self.u0[-1, -1, -1]) / self.dx2 +
            (self.u0[-1, -2, -1] - self.u0[-1, -1, -1]) / self.dy2 +
            (self.u0[-1, -1, -2] - self.u0[-1, -1, -1]) / self.dz2
        )

        # combination edges
        self.u[0, -1, 1:-1] = self.u0[0, -1, 1:-1]
        self.u[-1, 0, 1:-1] = self.u0[-1, 0, 1:-1]
        self.u[0, 1:-1, -1] = self.u0[0, 1:-1, -1]
        self.u[-1, 1:-1, 0] = self.u0[-1, 1:-1, 0]
        self.u[1:-1, -1, 0] = self.u0[1:-1, -1, 0]
        self.u[1:-1, 0, -1] = self.u0[1:-1, 0, -1]
        self.u[-1, -1, 0] = self.u0[-1, -1, 0]
        self.u[-1, 0, -1] = self.u0[-1, 0, -1]
        self.u[0, -1, -1] = self.u0[0, -1, -1]
        self.u[-1, 0, 0] = self.u0[-1, 0, 0]
        self.u[0, -1, 0] = self.u0[0, -1, 0]
        self.u[0, 0, -1] = self.u0[0, 0, -1]

        # insulated boundary, constant temp
        self.u[:, :, 0] = 25

        # add some heat
        if self.u[1:-1, 1:-1, self.bz[0]:self.bz[1]].max() < 1000:
            self.u += self.heater
        # u[50:60,1:-1,30:35] += SC*dt/(rho[1]*shc[1])
        # [70:80,1:-1,30:35] += SC*dt/(rho[1]*shc[1])
        self.u0 = self.u.copy()



if __name__=="__main__":
    root=tk.Tk()
    root.wm_title("Basic Diffusion for resistance Heating")
    MainApp(root)
    root.mainloop()