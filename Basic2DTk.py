import tkinter as tk
import scipy as sp
import matplotlib
matplotlib.use('TkAgg') # Change this as desired.
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2TkAgg
from matplotlib.backend_bases import key_press_handler
from pylab import plt,subplot,get_current_fig_manager,cm,show





class MainApp(tk.Frame):
    def __init__(self,parent,*args,**kwargs):
        tk.Frame.__init__(self,parent,*args,**kwargs)
        self.parent=parent
        self.running=False
        self.time=0
        self.initWidgets()
        self.SetInitialConditions()
        self.im = self.img.imshow(self.T[self.newVals,:,:], cmap=cm.hot, interpolation='nearest', origin='lower')
        self.fig.colorbar( self.im ) # Show the colorbar along the side
        self.TempPlot,=self.TempGraph.plot(self.T[self.newVals,20,:])
        self.TempGraph.set_ylim([0,400])
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
        #basic stuff, evenntuallt allow this to be modofied
        self.dx=.5    #in mm default 50 microns
        self.dz=0.05    #in mm default 50 microns
        self.kx=[7.59/50,7.63/50]    #x direction diffusiom in W/mmsk [no poly, poly]
        self.kz=[0.19/1000,0.23/1000]   #zz direction diffusion in J/mmsk [no poly,poly]
        self.crho=[15.28/1000,195.8/1000]  #heat capacity [graphite,poly] J/(Kmm3)

        #self.a=0.5
        self.q=.1

        #CHECK UNITS
        self.nx=int(3.5*25.4/self.dx)
        self.nz=int(6/self.dz)
        self.nt=500

        #precalculate some things
        self.dx2=self.dx*self.dx
        self.dz2=self.dz*self.dz
        #maximum dt limited by stability criteria.
        self.dt=self.crho[0]*self.dx2*self.dz2/(2*(self.kx[0]*self.dz2+self.kz[0]*self.dx2))/2
        print(self.dt)
# Start u and ui off as zero matrices:
        self.T = sp.full([2,self.nz,self.nx],25)  #old and new values for T switch each step
        self.Coeffs=sp.full([5,self.nz,self.nx],1/self.crho[0])  #this array holds the coeffs for each volume initiall just graphite
        #polymer block is cube 20x20 in side
        self.Coeffs[0,6:116,50:130]=1/self.crho[1]

        self.newVals=0      #this holds the values for the current step
        self.u=sp.zeros([self.nx,self.nz])
        self.ui=sp.zeros([self.nx,self.nz])

        #initial conditions
        #for i in range(self.nx):
         #   for j in range(self.nz):
        #        if ( ( (i*self.dx-0.4)**2+(j*self.dz-0.5)**2 <= 0.05)):
        #        # & ((i*self.dx-0.5)**2+(j*self.dz-0.5)**2>=.05) ):
        #            self.T[j,i,:] = .000001

                #initial conditions
        #verify stsabiliy



    def updateLoop(self):
        for n in range(10):
            self.UpdateFIgure()
            self.time+=self.dt

            self.update()
        self.lblTime['text']="Time={0:f}s".format(self.time)
        self.im.set_array(self.T[self.newVals,:,:])
        self.im.set_clim([0,400])
        self.TempPlot.set_ydata(self.T[self.newVals,20,:])
           # =self.TempGraph.plot(self.T[self.newVals,20,:])
        self.manager.canvas.draw()


    def UpdateFIgure(self):
        old=1-self.newVals
        xCoeff=self.dt*self.kx[0]/(self.dx2)
        zCoeff=self.dt*self.kz[0]/(self.dz2)
        #first we update the bulk elements, for all of these we use full size elements (DxDz)
        #we assume in this loop we have calculated all the relevant coefficients.
        #starts from temp at last point plus source term Note Power Everywhere so this applies to all volumes
        self.T[self.newVals,:,:] =self.T[old,:,:]+self.dt*self.q*self.Coeffs[0,:,:]
        #diffusion from z neighbors
        self.T[self.newVals,1:-1, 1:-1] +=zCoeff*self.Coeffs[0,1:-1,1:-1]*(self.T[old,2:, 1:-1] - 2*self.T[old,1:-1, 1:-1] + self.T[old,:-2, 1:-1])
        #diffusion from x neighbors
        self.T[self.newVals,1:-1, 1:-1] +=xCoeff*self.Coeffs[0,1:-1,1:-1]*(self.T[old,1:-1, 2:] - 2*self.T[old,1:-1, 1:-1] + self.T[old,1:-1, :-2])
        #now the boundaries.
        #top and bottom Note, already delat with last time plus source



        #sides
        #up and donwn no different
        self.T[self.newVals,1:-1, 0] +=zCoeff*self.Coeffs[0,1:-1,0]*(self.T[old,2:, 0] - 2*self.T[old,1:-1, 0] + self.T[old,:-2, 0])
        self.T[self.newVals,1:-1, -1:] +=zCoeff*self.Coeffs[0,1:-1,-1:]*(self.T[old,2:, -1:] - 2*self.T[old,1:-1, -1:] + self.T[old,:-2, -1:])
        #convetion and radiation losses are proportional to area and temp difference
        t0=25  #extenral temp
        ConvLoss=.01
        self.T[self.newVals,1:-1, 0] +=xCoeff*self.Coeffs[0,1:-1,0]*(self.T[old,1:-1, 1] - self.T[old,1:-1, 0])
        self.T[self.newVals,:, 0] +=ConvLoss*self.dt*self.Coeffs[0,:,0]*(t0-self.T[old,:, 0])
        self.T[self.newVals,1:-1, -1:] +=xCoeff*self.Coeffs[0,1:-1,-1:]*(self.T[old,1:-1, -2:-1] - self.T[old,1:-1, -1:])
        self.T[self.newVals,:, -1:] +=ConvLoss*self.dt*self.Coeffs[0,:,-1:]*(t0-self.T[old,:, -1:])

        #Top anddd Bottom
        #sides are  no different
        self.T[self.newVals,0,1:-1] +=xCoeff*self.Coeffs[0,0,1:-1]*(self.T[old,0,2:] - 2*self.T[old,0,1:-1] + self.T[old,0,:-2])
        self.T[self.newVals,-1:,1:-1 ] +=xCoeff*self.Coeffs[0,-1:,1:-1]*(self.T[old,-1:,2:] - 2*self.T[old,-1:,1:-1] + self.T[old,-1:,:-2])
        #convetion and radiation losses are propostional to area and temp difference
        t0=25  #extenral temp
        ConvLoss=.1
        self.T[self.newVals,0,1:-1] +=zCoeff*self.Coeffs[0,0,1:-1]*(self.T[old,1,1:-1] - self.T[old,0,1:-1])
        self.T[self.newVals,0,:] +=ConvLoss*self.dt*self.Coeffs[0,0,:]*(t0-self.T[old,0,: ])
        self.T[self.newVals,-1:,1:-1] +=zCoeff*self.Coeffs[0,-1:,1:-1]*(self.T[old, -2:-1,1:-1] - self.T[old,-1:,1:-1])
        self.T[self.newVals,-1:,:] +=ConvLoss*self.dt*self.Coeffs[0,-1:,:]*(t0-self.T[old,-1:,:])

        #finally the  diffusion at the corners
        self.T[self.newVals,0, 0] +=zCoeff*self.Coeffs[0,0,0]*(self.T[old,1, 0] - self.T[old,0, 0] )
        self.T[self.newVals,0, 0] +=xCoeff*self.Coeffs[0,0,0]*(self.T[old,0, 1] - self.T[old,0, 0] )
        self.T[self.newVals,-1:, 0] +=zCoeff*self.Coeffs[0,0,0]*(self.T[old,-2:-1, 0] - self.T[old,-1:, 0] )
        self.T[self.newVals,-1:, 0] +=xCoeff*self.Coeffs[0,0,0]*(self.T[old,-1:, 1] - self.T[old,-1:, 0] )
        self.T[self.newVals,0, -1:] +=zCoeff*self.Coeffs[0,0,0]*(self.T[old,1, -1:] - self.T[old,0, -1:] )
        self.T[self.newVals,0, -1:] +=xCoeff*self.Coeffs[0,0,0]*(self.T[old,0, -2:-1] - self.T[old,0, -1:] )
        self.T[self.newVals,-1:, -1:] +=zCoeff*self.Coeffs[0,0,0]*(self.T[old,-2:-1, 0] - self.T[old,-1:, -1:] )
        self.T[self.newVals,-1:, -1:] +=xCoeff*self.Coeffs[0,0,0]*(self.T[old,0, -2:-1] - self.T[old,-1:, -1:] )


        self.newVals=old

if __name__=="__main__":
    root=tk.Tk()
    root.wm_title("Basic Diffusion for resistance Heating")
    MainApp(root)
    root.mainloop()