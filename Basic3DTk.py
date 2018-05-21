import Tkinter as tk
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

        self.initWidgets()
        self.SetInitialConditions()
        self.im = self.img.imshow(self.T[self.newVals,:,:], cmap=cm.hot, interpolation='nearest', origin='lower')
        self.fig.colorbar( self.im ) # Show the colorbar along the side
        show()






    def initWidgets(self):
        self.fig = plt.figure(1)
        self.manager=get_current_fig_manager()
        self.img = subplot(2,1,1)
        self.TempGraph=subplot(2,1,2)
        x1=sp.linspace(0.0,5.0)
        y1=sp.cos(2*sp.pi*x1)*sp.exp(-x1)
        plt.plot(x1,y1)





        row=0
        self.grid()
        self.lblPower=tk.Label(self,text="Power")
        self.lblPower.grid(row=row,column=0)
        self.sclPower=tk.Scale(self,from_=0,to_=100000,orient=tk.HORIZONTAL)
        self.sclPower.grid(row=row,column=1,columnspan=3)




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
            self.q=self.sclPower.get()
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
            self.update()

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
        self.T0=25
        self.dx=0.05    #in mm default 50 microns
        self.zy=0.05
        self.dz=0.05    #in mm default 50 microns
        self.Lx=3.5*25.4
        self.Ly=3.5*25.4
        self.Lz=120*self.dx
        self.nx=int(self.Lx/self.dx)
        self.ny=int(self.Ly/self.dy)
        self.nz=int(self.Lz/self.dz)

        #create the array for temps, note old and new values alternate each step initializ at 25C
        self.T = sp.full([2,self.nz,self.ny,self.nx],self.T0)

        #set material matrix


        #use this to create the basic property matrices
        self.Material=sp.zeros_like(self.T[0,:,:,:],dtype=sp.bool_)
        self.Coeffs=sp.zeros([5,self.nz,self.ny,self.nx])
        self.k=sp.zeros([3,self.nz,self.ny,self.nx])
        #UpdateCoeffs


        #initiallized coefficient matrices




        #expelicit, we will generaly use the fastest time available initially
        self.dt=1  #seconds

        #material Properties for PA12 SHC=1.2     density=1.01 mg/mm3 Volume fraction ~15.2%
        # for Carbon SHC=072  denisty=1.79  volume fraction ~3.8%




        self.kx=[7.6,7.65]*1000     #x direction diffusiom in J/(mmks) [no poly, poly]
        self.kz=[.2,.3]*1000    #zz direction diffusion in J/(mmks) [no poly,poly]
        self.crho=[015,195]  #heat capacity [graphite,poly] J/(Kmm3)


        self.UpdateCoefficients()







        #self.a=0.5
        self.q=0

        #CHECK UNITS

        self.nt=500

        #precalculate some things
        self.dx2=self.dx*self.dx
        self.dz2=self.dz*self.dz
        self.dt=self.crho[0]*self.dx2*self.dz2/(2*(self.kx[0]*self.dz2+self.kz[0]*self.dx2))/2
# Start u and ui off as zero matrices:
        self.T = sp.full([2,self.nz,self.nx],25)  #old and new values for T switch each step
        self.Coeffs=sp.full([5,self.nz,self.nx],1/self.crho[0])  #this array holds the coeffs for each volume initiall just graphite
        #polymer block is cube 20x20 in side
        self.Coeffs[0,6:116,100:800]=1/self.crho[1]

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




    def UpdateCoefficients(self):
        #this subroutine updates the coefficients uses in the calculations
        #first the heat capacities for each cell (in J/ks)
        self.Coeffs[0,:,:,:]=1/(self.crho[0]*(self.dx*self.dy*self.dz))  #set all to be same as grapite
        self.Coeffs[0,100:200,100:200,10:110]=1/(self.chrho[1]*(self.dx*self.dy*self.dz)) #set a block to be value of graphite plus polymer
        #next the diffusivities in each direction
        self.k[2,:,:,:]=self.kx[0]
        self.k[1,:,:,:]=self.kx[0]
        self.k[0,:,:,:]=self.kz[0]
        self.k[2,100:200,100:200,10:110]=self.kx[1]
        self.k[1,100:200,100:200,10:110]=self.kx[1]
        self.k[0,100:200,100:200,10:110]=self.kz[1]
        #now some heavy lifting calculate harmonic means NOTE, last volume is redundant.
        self.Coeffs[1,:-1,:,:]=(2/self.dz*self.k[0,:-1,:,:]*self.k[0,1:,:,:]/(self.k[0,:-1,:,:]+self.k[0,1:,:,:]))
        self.Coeffs[2,:,-1,:]=(2/self.dy*self.k[1,:,:-1,:]*self.k[0,:,1:,:,]/(self.k[0,:,-1,:]+self.k[0,:,1:,:]))
        self.Coeffs[3,:,:,-1]=(2/self.dx*self.k[0,:,:,-1,]*self.k[0,:,:,1:]/(self.k[0,:,:,-1]+self.k[0,:,:,1:]))




    def updateLoop(self):
        #for n in range(100):

        self.UpdateFIgure()

        self.im.set_array(self.T[self.newVals,:,:])
        self.im.set_clim([0,400])

        self.manager.canvas.draw()


    def UpdateFIgure(self):
        old=1-self.newVals
        new=self.newVals
        #first add in last temp plus effect of source term
        self.T[new,:,:,:]=self.T[old,:,:,:]+self.q*self.Coeffs[0,:,:,:]














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
        #convetion and radiation losses are propostional to area and temp difference
        t0=25  #extenral temp
        ConvLoss=5000
        self.T[self.newVals,1:-1, 0] +=xCoeff*self.Coeffs[0,1:-1,0]*(self.T[old,1:-1, 1] - self.T[old,1:-1, 0])
        self.T[self.newVals,:, 0] +=-ConvLoss*self.dt*self.Coeffs[0,:,0]*(self.T[old,:, 0] -t0)
        self.T[self.newVals,1:-1, -1:] +=xCoeff*self.Coeffs[0,1:-1,-1:]*(self.T[old,1:-1, -2:-1] - self.T[old,1:-1, -1:])
        self.T[self.newVals,:, -1:] +=-ConvLoss*self.dt*self.Coeffs[0,:,-1:]*(self.T[old,:, -1:] -t0)

        #Top anddd Bottom
        #up and donwn no different
        self.T[self.newVals,0,1:-1] +=zCoeff*self.Coeffs[0,0,1:-1]*(self.T[old,0,2:] - 2*self.T[old,1,1:-1] + self.T[old,0,:-2])
        self.T[self.newVals,-1:,1:-1 ] +=zCoeff*self.Coeffs[0,-1:,1:-1]*(self.T[old,-1:,2:] - 2*self.T[old,-1:,1:-1] + self.T[old,-2:-1,1:-1])
        #convetion and radiation losses are propostional to area and temp difference
        t0=25  #extenral temp
        ConvLoss=5000
        self.T[self.newVals,0,1:-1] +=zCoeff*self.Coeffs[0,0,1:-1]*(self.T[old,1,1:-1] - self.T[old,0,1:-1])
        self.T[self.newVals,0,:] +=-ConvLoss*self.dt*self.Coeffs[0,0,:]*(self.T[old,0,: ] -t0)
        self.T[self.newVals,-1:,1:-1] +=zCoeff*self.Coeffs[0,-1:,1:-1]*(self.T[old, -2:-1,1:-1] - self.T[old,-1:,1:-1])
        self.T[self.newVals,-1:,:] +=-ConvLoss*self.dt*self.Coeffs[0,-1:,:]*(self.T[old,-1:,:] -t0)

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