# import modules that I'm using
import matplotlib
matplotlib.use('TKAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib import animation
import tkinter
from tkinter import *
from tkinter.ttk import *
import numpy as np
from utilities import*




#Make object for application
class App_Window(tkinter.Tk):
    
    dx = 300
    

    def __init__(self,parent):
        
        tkinter.Tk.__init__(self,parent)
        self.parent = parent
        self.geometry("800x700+100+100")
        self.initialize()
    
    
    def initialize(self):


        self.important_values()     
        self.layout()
        self.menues()   
        
        
    def menues(self):
        self.menu = tkinter.Menu(self)
        self.config(menu=self.menu)
        
        #Dimension
        self.subMenuConfig = tkinter.Menu(self.menu)
        self.menu.add_cascade(label="Edit", menu=self.subMenuConfig)
        self.subMenuConfig.add_command(label = "Configure", command=self.configer)
        
        #Potential
        self.subMenuPot = tkinter.Menu(self.menu)
        self.menu.add_cascade(label="Potential", menu=self.subMenuPot)
        self.subMenuPot.add_command(label="Box", command=self.potBox) 
        
        
    def layout(self):
        tkinter.Label(self,text = "PosX").grid()
        self.posX = tkinter.Entry(self)
        rnd = round(np.random.rand(),2)+np.random.randint(-5,5)
        self.posX.insert(END, rnd)
        self.posX.grid(row=0, column=1)
        
        tkinter.Label(self, text = "VarX").grid(row=1)
        self.varX = tkinter.Entry(self)
        rnd = round(np.random.rand(),2)+np.random.randint(3)
        self.varX.insert(END,rnd)
        self.varX.grid(row=1, column=1)
        
        tkinter.Label(self,text = "PosY").grid(row=0, column=2)
        self.posY = tkinter.Entry(self)
        rnd = round(np.random.rand(),2)+np.random.randint(-5,5)
        self.posY.insert(END,rnd)
        self.posY.grid(row=0, column=3)
        
        tkinter.Label(self, text = "VarY").grid(row=1, column=2)
        self.varY = tkinter.Entry(self)
        rnd = round(np.random.rand(),2)+np.random.randint(3)
        self.varY.insert(END,rnd)
        self.varY.grid(row=1, column=3)
        
        tkinter.Label(self,text = "MomentumX").grid(row=0, column=4)
        self.momX = tkinter.Entry(self)
        rnd = round(np.random.rand(),2)+np.random.randint(-5,5)
        self.momX.insert(END,rnd)
        self.momX.grid(row=0, column=5)
        
        tkinter.Label(self, text = "MomentumY").grid(row=1, column=4)
        self.momY = tkinter.Entry(self)
        rnd = round(np.random.rand(),2)+np.random.randint(-5,5)
        self.momY.insert(END,rnd)
        self.momY.grid(row=1, column=5)
        
        self.button1 = tkinter.Button(self,text="Apply",command=self.applyPsy0, width=10)
        self.button1.grid(row=3, column=2, sticky=W)     
        self.button2 = tkinter.Button(self,text="Test",command=self.tester, width=10)
        self.button2.grid(row=3, column=3, sticky=W)
        self.button3 = tkinter.Button(self,text="NewInput",command=self.OnButton1Click, width=10)
        self.button3.grid(row=3, column=4, sticky=W)
        
        self.fig = matplotlib.figure.Figure(figsize=(self.Ny/60,self.Nx/100),dpi=100)
        FigSubPlot = self.fig.add_subplot(111)
        firstimage = np.random.rand(self.Ny,self.Nx)
        self.im = FigSubPlot.imshow(firstimage)

        self.canvas = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(row=4, columnspan=8)
        self.canvas._tkcanvas.grid(row=4)
        self.resizable(True,False)
        self.update()
        
        
    def important_values(self):
        self.Dt = 0.01
        self.Tf = 2
        self.Nx = 512
        self.Ny = 512
        self.Dx = 0.1
        self.Dy = 0.1
        self.NModes = 50
        x   =(np.arange(self.Nx)-self.Nx/2)*self.Dx
        y   =(np.arange(self.Ny)-self.Ny/2)*self.Dy
        
        self.X, self.Y = np.meshgrid(x,y)
        self.Vim = np.zeros_like(self.X)
        
        self.fontTitle = ('Cambria', 19)
        
    def potBox(self):
        Nx = self.Nx
        Ny = self.Ny
        dx = self.Dx
        dy = self.Dy
        x   =(np.arange(Nx)-Nx/2)*dx
        y   =(np.arange(Ny)-Ny/2)*dy
        X,Y = np.meshgrid(x,y)
        boundX = X[0,0]*7.0/10
        boundY = Y[0,0]*7.0/10
        self.V = potential(self.X,self.Y,0.7,1E6)
        self.Vim = self.V/1E6
        self.im.set_data(self.Vim+self.im.get_array())
        self.canvas.show()
        self.update()
        
    def configer(self):
        edit = Toplevel(self)
        edit.wm_title("Config")
#        dimchange.geometry("400x200+100+100")
        
#        Time
        tkinter.Label(edit, text="Time Config", font=self.fontTitle).grid(row=0, columnspan=6)
        tkinter.Label(edit, text="Dt").grid(row=1, sticky=E)
        self.vdt = DoubleVar(edit, value=self.Dt)
        self.edt = tkinter.Entry(edit, width=10, justify=RIGHT, textvariable=self.vdt)
        self.edt.grid(row=1, column=1)
        tkinter.Label(edit, text="(0,5]").grid(row=1,column=2, sticky=W)
        
        tkinter.Label(edit, text="Tfinal").grid(row=1, column=3, sticky=E)
        self.vtf = DoubleVar(edit, value=self.Tf)
        self.etf = tkinter.Entry(edit, width=10, justify=RIGHT, textvariable=self.vtf)
        self.etf.grid(row=1, column=4)
        rldt = tkinter.Label(edit, text="(0,50]").grid(row=1,column=5, sticky=W)
        
        #Coordinates
        tkinter.Label(edit, text="Cordinate Config", font=self.fontTitle).grid(row=3, columnspan=6)
        tkinter.Label(edit, text="Dx").grid(row=4, sticky=E)
        self.vdx = DoubleVar(edit, value=self.Dx)
        self.edx = tkinter.Entry(edit, width=10, justify=RIGHT, textvariable=self.vdx)
        self.edx.grid(row=4, column=1)
        tkinter.Label(edit, text="(0-2]").grid(row=4,column=2, sticky=W)

        tkinter.Label(edit, text="Dy").grid(row=5, sticky=E)
        self.vdy = DoubleVar(edit, value=self.Dy)
        self.edy = tkinter.Entry(edit, width=10, justify=RIGHT, textvariable=self.vdy)
        self.edy.grid(row=5, column=1)
        tkinter.Label(edit, text="(0-2]").grid(row=5,column=2, sticky=W)
        
        tkinter.Label(edit, text="Nx").grid(row=4, column=3, sticky=E)
        self.vnx = DoubleVar(edit, value=self.Nx)
        self.eNx = tkinter.Entry(edit, width=10, justify=RIGHT, textvariable=self.vnx)
        self.eNx.grid(row=4, column=4)
        tkinter.Label(edit, text="[300-1200]").grid(row=4,column=5, sticky=W)
        
        tkinter.Label(edit, text="Ny").grid(row=5, column=3, sticky=E)
        self.vny = DoubleVar(edit, value=self.Ny)
        self.eNy = tkinter.Entry(edit, width=10, justify=RIGHT, textvariable=self.vny)
        self.eNy.grid(row=5, column=4)
        tkinter.Label(edit, text="[300-1200]").grid(row=5,column=5, sticky=W)
        
        # Print Modes
        tkinter.Label(edit, text="Print Modes", font=self.fontTitle).grid(row=6, columnspan=3)
        tkinter.Label(edit, text="N Modes").grid(row=7,column=0, sticky=E)
        self.vnm = DoubleVar(edit, value=self.NModes)
        self.eNMod = tkinter.Entry(edit, width=10, justify=RIGHT, textvariable=self.vnm)
        self.eNMod.grid(row=7, column=1)
        tkinter.Label(edit, text="[0-150]").grid(row=7,column=2, sticky=W)
        
        # Print Modes
        tkinter.Label(edit, text="Print Modes", font=self.fontTitle).grid(row=6, columnspan=3)
        tkinter.Label(edit, text="N Modes").grid(row=7,column=3, sticky=E)
        
        self.evar = tkinter.Entry(edit, width=10, justify=RIGHT)
        self.evar.grid(row=7, column=4)
        tkinter.Label(edit, text="[0-150]").grid(row=7,column=5, sticky=W)
        
        tkinter.Button(edit, text="cancel", command =lambda: self.closeWindow(edit)).grid(row=9, column = 4)
        tkinter.Button(edit, text="apply", command =lambda: self.applyEdit(edit)).grid(row=9, column = 5)

        



    def closeWindow(self,w):
        w.destroy()
        
    def applyEdit(self,w):
        goodEntry = True
        
        dt = self.vdt.get()
        if dt>0 and 5>dt:
            self.Dt = dt
            self.edt.configure(bg = "white")
        else:
            self.edt.configure(bg = "red")
            goodEntry = False

            

        tf = self.vtf.get()
        if tf>0 and 50>dt:
            self.etf.configure(bg = "white")
            self.Tf = tf
        else:
            self.etf.configure(bg = "red")
            goodEntry = False
            
        dx = self.vdx.get()
        if dx>0 and 2>dx:
            self.edx.configure(bg = "white")
            self.Dx = dx
        else:
            self.edx.configure(bg = "red")
            goodEntry = False

        dy = self.vdy.get()
        if dy>0 and 2>dy:
            self.edy.config(bg = "white")
            self.Dy = dy
        else:
            self.edy.config(bg = "red")
            goodEntry = False
            
        Nx = self.vnx.get()
        if Nx>299 and 1201>Nx:
            self.Nx = Nx
            self.eNx.configure(bg = "white")
        else:
            self.eNx.configure(bg = "red")
            goodEntry = False
            
            
        Ny = self.vny.get()
        if Nx>299 and 1201>Nx:
            self.Ny = Ny
            self.eNy.configure(bg = "white")
        else:
            self.eNy.configure(bg = "red")
            goodEntry = False 

        NM = self.vnm.get()
        if NM>0 and 150>NM:
            self.NModes = NM
            self.eNMod.configure(bg = "white")
        else:
            self.eNMod.configure(bg = "red")
            goodEntry = False

        if goodEntry:
            self.closeWindow(w)
            self.fig = matplotlib.figure.Figure(figsize=(self.Ny/60,self.Nx/100),dpi=100)
            FigSubPlot = self.fig.add_subplot(111)
            self.im = FigSubPlot.imshow(self.V/10E5)
            self.fig.colorbar(self.im)
            self.canvas = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(self.fig, master=self)
            self.canvas.show()
            self.canvas.get_tk_widget().grid(row=4, columnspan=8)
            self.canvas._tkcanvas.grid(row=4)
            self.resizable(False,False)
            self.update()
            
   
             
    def refreshFigure(self,x,y):
        self.line1.set_data(x,y)
        ax = self.canvas.figure.axes[0]
        ax.set_xlim(x.min(), x.max())
        ax.set_ylim(y.min(), y.max())        
        self.canvas.draw()
        
    def changePot(index, value, op):
        print("HELLO")
        
    def applyPsy0(self):
        self.goodapply = True
        
        try:
            PosX = float(self.posX.get())
            PosY = float(self.posY.get())
            VarX = float(self.varX.get())
            VarY = float(self.varY.get())
            k0x  = float(self.momX.get())
            k0y  = float(self.momY.get())
            
            if self.goodapply:
                self.psy0 = gaussian2d(self.X, self.Y, 1, PosX, PosY, VarX, VarY, k0x, k0y)
                psy = abs(self.psy0*self.psy0)
                self.im.set_data(psy+self.Vim)
                self.canvas.show()
                self.update()
        
        except ValueError:
            self.badInput()
        
    
    def badInput(self):
        badInPut = Toplevel(self)
        badInPut.wm_title("Bad Input")
        tkinter.Label(badInPut, text = "Bad Input").pack()
        tkinter.Button(badInPut, text = "Ok", command=lambda: self.closeWindow(badInPut)).pack()
        
    def tester(self):
        
#        print(self.V)
        X = self.X/2
        Y = self.Y/2
        
        s = schrodinger2d(self.X,self.Y,self.psy0,self.V, self.Dt)
        
        def anim(i):
            print(i)
            s.snapshot(20)
            self.im.set_data(s.real_psi)
            self.canvas.show()
            self.update()
            return self.im,
            
        T = np.arange(0,4,0.1)
        ani = animation.FuncAnimation(self.fig, anim,
                              interval=25, blit=True)
        
    
    def OnButton1Click(self):
        # file is opened here and some data is taken
        # I've just set some arrays here so it will compile alone
        x=[]
        y=[]
        for num in range(0,1000):x.append(num*.001+1)
        # just some random function is given here, the real data is a UV-Vis spectrum
        for num2 in range(0,1000):y.append(sc.math.sin(num2*.06)+sc.math.e**(num2*.001))
        X = np.array(x)
        Y = np.array(y)
        self.refreshFigure(X,Y)

if __name__ == "__main__":
    MainWindow = App_Window(None)
    MainWindow.mainloop()

