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
from POD import*
from time import strftime 
from movie import*



#Make object for application
class App_Window(tkinter.Tk):
    
    def __init__(self,parent):
        tkinter.Tk.__init__(self,parent)
        self.parent = parent
        self.wm_title("Building Sample")
        self.geometry("850x600+100+100")
        self.config(background = "#ecf2f9")
        self.initialize()
    
    def initialize(self):
        self.important_values()     
        self.layout()
        self.menues()        
        self.applyPsy0()
        
    def menues(self):
        
        self.potVar     = tkinter.StringVar(self)
        self.potVar.set("Box")
        potentials = ["Box", "Circle"]
        self.potMenue   = tkinter.OptionMenu(self,self.potVar, *potentials, command=self.potchange)
        self.potMenue.config(bg=color, height=2, width=4)
        self.potMenue["menu"].config(bg=color)
        self.potMenue.place(x=310,y=3)
  
        self.menu = tkinter.Menu(self)
        self.config(menu=self.menu, bg=color)
        
        #Dimension
#        self.subMenuConfig = tkinter.Menu(self.menu, bg=color)
#        self.menu.add_cascade(label="Edit", menu=self.subMenuConfig)
#        self.subMenuConfig.add_command(label = "Configure", command=self.configer)

        
    def potchange(self, var):
        X = self.X
        Y = self.Y
        
        if var=="Box":
            self.V  =(X*X)**(abs(X)/15)+(Y*Y)**(abs(Y)/15)
            self.strPot = "Box"
        elif var=="Circle":
            self.V  = np.sqrt((X*X)**2+(Y*Y)**2)
            self.strPot = "Circle"
        
        vimOld = self.Vim
        self.Vim = self.V/(np.amax(self.V))
        self.im.set_data(self.Vim+self.im.get_array()-vimOld)
        self.canvas.show()
        self.update()
             
        
    def layout(self):
        self.fig = matplotlib.figure.Figure(figsize=(4,4),dpi=140, facecolor="#ecf2f9")
        self.FigSubPlot = self.fig.add_subplot(111)
        self.FigSubPlot.tick_params(labelsize=6)
        firstimage = np.random.rand(self.Ny,self.Nx)
        self.im = self.FigSubPlot.imshow(firstimage, extent=[-self.Nx/20,self.Nx/20,self.Ny/20,-self.Ny/20], cmap = ccmap)


        self.canvas = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.show()
        
#        tkinter.Radiobutton(self,)
#        self.canvas.c
        self.canvas.get_tk_widget().place(y=10,x=0)
        self.canvas.get_tk_widget().configure(background="#ecf2f9", highlightcolor="#ecf2f9", highlightbackground="#ecf2f9")
        self.canvas._tkcanvas.place(y=10,x=0)
        self.canvas.mpl_connect('button_press_event', self.on_click)
        self.canvas.mpl_connect("motion_notify_event",self.on_motion)
        self.canvas.mpl_connect("button_release_event", self.on_release)
        self.canvas.mpl_connect("key_press_event", self.key_press)
#        self.fig.canvas.callbacks.connect('buttonbutton_press_event', self.on_click)
#        self.im.canvas.
        
        tkinter.Label(self,text = "PosX", bg=bLabel, fg=fLabel, font=fontLabel, width=6).place(y=5, x=60)
        self.posX = tkinter.Entry(self, width = 5, bg=bLabel, fg=fLabel)
        self.posX.place(y=5, x=104)
        
        tkinter.Label(self, text = "VarX", bg=bLabel, fg=fLabel, font=fontLabel, width=6).place(y=30, x=60)
        self.varX = tkinter.Entry(self, bg=bLabel, fg=fLabel, width = 5)
        self.varX.place(y=30, x=104)
        
        tkinter.Label(self,text = "PosY", bg=bLabel, fg=fLabel, font=fontLabel, width=6).place(y=5, x=140)
        self.posY = tkinter.Entry(self, bg=bLabel, fg=fLabel, width = 5)
        self.posY.place(y=5, x=181)
        
        tkinter.Label(self,text = "VarY", bg=bLabel, fg=fLabel, font=fontLabel, width=6).place(y=30, x=140)
        self.varY = tkinter.Entry(self, bg=bLabel, fg=fLabel, width = 5)
        self.varY.place(y=30, x=181)
        
        tkinter.Label(self,text = "MomX", bg=bLabel, fg=fLabel, font=fontLabel, width=6).place(y=5, x=225)
        self.momX = tkinter.Entry(self, bg=bLabel, fg=fLabel, width = 5)
        self.momX.place(y=5, x=270)
        
        tkinter.Label(self,text = "MomY", bg=bLabel, fg=fLabel, font=fontLabel, width=6).place(y=30, x=225)
        self.momY = tkinter.Entry(self, bg=bLabel, fg=fLabel, width = 5)
        self.momY.place(y=30, x=270)
        
        button1 = tkinter.Button(self,text="Apply",command=self.applyPsy0, width=7, height=2, bg=bcolor)
        button1.place(x=386,y=5)  
        
        button1 = tkinter.Button(self,text="Random",command=self.randomInput, width=7, height=2, bg=bcolor)
        button1.place(x=455,y=5)  
        
        button1 = tkinter.Button(self,text="Add Particle",command=self.addParticle, width=10, height=2, bg=bcolor)
        button1.place(x=530,y=5) 
        
        self.buttonTest = tkinter.Button(self,text="Test",command=self.tester, width=7, height=2, bg=bcolor)
        self.buttonTest.place(x=620, y=5) 
        
        button3 = tkinter.Button(self,text="NewInput",command=self.newInput, width=7, height=2, bg=bcolor)
        button3.place(x=690, y=5)
        
        button4 = tkinter.Button(self,text="Calculate",command=self.calculate, width=7, height=2, bg=bcolor)
        button4.place(x=760, y=5)
        self.FigSubPlot.add_artist(self.arrowX)
        self.FigSubPlot.add_artist(self.arrowY)
        self.potchange("Box")
        self.resizable(False,False)
        self.update()
                
        self.randomInput()
        
        
    def important_values(self):
        self.Dt = 0.01
        self.Tf = 2
        self.Nx = 2**8  
        self.Ny = 2**8
        self.Dx = 0.1
        self.Dy = 0.1
        self.NModes = 50
        x   =(np.arange(self.Nx)-self.Nx/2)*self.Dx
        y   =(np.arange(self.Ny)-self.Ny/2)*self.Dy
        self.running = False
        self.X, self.Y = np.meshgrid(x,y)
        self.Vim = np.zeros_like(self.X)
        self.V = None
        self.goodapply = False
        self.manual = False
        
        self.changeVarXBool = False
        self.changeVarYBool = False
        self.changePosBool = False
        self.changeMomBool = False
        
        self.arrows = []
        self.arrowcount = 0
        self.particles = []
        self.particlecount = 0
        self.nparticle = 0
        self.ninput = 0
        self.inputlist = []
        self.logbook = []
        self.fontTitle = ('Cambria', 19)
        self.strPot = ""
        self.arrowX = matplotlib.text.Annotation('', fontsize=1, xy=(0, 0),
                                    xytext=(1,1), 
                                    arrowprops=dict(arrowstyle="<-",
                                                    linewidth = 1.,
                                                    alpha =0,
                                                    color = 'orange')
                                                    )
        self.arrowY =  matplotlib.text.Annotation('', fontsize=1, xy=(0, 0),
                                    xytext=(1,1), 
                                    arrowprops=dict(arrowstyle="<-",
                                                    linewidth = 1.,
                                                    alpha =0,
                                                    color = 'orange')
                                                    )
                                                    
        global bLabel, bButton, fLabel, fButton, fontLabel, fontButton, color, bcolor, ccmap, path
        
        path = strftime("%d%b%Y_%H%M") + "/"
        if not os.path.exists(path):
            os.mkdir(path)
            
        ccmap = "Blues"
        bLabel = "#ecf2f9"
        color = "#ecf2f9"
        fLabel = "#29293d"
        bcolor = color
        bButton = "#ecf2f9"
        fontLabel = ('Cambria', 9)
        
    def configer(self):
        edit = Toplevel(self)
        edit.wm_title("Config")
        
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
            self.edt.configure(bg = "color")
        else:
            self.edt.configure(bg = "red")
            goodEntry = False

            

        tf = self.vtf.get()
        if tf>0 and 50>dt:
            self.etf.configure(bg = color)
            self.Tf = tf
        else:
            self.etf.configure(bg = "red")
            goodEntry = False
            
        dx = self.vdx.get()
        if dx>0 and 2>dx:
            self.edx.configure(bg = color)
            self.Dx = dx
        else:
            self.edx.configure(bg = "red")
            goodEntry = False

        dy = self.vdy.get()
        if dy>0 and 2>dy:
            self.edy.config(bg = color)
            self.Dy = dy
        else:
            self.edy.config(bg = "red")
            goodEntry = False
            
        Nx = self.vnx.get()
        if Nx>299 and 1201>Nx:
            self.Nx = Nx
            self.eNx.configure(bg = color)
        else:
            self.eNx.configure(bg = "red")
            goodEntry = False
            
            
        Ny = self.vny.get()
        if Nx>299 and 1201>Nx:
            self.Ny = Ny
            self.eNy.configure(bg = color)
        else:
            self.eNy.configure(bg = "red")
            goodEntry = False 

        NM = self.vnm.get()
        if NM>0 and 150>NM:
            self.NModes = NM
            self.eNMod.configure(bg = color)
        else:
            self.eNMod.configure(bg = "red")
            goodEntry = False

        if goodEntry:
            self.closeWindow(w)
            self.fig = matplotlib.figure.Figure(figsize=(self.Ny/60,self.Nx/100),dpi=100)
            self.FigSubPlot = self.fig.add_subplot(111)
            self.im = self.FigSubPlot.imshow(self.V/10E5, cmap = ccmap)
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
        
        
    def randomInput(self):
        self.posX.delete(0,END)
        rnd = round(np.random.rand()+np.random.randint(-10,10),2)
        self.posX.insert(END, rnd)
        
        self.varX.delete(0,END)
        rnd = round(np.random.rand()+np.random.randint(5),2)
        self.varX.insert(END,rnd)
        
        self.posY.delete(0,END)
        rnd = round(np.random.rand()+np.random.randint(-10,10),2)
        self.posY.insert(END,rnd)
        
        self.varY.delete(0,END)
        rnd = round(np.random.rand()+np.random.randint(5),2)
        self.varY.insert(END,rnd)
        
        self.momX.delete(0,END)
        rnd = round(np.random.rand()+np.random.randint(-5,5),2)
        self.momX.insert(END,rnd)
        
        self.momY.delete(0,END)
        rnd = round(np.random.rand()+np.random.randint(-5,5),2)
        self.momY.insert(END,rnd)

#        self.applyPsy0()

        
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
                self.particles = []
                self.particles.append([PosX, PosY, VarX, VarY, k0x, k0y])
                self.psy0 = gaussian2d(self.X, self.Y, 1, PosX, PosY, VarX, VarY, k0x, k0y)
                psy = abs(self.psy0*self.psy0)
                self.im.set_data(psy+self.Vim)
                if self.arrowcount:
                    for arrow in self.arrows:
                        arrow.remove()
                    self.arrows = []
                    self.arrowcount = 0
                    
                self.arrows.append(matplotlib.text.Annotation('', fontsize=20, xy=(PosX, PosY),
                                    xytext=(PosX+k0x,PosY+k0y), 
                                    arrowprops=dict(arrowstyle="<-",
                                                    linewidth = 1.,
                                                    alpha =0.8,
                                                    color = 'orange')
                                                    ))
                self.arrowcount += 1
                                                    
                self.FigSubPlot.add_artist(self.arrows[0])
                self.canvas.show()
                self.update()
        
        except ValueError:
            self.badInput()
            
        
    def addParticle(self):
        self.goodapply = True
        
        self.arrowY.set_visible(False)
        self.arrowX.set_visible(False)
        try:
            PosX = float(self.posX.get())
            PosY = float(self.posY.get())
            VarX = float(self.varX.get())
            VarY = float(self.varY.get())
            k0x  = float(self.momX.get())
            k0y  = float(self.momY.get())
            
            if self.goodapply:
                self.particles.append([PosX, PosY, VarX, VarY, k0x, k0y])
                self.psy0 += gaussian2d(self.X, self.Y, 1, PosX, PosY, VarX, VarY, k0x, k0y)
                psy = abs(self.psy0*self.psy0)
                self.im.set_data(psy+self.Vim)
#                    if self.arrowcount:
#                        self.arrow.remove()
                    
                self.arrows.append(matplotlib.text.Annotation('', fontsize=20, xy=(PosX, PosY),
                                    xytext=(PosX+k0x,PosY+k0y), 
                                    arrowprops=dict(arrowstyle="<-",
                                                    linewidth = 1.,
                                                    alpha =0.8,
                                                    color = 'orange')
                                                    )
                                            )
                
                self.FigSubPlot.add_artist(self.arrows[self.arrowcount])
                self.arrowcount += 1
                                                    
                self.canvas.show()
                self.update()
        
        except ValueError:
            self.badInput()
            
        
    def rebuildParticles(self):
        
        self.psy0 = 0*self.psy0  
        psy = 0
        for arrow in self.arrows:
            arrow.set_visible(False)
        self.arrows = []
        for particle in self.particles:
            PosX = particle[0]
            PosY = particle[1]
            VarX = particle[2]
            VarY = particle[3]
            k0x  = particle[4]
            k0y  = particle[5]
            
            self.psy0 += gaussian2d(self.X, self.Y, 1, PosX, PosY, VarX, VarY, k0x, k0y)
            psy       += abs(self.psy0*self.psy0)

            self.arrows.append(matplotlib.text.Annotation('', fontsize=20, xy=(PosX, PosY),
                    xytext=(PosX+k0x,PosY+k0y), 
                    arrowprops=dict(arrowstyle="<-",
                                    linewidth = 1.,
                                    alpha =0.8,
                                    color = 'orange')
                                    )
                            )


        self.im.set_data(psy+self.Vim)  
        for arrow in self.arrows:
            self.FigSubPlot.add_artist(arrow)
                                            
        self.canvas.show()
        self.update()
    
    def badInput(self):
        badInPut = Toplevel(self)
        badInPut.wm_title("Bad Input")
        tkinter.Label(badInPut, text = "Bad Input").pack()
        tkinter.Button(badInPut, text = "Ok", command=lambda: self.closeWindow(badInPut)).pack()
        
    def tester(self):
        
        
        self.arrowY.set_visible(False)
        self.arrowX.set_visible(False)
        
        if not self.running:
            self.running = True
            self.buttonTest["text"] = "Stop"
            for arrow in self.arrows:
                arrow.set_visible(False)
            self.FigSubPlot.tick_params(left="off", top="off", right="off", bottom="off")
        else:
            self.running = False
            self.buttonTest["text"] = "Test"
            self.rebuildParticles()
            for arrow in self.arrows:
                arrow.set_visible(True)    
            self.FigSubPlot.tick_params(left="on", top="on", right="on", bottom="on")

        s = schrodinger2d(self.X,self.Y,self.psy0,self.V, self.Dt)
        
        def anim(i):
            if self.running:
                print(i)
                s.snapshot(20)
                self.im.set_data(s.real_psi)
                self.canvas.show()
                self.update()
                return self.im,
         

        animation.FuncAnimation(self.fig, anim,
                              interval=25, blit=True)


    def key_press(self,event):
        if event.key == "delete" and self.manual:
            del self.particles[self.nparticle]
            self.particlecount -= 1
            self.arrowY.set_visible(False)
            self.arrowX.set_visible(False)
            self.nparticle = 0
            self.manual = False
            self.rebuildParticles()
        
    def on_click(self,event):
#        if self.manual:
#            print("so")
##            self.changeVar(event)
#        else:
        if event.dblclick:
            xc = event.xdata
            yc = event.ydata
            self.arrowY.set_visible(False)
            self.arrowX.set_visible(False)
            self.canvas.show()
            self.manual = False
            index = 0
            for particle in self.particles:
                distx = abs(xc-particle[0])
                disty = abs(yc-particle[1])
                if distx<0.4 and disty<0.4:
                    PosX = particle[0]
                    PosY = particle[1]
                    VarX = particle[2]
                    VarY = particle[3]
                    MomX = particle[4]
                    MomY = particle[5]
                    self.posX.delete(0,END)
                    self.posX.insert(END,PosX)
                    self.posY.delete(0,END)
                    self.posY.insert(END,PosY)
                    self.varX.delete(0,END)
                    self.varX.insert(END,VarX)
                    self.varY.delete(0,END)
                    self.varY.insert(END,VarY)
                    self.momX.delete(0,END)
                    self.momX.insert(END,MomX)
                    self.momY.delete(0,END)
                    self.momY.insert(END,MomY)
                    self.addXYvar(VarX,VarY,PosX,PosY)
                    self.nparticle = index
                    self.manual = True
                    self.changeVar(event)
                    break
                index += 1
        else:
            if self.manual:
                self.changeVar(event)
                self.changePos(event)
                self.changeMom(event)
                
    def addXYvar(self,VarX,VarY,PosX,PosY):
        self.startX = PosX-VarX
        self.arrowY.set_visible(False)
        self.arrowX.set_visible(False)
        limX = self.FigSubPlot.get_xlim()[1]
        limY = self.FigSubPlot.get_ylim()[0]
        if self.startX<-limX:
            self.startX=-limX
        self.startY = PosY-VarY
        if self.startY<-limY:
            self.startY=-limY
        self.endX = PosX+VarX
        if self.endX>limX:
            self.endX=limX
        self.endY = PosY+VarY
        if self.endY>limX:
            self.endY=limY
        
        
        self.arrowX = matplotlib.text.Annotation('', fontsize=20, xy=(self.startX,PosY),
                                            xytext=(self.endX,PosY),
                                            arrowprops=dict(arrowstyle="-",
                                                            linewidth = 1.5,
                                                            alpha =0.5,
                                                            color = 'orange')
                                                            )
        self.arrowY = matplotlib.text.Annotation('', fontsize=20, xy=(PosX,self.startY),
                                            xytext=(PosX,self.endY),
                                            arrowprops=dict(arrowstyle="-",
                                                            linewidth = 1.5,
                                                            alpha =0.5,
                                                            color = 'orange')
                                                            )
        if self.changeVarXBool:
            self.varX.delete(0,END)
            self.varX.insert(END,np.round(VarX,decimals=2))
        elif self.changeVarYBool:
            self.varY.delete(0,END)
            self.varY.insert(END,np.round(VarY,decimals=2))
        self.FigSubPlot.add_artist(self.arrowX)
        self.FigSubPlot.add_artist(self.arrowY)
        self.canvas.show()
        
    def changeMom(self,event):
        xc = event.xdata
        yc = event.ydata
        MomXStart = self.particles[self.nparticle][4] + self.particles[self.nparticle][0]
        MomYStart = self.particles[self.nparticle][5] + self.particles[self.nparticle][1]
        distX = abs(xc-MomXStart)
        distY = abs(yc-MomYStart)
        if distX<0.4 and distY<0.4:
            self.changeMomBool = True
    
    def changeVar(self,event):
        xc = event.xdata
        yc = event.ydata
        distX = min([abs(xc-self.startX),abs(xc-self.endX)])
        distY = min([abs(yc-self.startY),abs(yc-self.endY)])
        PosX = self.particles[self.nparticle][0]
        PosY = self.particles[self.nparticle][1]
        
        if distX<0.4 and abs(yc-PosY)<0.4:
            self.changeVarXBool = True
        elif distY<0.4 and abs(abs(xc)-abs(PosX))<0.4:
            self.changeVarYBool = True
            
    def changePos(self,event):
        xc = event.xdata
        yc = event.ydata
        PosX = self.particles[self.nparticle][0]
        PosY = self.particles[self.nparticle][1]
        distX = abs(PosX-xc)
        distY = abs(PosY-yc)
        
        if distX<0.4 and distY<0.4:
            self.changePosBool = True
    
    def addPos(self,xc,yc):
        self.posX.delete(0,END)
        self.posX.insert(END,np.round(xc,decimals=2))
        self.posY.delete(0,END)
        self.posY.insert(END,np.round(yc,decimals=2))
        VarX = self.particles[self.nparticle][2]
        VarY = self.particles[self.nparticle][3]
        self.addXYvar(VarX,VarY,xc,yc)
        
    def addMom(self,xc,yc,PosX,PosY):
        MomX = xc-PosX
        MomY = yc-PosY
        self.momX.delete(0,END)
        self.momX.insert(END,np.round(MomX,decimals=2))
        self.momY.delete(0,END)
        self.momY.insert(END,np.round(MomY,decimals=2))

    def on_motion(self, event):
        xc = event.xdata
        yc = event.ydata
        PosX = self.particles[self.nparticle][0]
        PosY = self.particles[self.nparticle][1]
        if self.changeVarXBool:
            VarX = abs(PosX-xc)
            VarY = self.particles[self.nparticle][3]
            self.addXYvar(VarX,VarY,PosX,PosY)
        elif self.changeVarYBool:
            VarY = abs(PosY-yc)
            VarX = self.particles[self.nparticle][2]
            self.addXYvar(VarX,VarY,PosX,PosY)
        elif self.changePosBool:
            self.addPos(xc,yc)
        elif self.changeMomBool:
            self.addMom(xc,yc,PosX,PosY)
        
    def on_release(self,event):
        if self.changeVarXBool:
            self.particles[self.nparticle][2] = float(self.varX.get())
            self.changeVarXBool = False
            self.rebuildParticles()
        elif self.changeVarYBool:
            self.particles[self.nparticle][3] = float(self.varY.get())
            self.changeVarYBool = False
            self.rebuildParticles()
        elif self.changePosBool:
            self.particles[self.nparticle][0] = float(self.posX.get())
            self.particles[self.nparticle][1] = float(self.posY.get())
            self.changePosBool = False
            self.rebuildParticles()
        elif self.changeMomBool:
            self.changeMomBool = False
#            self.
            self.particles[self.nparticle][4] = float(self.momX.get())
            self.particles[self.nparticle][5] = float(self.momY.get())
            self.rebuildParticles()
        
    def newInput(self):
        
        fig1 = matplotlib.figure.Figure(figsize=(4.7,4.7),dpi=20, facecolor="#ecf2f9")
        FigSubPlot1 = fig1.add_subplot(111)
        FigSubPlot1.axes.get_xaxis().set_visible(False)
        FigSubPlot1.axes.get_yaxis().set_visible(False)
        psy = abs(self.psy0*self.psy0)
        im1 = FigSubPlot1.imshow(psy+self.Vim, extent=[-self.Nx/20,self.Nx/20,self.Ny/20,-self.Ny/20], cmap = ccmap)
        im1.set_data(psy+self.Vim)
        
#        Doing this because of need for biger linewidth
        for particles in self.particles:
            PosX = particles[0]
            PosY = particles[1] 
            k0x  = particles[4]
            k0y  = particles[5]
            
            arrow = matplotlib.text.Annotation('', fontsize=20, xy=(PosX, PosY),
                                xytext=(PosX+k0x,PosY+k0y),
                                arrowprops=dict(arrowstyle="<-",
                                                linewidth = 5.,
                                                alpha =0.8,
                                                color = 'orange')
                                                )
            FigSubPlot1.add_artist(arrow)
        
        shift = self.ninput%12
        
            
        if self.ninput%3==0:
            xc = 520
            yc = 61+int(shift/3)*121
        elif self.ninput%3==1:
            xc = 621
            yc = 61+int(shift/3)*121
        elif self.ninput%3==2:
            xc = 722
            yc = 61+int(shift/3)*121
            

        canvas1 = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(fig1, master=self)
        canvas1.get_tk_widget().place(y=yc, x=xc)
        canvas1.get_tk_widget().configure(background="#ecf2f9", highlightcolor="#ecf2f9", highlightbackground="#ecf2f9")
        canvas1._tkcanvas.place(y=yc, x=xc)
        self.update()
        canvas1.show()
        self.logbook.append([self.ninput, self.V, self.strPot, self.particles])
        print(self.particles)
        self.inputlist.append([self.psy0,self.V])
        self.ninput += 1
        
    def calculate(self):
        self.infoW = Toplevel(self)
        self.infoW.geometry("200x100+100+150")
        self.infoW.wm_title("Calculating")
        self.infoW.config(background = color)
        label1 = tkinter.Label(self.infoW, compound = LEFT, text = "Number of Snapshots per Input", bg=color)
        label1.grid()
        
        self.nSnap = tkinter.Entry(self.infoW, width=10, justify=LEFT)
        self.nSnap.insert(END, 40)
        self.nSnap.grid(row = 0, column = 1)
        
        label1 = tkinter.Label(self.infoW, compound = LEFT, text = "Time between a Snapshot", bg=color)
        label1.grid(row = 1)
        self.btwSnap = tkinter.Entry(self.infoW, width=10, justify=LEFT)
        self.btwSnap.insert(END, 5)
        self.btwSnap.grid(row = 1, column = 1)
        
        label1 = tkinter.Label(self.infoW, compound = LEFT, text = "Create Movie from the Input", bg=color)
        label1.grid(row = 2, column=0)
        self.makeMovie = tkinter.IntVar(self)
        chk = tkinter.Checkbutton(self.infoW, variable=self.makeMovie, onvalue = 1, offvalue = 0)
        chk.grid(row=2, column=1)
        
        
        tkinter.Button(self.infoW, text="Apply", command=self.calculate2, bg=color).grid(row=3, column=1, sticky=E)
        
        
        
        
    def calculate2(self):
        mkmovie = self.makeMovie.get()
        self.infoW.withdraw()
        progWindow = Toplevel(self)
        progWindow.geometry("200x45+200+200")
        progWindow.wm_title("Preparing...")
        progWindow.config(background = color)
        progWindow.update_idletasks()
        label = tkinter.Label(progWindow, text = "Calculating", bg=color)
        label.pack()
        label.update_idletasks()
        bar = Progressbar(progWindow)
        bar.pack()
        bar["value"] = 0
        bar["maximum"] = 100
        bar.update_idletasks()
        snapshots = int(self.nSnap.get())
        betweenSnapshots = int(self.btwSnap.get())
        
        datapath = path + "DataSet/"
        if not os.path.exists(datapath):
            os.mkdir(datapath)
                
        string = "#Input \tPosX \tPosY \tVarX \tVarY \tMomX \tMomY \tV \n"
        file = open((datapath + "Logbook_Input.txt"), "w")
        file.write(string)
        
        for i in self.logbook:
            
            for particle in i[3]:
                string = str(i[0]) + "\t" + str(particle[0]) + "\t" + str(particle[1]) + "\t" + str(particle[2]) + "\t"
                string+= str(particle[3])+ "\t" + str(particle[4]) + "\t" + str(particle[5]) + "\t" + str(i[2]) +"\n"
                file.write(string)
        file.close()

        
        self.A = np.zeros((((self.ninput+1)*snapshots),self.Nx*self.Ny), dtype=np.complex)
        
        Movie_real = np.zeros((self.Ny,self.Nx,snapshots))
        Movie_imag = np.zeros((self.Ny,self.Nx,snapshots))
        Movie_wave = np.zeros((self.Ny,self.Nx,snapshots))
        index = 0
        
        for i in range(self.ninput):
            label["text"] = "Calculating the " +str(i) + "th Input"
            label.update_idletasks()
            psy0 = self.inputlist[i][0]
            V = self.inputlist[i][1]
            schrodinger = schrodinger2d(self.X,self.Y,psy0,V, self.Dt)

                
            for s in range(snapshots):    
                print(s)
                Movie_real[:,:,s] = schrodinger.psi_x.real
                Movie_imag[:,:,s] = schrodinger.psi_x.imag
                Movie_wave[:,:,s] = schrodinger.real_psi
                self.A[index,:] = schrodinger.psi_x.flatten()
                schrodinger.snapshot(betweenSnapshots)
                index += 1  
            if mkmovie:
                label["text"] = "Creating Movies for " +str(i) + "th Input"
                label.update_idletasks()
                
                inputpath = datapath + str(i) + "th Input/" 
                if not os.path.exists(inputpath):
                    os.mkdir(inputpath)
                ani_frame(Movie_real,inputpath+"Real.mp4","Realpart")
                ani_frame(Movie_imag,inputpath+"Imag.mp4","Imaginarypart")
                ani_frame(Movie_wave,inputpath+"AbsSquare.mp4","AbsolutSquare")
            
            bar["value"] = (i+1)*100/(self.ninput)
            bar.update_idletasks()
        
        label["text"] = "Calculating the SVD"
        self.Modes, Atemps, self.sig = SVD_Modes(self.A, self.Nx, self.Ny)
        sig_print(self.sig,datapath+"Sigma.jpg")
        label["text"] = "DONE"
        self.closeWindow(progWindow)
        self.infoW.destroy()
        self.podgo()


        
    def podgo(self):    
        self.hide()
        podFrame = Pod(self.Modes, self.sig, path, self.logbook)
        handler = lambda: self.onCloseOtherFrame(podFrame)
        btn = tkinter.Button(podFrame, text="Close", command=handler, width=10, height=2, bg="#e5ecff")
        btn.place(x=1185, y=5)
        
        
    def hide(self):
        self.withdraw()
        
    def onCloseOtherFrame(self, otherFrame):
        otherFrame.destroy()
        self.destroy()
        

        
if __name__ == "__main__":
    MainWindow = App_Window(None)
    MainWindow.mainloop()

