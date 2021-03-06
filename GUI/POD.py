# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 20:20:19 2015

@author: TinoValentin
"""

import matplotlib
matplotlib.use('TKAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import tkinter
from tkinter import *
from tkinter.ttk import *
import numpy as np
from utilities import*
from scipy.integrate import ode
from movie import*

class Pod(tkinter.Toplevel):
    
    def __init__(self, Modes, sig, path_input, logbook):
        tkinter.Toplevel.__init__(self)
        self.Modes = Modes
        self.sig = sig
        self.maxMode = int(len(sig))
        self.logbook = logbook
        global color, path
        self.evaluation = 0
        color = "#e5ecff"
        path = path_input
        self.particles = []
        self.nparticle = 0
        self.particlecount = 0
        
        for particle in logbook[0][3]:
            self.particlecount += 1
            PosX = particle[0]
            PosY = particle[1]
            VarX = particle[2]
            VarY = particle[3]
            k0x  = particle[4]
            k0y  = particle[5]
            self.particles.append([PosX,PosY,VarX,VarY,k0x,k0y])
        
        self.config(background = color)
        self.wm_title("Evaluating Pod")
        self.geometry("1300x720+100+100")
        self.initialize()

    
    def initialize(self):
        self.important_values()  
        self.menues()    
        self.layout()  
#        self.logentry(0)
        self.potchange("Box")
        
    def menues(self):                
        self.potVar     = tkinter.StringVar(self)
        self.potVar.set("Box")
        potentials = ["Box", "Circle"]
        self.potMenue   = tkinter.OptionMenu(self,self.potVar, *potentials, command=self.potchange)
        self.potMenue.config(bg=color, height=2, width=4)
        self.potMenue["menu"].config(bg=color)
        self.potMenue.place(x=385,y=3)

        self.logVar = tkinter.StringVar(self)

    
        inputs = []
        for i in self.logbook:
            inputs.append(str(i[0])+"th Input")
        
        self.logVar.set("0th Input")
        self.logMenue   = tkinter.OptionMenu(self,self.logVar, *inputs, command=self.logentry)
        self.logMenue.config(bg=color, height=2, width=8)
        self.logMenue["menu"].config(bg=color)
        self.logMenue.place(x=460,y=3)
        
        
        
    def potchange(self, var):
        X = self.X
        Y = self.Y
        
        if var=="Box":
            self.V  =(X*X)**(abs(X)/15)+(Y*Y)**(abs(Y)/15)
            self.strPot = "BoxPotential"
        elif var=="Circle":
            self.V  = np.sqrt((X*X)**2+(Y*Y)**2)
            self.strPot = "CirclePotential"
            
        
        vimOld = self.Vim
        self.Vim = self.V/(np.amax(self.V))
        self.impod.set_data(self.Vim+self.impod.get_array()-vimOld)
        self.canvaspod.show()
        self.update()
             
        
    def layout(self):
        
        self.figpod = matplotlib.figure.Figure(figsize=(4.8,4.8),dpi=140, facecolor=color)
        self.figpodSubPlot = self.figpod.add_subplot(111)
        self.figpodSubPlot.tick_params(labelsize=6, colors="blue")
        firstimage = np.random.rand(self.Ny,self.Nx)
        self.impod = self.figpodSubPlot.imshow(firstimage, extent=[-self.Nx/20,self.Nx/20,self.Ny/20,-self.Ny/20], cmap = "ocean_r")

        self.canvaspod = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(self.figpod, master=self)
        self.canvaspod.show()
        self.canvaspod.get_tk_widget().place(y=50,x=0)
        self.canvaspod.get_tk_widget().configure(background=color, highlightcolor=color, highlightbackground=color)
        self.canvaspod._tkcanvas.place(y=50,x=0)
        self.canvaspod.mpl_connect('button_press_event', self.on_click)
        self.canvaspod.mpl_connect("motion_notify_event",self.on_motion)
        self.canvaspod.mpl_connect("button_release_event", self.on_release)
        self.canvaspod.mpl_connect("key_press_event", self.key_press)
        self.resizable(True,True)
        self.update()
        
        self.figsig = matplotlib.figure.Figure(figsize=(4.6,4.6),dpi=140, facecolor=color)
        self.figsigSubPlot = self.figsig.add_subplot(111)
        self.figsigSubPlot.tick_params(labelsize=6, colors="blue")
        self.figsigSubPlot.set_ylabel(r'$\sigma$', fontsize=9)
        self.figsigSubPlot.set_xlabel(r'$\aleph \sigma$', fontsize=9)
        self.imsig = self.figsigSubPlot.plot(self.sig)
        
        self.canvassig = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(self.figsig, master=self)
        self.canvassig.show()
        self.canvassig.get_tk_widget().place(y=58,x=650)
        self.canvassig.get_tk_widget().configure(background=color, highlightcolor=color, highlightbackground=color)
        self.canvassig._tkcanvas.place(y=58,x=650)

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
        
#        tkinter.Label(self,text = "How\n many\n Modes", bg=bLabel, fg=fLabel, font=fontLabel, width=6).place(y=3, x=305)
#        self.nModes = tkinter.Entry(self, bg=bLabel, fg=fLabel, width = 3)
#        rnd = np.random.randint(1,10)
#        self.nModes.insert(END,rnd)
#        self.nModes.place(y=17, x=355)
        
        button0 = tkinter.Button(self,text="Add Particle", bg=bLabel, fg=fLabel,command=self.addParticle, width=7, height=2)
        button0.place(x=315, y=5)  
        
        button1 = tkinter.Button(self,text="Apply", bg=bLabel, fg=fLabel,command=self.applyPsy0, width=7, height=2)
        button1.place(x=560, y=5)  
             
        button5 = tkinter.Button(self,text="Show Sigma", bg=bLabel, fg=fLabel,command=self.showSig, width=10, height=2)
        button5.place(x=730, y=5) 
        
        button5 = tkinter.Button(self,text="Show Mode", bg=bLabel, fg=fLabel,command=self.showMod, width=10, height=2)
        button5.place(x=820, y=5)
        self.mode = tkinter.Entry(self, bg=bLabel, fg=fLabel, width = 5)
        self.mode.place(x=910, y=15)
          
        button4 = tkinter.Button(self,text="Calculate", bg="#d3ddf9", fg=fLabel,command=self.calculate, width=9, height=2)
        button4.place(x=640, y=5)
        
        self.resizable(True,True)

        
        
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
        self.psy0 = gaussian2d(self.X, self.Y, 1, 1, 1, 1, 1, 1, 1)
        self.Vim = np.zeros_like(self.X)
        self.V = None
        self.arrowcount = 0
        self.arrows = []
        self.ninput = 0
        self.inputlist = []
        self.fontTitle = ('Cambria', 19)

        self.goodapply = False
        self.manual = False
        
        self.changeVarXBool = False
        self.changeVarYBool = False
        self.changePosBool = False
        self.changeMomBool = False
        


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
                                                    
        
        global bLabel, bButton, fLabel, fButton, fontLabel, fontButton
        color = "#e5ecff"
        bLabel = color
        bButton = color
        fLabel = "black"
        fontLabel = ('Cambria', 9)
        

    def closeWindow(self,w):
        w.destroy()
             
    def refreshFigure(self,x,y):
        self.line1.set_data(x,y)
        ax = self.canvaspod.figure.axes[0]
        ax.set_xlim(x.min(), x.max())
        ax.set_ylim(y.min(), y.max())        
        self.canvaspod.draw()
        
#    def changePot(index, value, op):
#        print("HELLO")
        
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
                self.impod.set_data(psy+self.Vim)
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
                                                    
                self.figpodSubPlot.add_artist(self.arrows[0])
                self.canvaspod.show()
                self.update()
        
        except ValueError:
            self.badInput()
            
    def logentry(self, var):
        entry = int(var[0])
        log   = self.logbook[entry]
        self.particles = []
        self.particlecount = 0
        self.psy0 = 0*self.psy0
        psy = 0*abs(self.psy0*self.psy0)
        
        for arrow in self.arrows:
            arrow.set_visible(False)
        self.arrows = []
        
        for particle in log[3]:
            self.particlecount += 1
            PosX = particle[0]
            PosY = particle[1]
            VarX = particle[2]
            VarY = particle[3]
            k0x  = particle[4]
            k0y  = particle[5]
            self.particles.append([PosX,PosY,VarX,VarY,k0x,k0y])
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


        self.impod.set_data(psy+self.Vim)  
        for arrow in self.arrows:
            self.figpodSubPlot.add_artist(arrow)
                                            
        self.canvaspod.show()
        self.update()

            
#            
#            
#        strPot = log[2]
#        
#        self.posX.delete(0,END)
#        self.posX.insert(END,PosX)
#        self.posY.delete(0,END)
#        self.posY.insert(END,PosY)
#        self.varX.delete(0,END)
#        self.varX.insert(END,VarX)
#        self.varY.delete(0,END)
#        self.varY.insert(END,VarY)        
#        self.momX.delete(0,END)
#        self.momX.insert(END,k0x)
#        self.momY.delete(0,END)
#        self.momY.insert(END,k0y)
#        self.V = log[1]
#        self.potVar.set(str(strPot))
#        
#        self.psy0 = gaussian2d(self.X, self.Y, 1, PosX, PosY, VarX, VarY, k0x, k0y)
#        psy = abs(self.psy0*self.psy0)
#        self.impod.set_data(psy+self.Vim)
#        for arrow in self.arrows:
#            
#        if self.arrowexist:
#            self.arrow.remove()
#            
#        self.arrow = matplotlib.text.Annotation('', fontsize=20, xy=(PosX, PosY),
#                            xytext=(PosX+k0x,PosY+k0y), 
#                            arrowprops=dict(arrowstyle="<-",
#                                            linewidth = 1.,
#                                            alpha =0.8,
#                                            color = 'blue')
#                                            )
#        self.arrowexist = True
#                                            
#        self.figpodSubPlot.add_artist(self.arrow)
#        self.potchange(strPot)
#        self.canvaspod.show()
#        self.update()
        
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
                
                self.figpodSubPlot.add_artist(self.arrows[self.arrowcount])
                self.arrowcount += 1
                                                    
                self.canvaspod.show()
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


        self.impod.set_data(psy+self.Vim)  
        for arrow in self.arrows:
            self.figpodSubPlot.add_artist(arrow)
                                            
        self.canvaspod.show()
        self.update()
        
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

        if event.dblclick:
            xc = event.xdata
            yc = event.ydata
            self.arrowY.set_visible(False)
            self.arrowX.set_visible(False)
            self.canvaspod.show()
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
        limX = self.figpodSubPlot.get_xlim()[1]
        limY = self.figpodSubPlot.get_ylim()[0]
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
        self.figpodSubPlot.add_artist(self.arrowX)
        self.figpodSubPlot.add_artist(self.arrowY)
        self.canvaspod.show()
        
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
#        print("Index ", self.nparticle)
#        print("Länge ", len(self.particles))
#        print("Particle ", self.particles[self.nparticle])
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
            self.particles[self.nparticle][4] = float(self.momX.get())
            self.particles[self.nparticle][5] = float(self.momY.get())
            self.rebuildParticles()
        
            
    def showSig(self):
        self.figsig = matplotlib.figure.Figure(figsize=(4.6,4.6),dpi=140, facecolor=color)
        self.figsigSubPlot = self.figsig.add_subplot(111)
        self.figsigSubPlot.tick_params(labelsize=6, colors="blue")
        self.figsigSubPlot.set_ylabel(r'$\sigma$', fontsize=9)
        self.figsigSubPlot.set_xlabel(r'$\aleph \sigma$', fontsize=9)
#        self.figsigSubPlot.label_params()
        self.imsig = self.figsigSubPlot.plot(self.sig)
        
        self.canvassig = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(self.figsig, master=self)
        self.canvassig.show()
        self.canvassig.get_tk_widget().place(y=58,x=650)
        self.canvassig.get_tk_widget().configure(background=color, highlightcolor=color, highlightbackground=color)
        self.canvassig._tkcanvas.place(y=58,x=650)
        self.canvassig.tag_lower()
        
    def showMod(self):
        nMode = int(self.mode.get())
        
        self.figsig = matplotlib.figure.Figure(figsize=(4.8,4.8),dpi=140, facecolor=color)
        self.figsigSubPlot = self.figsig.add_subplot(111)
        self.figsigSubPlot.tick_params(labelsize=6, colors="blue")
        self.imsig = self.figsigSubPlot.imshow(self.Modes[nMode].real, extent=[-self.Nx/20,self.Nx/20,self.Ny/20,-self.Ny/20], cmap = "ocean_r")

        
        self.canvassig = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(self.figsig, master=self)
        self.canvassig.show()
        self.canvassig.get_tk_widget().place(y=50,x=650)
        self.canvassig.get_tk_widget().configure(background=color, highlightcolor=color, highlightbackground=color)
        self.canvassig._tkcanvas.place(y=50,x=650)
        
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

        self.applyPsy0()        
        
        
    def calculate(self):
        
        self.infoW = Toplevel(self)
        self.infoW.geometry("200x100+100+150")
        self.infoW.wm_title("Calculating")
        self.infoW.config(background = color)
        label1 = tkinter.Label(self.infoW, compound = LEFT, text = "Number of Snapshots", bg=color)
        label1.grid()
        
        self.nSnap = tkinter.Entry(self.infoW, width=10, justify=LEFT)
        self.nSnap.insert(END, 40)
        self.nSnap.grid(row = 0, column = 1)
        
        label1 = tkinter.Label(self.infoW, compound = LEFT, text = "Time between a Snapshot", bg=color)
        label1.grid(row = 1)
        self.btwSnap = tkinter.Entry(self.infoW, width=10, justify=LEFT)
        self.btwSnap.insert(END, 5)
        self.btwSnap.grid(row = 1, column = 1)
        
        tkinter.Label(self.infoW,text = "How\n many\n Modes", bg=bLabel, fg=fLabel, font=fontLabel, width=6).grid(row = 2)
        self.nModes = tkinter.Entry(self.infoW, bg=bLabel, fg=fLabel, width = 3)
        rnd = np.random.randint(1,10)
        self.nModes.insert(END,rnd)
        self.nModes.grid(row=2, column=1)
        
        tkinter.Button(self.infoW, text="Apply", command=self.calculate2, bg=color).grid(row=3, column=1, sticky=E)
        
        
    def calculate2(self):
        self.infoW.withdraw()
        progWindow = Toplevel(self)
        progWindow.geometry("200x45+500+500")
        progWindow.wm_title("Calculating")
        progWindow.config(background = color)
        progWindow.update_idletasks()
        label = tkinter.Label(progWindow, text = "Preparing...", bg=color)
        label.pack()
        label.update_idletasks()
        bar = Progressbar(progWindow)
        bar.pack()
        bar["value"] = 0
        bar["maximum"] = 100
        bar.update_idletasks()
        
        evalpath = path + str(self.evaluation) + "th evaluation" + "/"
        if not os.path.exists(evalpath):
            os.mkdir(evalpath)
            
#        evaluation += 1
        
        nModes = int(self.nModes.get())
        if nModes>self.maxMode:
            print("Input for Modes is too big - taking Maximum: ", self.maxMode, " instead.")
            nModes = self.maxMode

        string = "#Particle \tPosX \tPosY \tVarX \tVarY \tMomX \tMomY \tV \t\t#Modes \n"
        print(string)
        file = open((evalpath + "Logbook.txt"), "w")
        file.write(string)
        pnum = 0
        for particle in self.particles:
            
            PosX = float(self.posX.get())
            PosY = float(self.posY.get())
            VarX = float(self.varX.get())
            VarY = float(self.varY.get())
            k0x  = float(self.momX.get())
            k0y  = float(self.momY.get())
            Modes = self.Modes
    
            string =str(pnum) + "\t" + str(PosX) + "\t" + str(PosY) + "\t" + str(VarX) + "\t" + str(VarY) 
            string+= "\t" + str(k0x) + "\t" + str(k0y) + "\t" + str(self.strPot) + "\t" + str(nModes) + "\n"
            file.write(string)
            pnum += 1
            
        file.close()
        
        
        V = self.V
        Nx = self.Nx
        Ny = self.Ny
        dx = self.Dx
        dy = self.Dy
        hbar = 1
        m = 1
        t0 = 0
        dt = 0.01
        steps = int(self.nSnap.get())
        betweenSnapshots = int(self.btwSnap.get())
        Tf = dt*steps*betweenSnapshots
        
        psy0 = self.psy0
        a0 = self.getA0(nModes)
        
        
        # Calculating the Laplacian of every POD-Mode   
        DModes = np.zeros((nModes,Ny,Nx), np.complex)   
        Coef = np.zeros((nModes,Ny,Nx), np.complex)   
            
        
        for i in range(nModes):
            DModes[i,:,:] = laplacian(Modes[i,:,:],dx,dy)
            Coef[i,:,:]   = (DModes[i,:,:]*hbar/(2*m)-V*Modes[i,:,:]/hbar)
        
        H = np.zeros((nModes,nModes), dtype=np.complex)
        
        for k in range(nModes):
            for l in range(nModes):
                H[k,l] = (np.conjugate(Modes[k])*Coef[l]).sum()
                
                
        def fun(t,y):
            """
            returns:        da_k/dt =  Σ_i a_i H[k,i] (here with hbar=1)
            the summation of i is performed in the middle loop
            the integral is performed in the last line before returning alpha
            alpha in the return is a vektor of length cut
            """
    
            alpha = np.zeros(nModes, dtype=np.complex) 
            
            for k in range(nModes):                
                temp = 0+0j
                for i in range(nModes):
                    temp   += y[i]*H[k,i]
                    
                alpha[k] = 1j*temp
            return alpha
        
        solution = []
        POD_A = []

# Solving Ode    
        r = ode(fun).set_integrator('zvode', method='bdf', with_jacobian=False)
        r.set_initial_value(a0, t0)
        
        poddt = betweenSnapshots*dt
        
        label["text"] = "Start with Galerkin Ode"
        bar["value"] = 0
        bar.update_idletasks()
        label.update_idletasks()
        
        while r.successful() and r.t < Tf:
            POD_A.append(r.y)
            r.integrate(r.t+poddt)
            solution.append([r.t,r.y])
        
        label["text"] = "Calculating Leapfrog"
        bar["value"] = 33
        bar.update_idletasks()
        label.update_idletasks()
        
        LeapfrogMovie_real = np.zeros((Ny,Nx,steps))
        LeapfrogMovie_imag = np.zeros((Ny,Nx,steps))
        LeapfrogMovie_wave = np.zeros((Ny,Nx,steps))
        schrodinger = schrodinger2d(self.X,self.Y,psy0,V, self.Dt)
        
        for s in range(steps):  
            print(s)
            LeapfrogMovie_real[:,:,s] = schrodinger.psi_x.real
            LeapfrogMovie_imag[:,:,s] = schrodinger.psi_x.imag
            LeapfrogMovie_wave[:,:,s] = schrodinger.real_psi
            schrodinger.snapshot(betweenSnapshots)
            
        
#        POD_Atemps = np.zeros((snapshots,cut), dtype=np.complex)
        
        ModesFlat = np.zeros((nModes,Nx*Ny), dtype=np.complex)
        for i in range(nModes):
            ModesFlat[i] = Modes[i].flatten()
        

        PodMovie_real = np.zeros((Ny,Nx,steps))
        PodMovie_imag = np.zeros((Ny,Nx,steps))
        PodMovie_wave = np.zeros((Ny,Nx,steps))
        
        DifMovie_real = np.zeros((Ny,Nx,steps))
        DifMovie_imag = np.zeros((Ny,Nx,steps))
        DifMovie_wave = np.zeros((Ny,Nx,steps))
        
        Ap = np.dot(POD_A,ModesFlat)
        

        for i in range(steps):
            psi           = np.reshape(Ap[i,:],(Ny,Nx))
            PodMovie_real[:,:,i] = psi.real
            PodMovie_imag[:,:,i] = psi.imag
            PodMovie_wave[:,:,i] = abs(psi*psi)
            
            DifMovie_real[:,:,i] = LeapfrogMovie_real[:,:,i] - PodMovie_real[:,:,i]
            DifMovie_imag[:,:,i] = LeapfrogMovie_imag[:,:,i] - PodMovie_imag[:,:,i]
            DifMovie_wave[:,:,i] = LeapfrogMovie_wave[:,:,i] - PodMovie_wave[:,:,i]
#            print("POD", PodMovie_real[:,:,i])
#            print("Leap", LeapfrogMovie_real[:,:,i])
#            print("Dif", DifMovie_real[:,:,i])
            
#        path = ""
        label["text"] = "Creating POD Movie"
        bar["value"] = 66
        bar.update_idletasks()
        label.update_idletasks()
        
        ani_frame(PodMovie_real,evalpath+"POD-Real.mp4", "POD Realpart")
        ani_frame(PodMovie_imag,evalpath+"POD-Imag.mp4", "POD Imaginarypart")
        ani_frame(PodMovie_wave,evalpath+"POD-AbsSquare.mp4", "POD AbsolutSquare")
        
        label["text"] = "Creating LeapFrog Movie"
        bar["value"] = 77
        bar.update_idletasks()
        label.update_idletasks()

        ani_frame(LeapfrogMovie_real,evalpath+"Leapfrog-Real.mp4", "Leapfrog Realpart")
        ani_frame(LeapfrogMovie_imag,evalpath+"Leapfrog-Imag.mp4", "LeapfrogImaginarypart")
        ani_frame(LeapfrogMovie_wave,evalpath+"Leapfrog-AbsSquare.mp4", "Leapfrog AbsolutSquare")
        
        label["text"] = "Creating the Difference Movie"
        bar["value"] = 88
        bar.update_idletasks()
        label.update_idletasks()

        ani_frame(DifMovie_real,evalpath+"Dif-Real.mp4", "Dif Realpart")
        ani_frame(DifMovie_imag,evalpath+"Dif-Imag.mp4", "Dif Imaginarypart")
        ani_frame(DifMovie_wave,evalpath+"Dif-AbsSquare.mp4", "Dif AbsolutSquare")
        
        self.evaluation += 1        
        
        label["text"] = "Done"
        bar["value"] = 100
        bar.update_idletasks()
        label.update_idletasks()
        self.infoW.destroy()
        progWindow.destroy()
        
        
        
    def getA0(self,nModes):
        a0 = np.zeros(nModes, dtype=np.complex)
        
        for i in range(nModes):
            a0[i]  = np.sum(self.psy0*np.conjugate(self.Modes[i]))
        
        return a0
        
        
        
    
if __name__ == "__main__":
    Modes = np.random.rand(10,256,256)
    sig = np.random.rand(30)
    path_input = "hier"
    logbook = []
    logbook.append([1,"sa", 2, [1,2,3,4,5,4]])
    MainWindow = Pod(Modes, sig, path_input, logbook)
    MainWindow.mainloop()