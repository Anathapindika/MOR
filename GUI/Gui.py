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
        self.geometry("750x600+100+100")
        self.config(background = "#ecf2f9")
        self.initialize()
    
    def initialize(self):
        self.important_values()     
        self.layout()
        self.menues()   
        
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
        
        self.buttonTest = tkinter.Button(self,text="Test",command=self.tester, width=7, height=2, bg=bcolor)
        self.buttonTest.place(x=530, y=5) 
        
        button3 = tkinter.Button(self,text="NewInput",command=self.newInput, width=7, height=2, bg=bcolor)
        button3.place(x=600, y=5)
        
        button4 = tkinter.Button(self,text="Calculate",command=self.calculate, width=7, height=2, bg=bcolor)
        button4.place(x=670, y=5)
        

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
        self.arrowexist = False
        self.ninput = 0
        self.inputlist = []
        self.logbook = []
        self.fontTitle = ('Cambria', 19)
        self.strPot = ""
        
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

        self.applyPsy0()

        
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
                if self.arrowexist:
                    self.arrow.remove()
                    
                self.arrow = matplotlib.text.Annotation('', fontsize=20, xy=(PosX, PosY),
                                    xytext=(PosX+k0x,PosY+k0y), 
                                    arrowprops=dict(arrowstyle="<-",
                                                    linewidth = 1.,
                                                    alpha =0.8,
                                                    color = 'orange')
                                                    )
                self.arrowexist = True
                                                    
                self.FigSubPlot.add_artist(self.arrow)
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
        
        if self.arrowexist:
            self.arrow.remove()
            self.arrowexist = False

        
        if not self.running:
            self.running = True
            self.buttonTest["text"] = "Stop"
        else:
            self.running = False
            self.buttonTest["text"] = "Test"
            self.applyPsy0()
            

        s = schrodinger2d(self.X,self.Y,self.psy0,self.V, self.Dt)
        
        def anim(i):
            if self.running:
                print(i)
                s.snapshot(20)
                self.im.set_data(s.real_psi)
                self.canvas.show()
                self.update()
                return self.im,
         

        T = np.arange(0,4,0.1)
        
        ani = animation.FuncAnimation(self.fig, anim,
                              interval=25, blit=True)
        
    
    def newInput(self):
        self.applyPsy0()
        self.fig1 = matplotlib.figure.Figure(figsize=(4.7,4.7),dpi=20, facecolor="#ecf2f9")
        self.FigSubPlot1 = self.fig1.add_subplot(111)
        self.FigSubPlot1.axes.get_xaxis().set_visible(False)
        self.FigSubPlot1.axes.get_yaxis().set_visible(False)
        psy = abs(self.psy0*self.psy0)
        self.im1 = self.FigSubPlot1.imshow(psy+self.Vim, extent=[-self.Nx/20,self.Nx/20,self.Ny/20,-self.Ny/20], cmap = ccmap)
        self.im1.set_data(psy+self.Vim)
        PosX = float(self.posX.get())
        PosY = float(self.posY.get())
        k0x  = float(self.momX.get())
        k0y  = float(self.momY.get())
        VarX  = float(self.varX.get())
        VarY = float(self.varY.get())
            
        arrow = matplotlib.text.Annotation('', fontsize=20, xy=(PosX, PosY),
                            xytext=(PosX+k0x,PosY+k0y),
                            arrowprops=dict(arrowstyle="<-",
                                            linewidth = 5.,
                                            alpha =0.8,
                                            color = 'orange')
                                            )
        
        shift = self.ninput%8
        
            
        if self.ninput%2==0:
            xc = 520
            yc = 61+int(shift/2)*121
        else:
            xc = 633
            yc = 61+int(shift/2)*121
            
        self.FigSubPlot1.add_artist(arrow)
        self.canvas1 = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(self.fig1, master=self)
        self.canvas1.get_tk_widget().place(y=yc, x=xc)
        self.canvas1.get_tk_widget().configure(background="#ecf2f9", highlightcolor="#ecf2f9", highlightbackground="#ecf2f9")
        self.canvas1._tkcanvas.place(y=yc, x=xc)
        self.update()
        self.canvas1.show()
        self.logbook.append([self.ninput, self.V, PosX, PosY, VarX, VarY, k0x, k0y, self.strPot])
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
        print(mkmovie)
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
            string = str(i[0]) + "\t" + str(i[2]) + "\t" + str(i[3]) + "\t" + str(i[4]) 
            string+= "\t" + str(i[5]) + "\t" + str(i[6]) + "\t" + str(i[7]) + "\t" + str(i[8]) + "\n"
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

