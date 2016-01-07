# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 10:25:09 2015

@author: TinoValentin
"""
import numpy as np
import matplotlib.pylab as plt 
from matplotlib import animation
from movie import*
from scipy.integrate import ode
import os
from scipy.special import hermite
from scipy.misc import factorial
from scipy.ndimage.interpolation import shift
import math

def eigenstates(nx,ny,wx,wy,X,Y,m=1,hbar=1):
    betax = np.sqrt(m*wx/hbar)
    resx  = np.sqrt(betax)*np.exp(-(m*wx*X**2)/(2*hbar))*hermite(nx)(betax*X)*1.0/np.sqrt(2**nx*factorial(nx))
    
    betay = np.sqrt(m*wy/hbar)
    resy  = np.sqrt(betay)*np.exp(-(m*wy*Y**2)/(2*hbar))*hermite(ny)(betay*Y)*1.0/np.sqrt(2**ny*factorial(ny))
    
    res   = resx*resy*1/np.sqrt(np.pi)+0j
    
    normalizer = abs(res*res).sum()
    res = res/np.sqrt(normalizer)
    
    return res


#Laplacian Module - was more exakt then the fouriertransform
def laplacian(Array2d,dx,dy):

    x  = np.gradient(Array2d,dx)[1]    
    xx = np.gradient(x,dx)[1]
 
    y  = np.gradient(Array2d,dy)[0]
    yy = np.gradient(y,dy)[0]

    return xx+yy
    
#Getting initialValues for ODE
def getA0(nModes,Modes,psy0):
    a0 = np.zeros(nModes, dtype=np.complex)
    
    for i in range(nModes):
        a0[i]  = np.sum(psy0*np.conjugate(Modes[i]))
    
    return a0
    
#Initial Gaussian - k0x k0y are the initial Momentum
def gaussian2d(x, y, x0, y0, varx, vary, k0x, k0y):
    Coef = 1
#    Coef = 1./(2*np.pi*varx*vary)
    psi = (Coef * np.exp(-0.5 * (((x - x0) * 1. / varx) ** 2 
        + ((y - y0) * 1. / vary) ** 2) + (1j * x * k0x + 1j * y *k0y)))
    normalizer = abs(psi*psi).sum()
    psi = psi/np.sqrt(normalizer)
    return psi

#Simple BoxPotential - boundary is the length of the boundarywall - pot the value of the potential 
def potential(X, Y, boundary, pot):
    
    V = np.zeros_like(X)
    boundx = abs(X[0,0]*boundary)
    boundy = abs(Y[0,0]*boundary)
    V[X>boundx]  = pot
    V[X<-boundx] = pot
    V[Y>boundy]  = pot
    V[Y<-boundy] = pot
    
    return V

#Singular Value decomposition - returns the TimeCoefficiant (Atemps) and more importatntly the POD-Modes
def SVD_Modes(A, Mpath, Nx, Ny, nModes = 100): 
    """
    Takes the Data Matrix A with shape (Nt,Nx*Ny) 
    The Path for printing Sigma Plot
    Prints the first nModes 
    """  
    
    At      = np.conjugate(A.T)   
    AAt     = np.dot(A,At)
    sig,U   = np.linalg.eigh(AAt)
    
    sorter  = sig.argsort()[::-1] 
    sig     = np.sqrt(sig)
    sig     = sig[sorter]
    U       = U[:,sorter]
    
        
    if math.isnan(np.amin(sig)):
        index = 0
        while not math.isnan(sig[index]):
            index +=1
        print("The SVD did converge too fast - had to cut the", len(sig), " entries to ", index, "entries")
        sig = sig[:index]
        U   = U[:,:index]


    path = Mpath + "sig.jpg"
    i = sig_print(sig,path)
    
    Ut      = np.conjugate(U.T)
    Sig     = np.diag(sig)
    
    cuter = len(sig)
    if nModes>cuter:
        nModes=cuter
        print("can only save ", cuter, " modes")
        
            
            
        print("Minimaler SigmaWert ",np.amin(sig), " Index", index)
        
    UtA     = np.dot(Ut,A)
    V       = np.zeros_like(UtA, dtype=np.complex)
    Modes = np.zeros((len(V[:,0]),Ny,Nx), dtype=np.complex)
    
    for i in range(len(V[:,0])):
        V[i,:] = UtA[i,:]/linalg.norm(UtA[i,:])
        Modes[i,:,:] = np.reshape(V[i,:],(Ny,Nx))
   

    if nModes:        
        Mpath = Mpath + "Modes/"
        if not os.path.exists(Mpath):
            os.mkdir(Mpath)
        
        for i in range(nModes):
            imtemp = np.real(Modes[i,:,:])
            print_frame(imtemp,Mpath+"mode_"+str(i)+".png")
            
    Atemps = np.dot(U,Sig)
    Ap     = np.dot(Atemps,V)

    #Assert if SVD failled
    assert(np.allclose(A,Ap))
    print("SVD succses:")
    print("True")
    return Modes, Atemps


#Schrodinger Approxiation for initial Data
class schrodinger2d(object):  

    """
    initial Values are The Meshgrid X,Y, the initial waveFunction, the Potential
    and dt.
    """      
   
    def __init__(self, X, Y, psi0, V, dt
                 , hbar=1, m=1, t0=0.0):
                  
        self.dt     = dt
        self.X      = X
        self.dx     = X[0,1]-X[0,0]
        self.Y      = Y
        self.dy     = Y[1,0]-Y[0,0]
        self.V      = V
        self.hbar   = hbar
        self.m      = m
        self.t      = t0
        self.Nx     = len(X)
        self.Ny     = len(Y)        
        self.dkx    = 2 * np.pi / (self.Nx * self.dx)
        self.dky    = 2 * np.pi / (self.Ny * self.dy)
        
        kx                  = (np.arange(self.Nx) - self.Nx/2) * self.dkx
        ky                  = (np.arange(self.Ny) - self.Ny/2) * self.dky
        self.Kx, self.Ky    = np.meshgrid(kx,ky) 
        self.Kx     = np.fft.fftshift(self.Kx)
        self.Ky     = np.fft.fftshift(self.Ky)

        self.steps()         
        
        self.psi_x    = psi0
        
        self.real_psi = np.abs(self.psi_x*self.psi_x)
        
     
    def steps(self):
        self.evolve_x_half      = np.exp(-(1j * (self.dt / 2) * V / hbar))
        
        self.evolve_x           = self.evolve_x_half * self.evolve_x_half
        
        self.evolve_k           = np.exp(-(1j * self.hbar * (self.Kx * self.Kx
                                + self.Ky * self.Ky) * self.dt / (2 * self.m)))
                                
                            
    def x_to_k(self):
        self.psi_k   = np.fft.fft2(self.psi_x)
        
    def k_to_x(self):
        self.psi_x    = np.fft.ifft2(self.psi_k)
    
        
    def snapshot(self, nSteps):
        """
        Main function where leapfrog is calculated 
        Takes "nSteps" - the value how many steps of timelength dt are evaluated
        befor giving back the self.psi_x
        """


        self.psi_x    *= self.evolve_x_half
        
        for i in range(nSteps-1):
            self.x_to_k()
            self.psi_k *= self.evolve_k
            self.k_to_x()
            self.psi_x *= self.evolve_x
        
        self.x_to_k()
        self.psi_k *= self.evolve_k
        self.k_to_x()
        
        self.psi_x    *= self.evolve_x_half
        
        self.real_psi = abs(self.psi_x*self.psi_x)
        self.t      += self.dt*nSteps
        
                                        
        
if __name__ == "__main__":
        
#    Initiating the Grid
    Nx          = pow(2,8)
    Ny          = pow(2,8)
    dx          = 0.08
    dy          = 0.08


#    Initiating TimeValues    
    t0          = 0.0
    dt          = 0.0001
    
    hbar        = 1
    m           = 1
#   Initiating the Potential  
    wx          = 1
    wy          = 2
    boundary    = 0.7
    pot         = 10E6  

    evaluation = 10
    snapshots = 30
    betweenSnp = 300
    A           = np.zeros((snapshots*evaluation,Nx*Ny), dtype=np.complex)
    Aindex      = 0
    x           = (np.arange(Nx) - Nx/2) * dx # + np.random.randint(-Nx/5,Nx/5)*dx
    y           = (np.arange(Ny) - Ny/2) * dy # + np.random.randint(-Ny/5,Ny/5)*dy
    X, Y        = np.meshgrid(x,y)
    V  = ((X*wx)**2+(Y*wy)**2)/2
        
    for evals in range(evaluation):
        print(evals)
        #    Calculating the initail wavefunction psi0- which are Eigenstates

#        psi0        = eigenstates(2*evals,2*evals,wx,wy,X,Y)

        x0          = np.random.rand()+np.random.randint(-5,4)
        y0          = np.random.rand()+np.random.randint(-5,4)
        varx        = 0.5+np.random.randint(3)+np.random.rand()
        vary        = 0.5+np.random.randint(3)+np.random.rand()
        k0x         = np.random.rand()+np.random.randint(-5,4)
        k0y         = np.random.rand()+np.random.randint(-5,4)
        psi0        = gaussian2d(X,Y,x0,y0,varx,vary,k0x,k0y)
        
     #  Creating Schrodinger Object - this should finnaly be in a loop with several Inital wavefunction
        s = schrodinger2d(X = X,
                                    Y = Y,
                                    psi0 = psi0,
                                    V = V,
                                    dt = dt,
                                    hbar = hbar,
                                    m = 1,
                                    t0 = t0)
         
    # Performing the Calculation of the Wavefunction - A will be the DataMatrix - as above this had to bigger in the end
        Movie1      = np.zeros((Ny,Nx,snapshots))
       
        print("Start with", str(evals), "th Input Leapfrog for Modes")     
        for i in range(snapshots):
            Movie1[:,:,i]   = s.real_psi
            A[Aindex,:]          = s.psi_x.flatten()
            Aindex              += 1
            s.snapshot(betweenSnp)
        normalizer = np.amax(Movie1)
        Movie1 *= 1./normalizer
        print("HIER ", np.amax(Movie1))
        ani_frame(Movie1,"leapfrog_" + str(evals) + "Input.mp4","leapfrog_input")
        

# Calculating the Modes    
    Mpath = "Modes/"
    if not os.path.exists(Mpath):
        os.mkdir(Mpath)
    print("Start with SVD")
    Modes, Atemps = SVD_Modes(A,Mpath,Nx,Ny, nModes=100)
    
    Epath = "Eigenstates/"
    if not os.path.exists(Epath):
        os.mkdir(Epath)
    #Printing Eigenstates
#    Modes = np.zeros((100,Ny,Nx), dtype=np.complex)
#    for i in range(10):
#        for k in range(1,100,10):
#            imtemp = eigenstates(i,k,wy,wy,X,Y)
#            Modes[i,:,:] = imtemp/10
#            print_frame(imtemp,Epath+"State_nx="+str(i)+"_ny="+str(k)+".jpg")
    
    dtpod = dt*betweenSnp 
    evaluation = 5

    for evals in range(evaluation):

        # Random initial Value (quoted for gaussian input)
        x0          = np.random.rand()+np.random.randint(-2,3)
        y0          = np.random.rand()+np.random.randint(-2,3)
        varx        = 1+np.random.rand()
        vary        = 1+np.random.rand()
        k0x         = np.random.rand()+np.random.randint(-2,3)
        k0y         = np.random.rand()+np.random.randint(-2,3)
        psi0        = gaussian2d(X,Y,x0,y0,varx,vary,k0x,k0y)

        
        longer = 1
        #Eigenstate input
#        psi0        = eigenstates(evals+np.random.randint(10),2*evals+np.random.randint(10),wx,wy,X,Y)
        s = schrodinger2d(X = X,
                            Y = Y,
                            psi0 = psi0,
                            V = V,
                            dt = dt,
                            hbar = hbar,
                            m = 1,
                            t0 = t0)
        Movie1      = np.zeros((Ny,Nx,snapshots*longer))
           
        print("Start with Leapfrog for evaluation")     
        for i in range(snapshots*longer):
            Movie1[:,:,i]   = s.real_psi
            s.snapshot(betweenSnp)

        
        # For simplicity and calculating time just taking the first "cut" Modes
        cut= 80
        if cut>len(Modes):
            cut=len(Modes)-1
            print("using Maximal Modes: ", cut)
        
        path = "evaluation"+str(evals)+"/"
        if not os.path.exists(path):
            os.mkdir(path)
        
        normalizer = np.amax(Movie1)
        Movie1 *= 1./normalizer
        ani_frame(Movie1,path + "leapfrog_evaluation.mp4","leapfrog_evaluation")
     
    # Calculating the Laplacian of every POD-Mode   
        DModes = np.zeros((cut,Ny,Nx), np.complex)   
        Coef = np.zeros((cut,Ny,Nx), np.complex)   
            
        
        for i in range(cut):
            DModes[i,:,:] = laplacian(Modes[i,:,:],dx,dy)
            Coef[i,:,:]   = (DModes[i,:,:]*hbar/(2*m)-V*Modes[i,:,:]/hbar)
        
        H = np.zeros((cut,cut), dtype=np.complex)
        
        for k in range(cut):
            for l in range(cut):
                H[k,l] = (np.conjugate(Modes[k])*Coef[l]).sum()
         
        
    # Galerkin Eq function for the ode   
        def fun(t,y):
            """
            returns:        da_k/dt =  Σ_i a_i H[k,i] (here with hbar=1)
            the summation of i is performed in the middle loop
            the integral is performed in the last line before returning alpha
            alpha in the return is a vektor of length cut
            """
    
            alpha = np.zeros(cut, dtype=np.complex) 
            
            for k in range(cut):                
                temp = 0+0j
                for i in range(cut):
                    temp   += y[i]*H[k,i]
                    
                alpha[k] = 1j*temp
            return alpha
        
        solution = []
        POD_A = []
    # Solving Ode    
        y0 = getA0(cut,Modes,psi0)
        r = ode(fun).set_integrator('zvode', method='bdf', with_jacobian=False)
        r.set_initial_value(y0, t0)
        
        steps = snapshots
        Tf = dtpod*steps*longer
        print("Start with Galerkin Ode")
        while r.successful() and r.t < Tf:
            POD_A.append(r.y)
            r.integrate(r.t+dtpod)
            solution.append([r.t,r.y])
        
    
        print("Done with Galerkin Ode")
        
    #    POD_Atemps = np.zeros((snapshots,cut), dtype=np.complex)
        
        ModesFlat = np.zeros((cut,Nx*Ny), dtype=np.complex)
        for i in range(cut):
            ModesFlat[i] = Modes[i].flatten()
        
    
        Ap = np.dot(POD_A,ModesFlat)
        Movie2 = np.zeros_like(Movie1)
        Movie3 = np.zeros_like(Movie1)
        for i in range(snapshots*longer):
            psi           = np.reshape(Ap[i,:],(Ny,Nx))
            Movie2[:,:,i] = abs(psi*psi)
            
        Movie2 *= 1./normalizer
        ani_frame(Movie2,path+"POD.mp4", "POD_with"+str(cut))
        difMovie = Movie2-Movie1
        ani_frame(difMovie,path+"Difference.mp4", "difference_with"+str(cut)+" Modes")  
        
#        for i in range(snapshots):
#            Movie3[:,:,i] = (Movie1[:,:,i] - Movie2[:,:,i])/Movie2[:,:,i]
#        ani_frame(Movie3, path+"Error.mp4", path+"error_with"+str(cut))
      
#    Movie2 = 
#    ani_frame(Movie1,"leapfrog.mp4","leapfrog")
###Plottting Values
##    POD_A = np.asarray(POD_A)
##    if cut==1:
##        tt = np.arange(steps)*dt 
##        plt.plot(tt,POD_A, 'rx', label="POD")
##        plt.plot(tt,Atemps[0:steps,0], label="Expected (from SVD)")
##        plt.legend(loc = "upper right")
##        plt.show()
##    elif cut==2:
##        tt = np.arange(steps)*dt 
##        plt.plot(tt,POD_A[:,0], 'rx', label="POD-Mode 1")
##        plt.plot(tt,Atemps[0:steps,0], label="Expected (from SVD) Mode 1")
##        plt.plot(tt,POD_A[:,1], 'gx', label="POD-Mode 2")
##        plt.plot(tt,Atemps[0:steps,1], label="Expected (from SVD) Mode 2")
##        plt.legend(loc = "upper right")
##        plt.show()
##    
