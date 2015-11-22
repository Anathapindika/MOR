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

#Laplacian Module - was more exakt then the fouriertransform
def laplacian(Array2d,dx,dy):

    x  = np.gradient(Array2d,dx)[1]    
    xx = np.gradient(x,dx)[1]
 
    y  = np.gradient(Array2d,dy)[0]
    yy = np.gradient(y,dy)[0]

    return xx+yy
    
    
#Initial Gaussian - k0x k0y are the initial Momentum
def gaussian2d(x, y, Coef, x0, y0, varx, vary, k0x, k0y):
    
    return (Coef * np.exp(-0.5 * (((x - x0) * 1. / varx) ** 2 
    + ((y - y0) * 1. / vary) ** 2) + 1j * x * k0x + 1j * y *k0y))

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


    path = Mpath + "sig.jpg"
    i = sig_print(sig,path)
    
    Ut      = np.conjugate(U.T)
    Sig     = np.diag(sig)
   
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

        self.steps()         
        
        self.psi_x    = psi0
        
        self.real_psi = np.abs(self.psi_x*self.psi_x)
        
     
    def steps(self):
        self.evolve_x_half      = np.exp(-(1j * (self.dt / 2) * V / hbar))
        
        self.evolve_x           = self.evolve_x_half * self.evolve_x_half
        
        self.evolve_k           = np.exp(-(1j * self.hbar * np.fft.fftshift(self.Kx * self.Kx
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
        
        self.evolve_x_half
        self.x_to_k()
        
        self.real_psi = abs(self.psi_x*self.psi_x)
        self.t      += self.dt*nSteps
        
                                        
        
if __name__ == "__main__":
        
#    Initiating the Grid
    Nx          = pow(2,9)
    Ny          = pow(2,9)
    dx          = 0.1
    dy          = 0.1
    x           = (np.arange(Nx) - Nx/2) * dx
    y           = (np.arange(Ny) - Ny/2) * dy
    X, Y        = np.meshgrid(x,y)

#    Initiating TimeValues    
    t0          = 0.0
    dt          = 0.01
 
#    Calculating the initail wavefunction psi0 
    x0          = 0
    y0          = 0
    varx        = 5
    vary        = 5
    k0x         = 1
    k0y         = 0
    Coef        = 2
    psi0        = gaussian2d(X,Y,Coef,x0,y0,varx,vary,k0x,k0y)
    
    hbar        = 1
    m           = 1
#   Initiating the Potential   
    boundary    = 0.7
    pot         = 10E8   
    V           = potential(X,Y,boundary,pot)
   
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
    snapshots = 50
    betweenSnp = 40
    A           = np.zeros((snapshots,len(psi0.flatten())))
   
    print("Start with Leapfrog")     
    for i in range(snapshots):
        A[i,:] = s.psi_x.flatten()
        s.snapshot(betweenSnp)
        
# Calculating the Modes    
    Mpath = "Modes/"
    if not os.path.exists(Mpath):
        os.mkdir(Mpath)
    print("Start with SVD")
    Modes, Atemps = SVD_Modes(A,Mpath,Nx,Ny, nModes=10)
    
# For simplicity and calculating time just taking the first "Cut" Modes
    Cut = 5
# Inserting the same initial Value for the galerkin Eq as above the calculated Wave with leapfrog!
    y0 = Atemps[0,:Cut]
 
# Calculating the Laplacian of every POD-Mode   
    DModes = np.zeros((Cut,Ny,Nx), np.complex)   
    for i in range(Cut):
        DModes[i,:,:] = laplacian(Modes[i,:,:],dx,dy)
        
# Galerkin Eq function for the ode   
    def fun(t,y):
        """
        returns:        da_k/dt = ∫ Σ_i a_i (Δφ_i / 2 - V φ_i) φ_k* dx (here with hbar=1)
        the summation of i is performed in the first loop
        the integral is performed in the last line before returning alpha
        alpha in the return is a vektor of length Cut
        """
        print(t)
        alpha = np.zeros(Cut, dtype=np.complex) 
        for k in range(Cut):                
            temp = np.zeros((Ny,Nx), dtype=np.complex)
            
            for i in range(Cut):
                temp   += y[i]*(DModes[i,:,:]*hbar/2-V*Modes[i,:,:]/hbar)  
                
            temp = 1j*temp*np.conjugate(Modes[k,:,:])        
            alpha[k] = temp.sum()
        return alpha
    
    solution = []
# Solving Ode    
    r = ode(fun).set_integrator('zvode', method='bdf', with_jacobian=False)
    r.set_initial_value(y0, t0)
    
    Tf = 0.1
    print("Start with Ode")
    while r.successful() and r.t < Tf:
        r.integrate(r.t+dt)
        solution.append([r.t,r.y])


