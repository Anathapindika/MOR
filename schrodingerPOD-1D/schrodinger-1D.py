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
def laplacian(Array,dx):

    x  = np.gradient(Array,dx)    
    xx = np.gradient(x,dx)
    return xx
    
    
#Initial Gaussian - k0x k0y are the initial Momentum
def gaussian1d(x, y, Coef, x0, y0, varx, vary, k0x, k0y):
    
    return (Coef * np.exp(-0.5 * (((x - x0) * 1. / varx) ** 2 
    ) + 1j * x * k0x ))


#Singular Value decomposition - returns the TimeCoefficiant (Atemps) and more importatntly the POD-Modes
def SVD_Modes(A, Mpath, N, nModes = 10): 
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
    print(shape(UtA))
    print(N)
    print(shape(A))
    Modes = np.zeros((len(V[:,0]),N), dtype=np.complex)
    
    for i in range(len(V[:,0])):
        V[i,:] = UtA[i,:]/linalg.norm(UtA[i,:])
        Modes[i,:] = V[i,:]
   

    if nModes:        
        Mpath = Mpath + "Modes/"
        if not os.path.exists(Mpath):
            os.mkdir(Mpath)
        
        for i in range(nModes):
            imtemp = np.real(Modes[i])
            print_frame(imtemp,Mpath+"mode_"+str(i)+".png")
            
    Atemps = np.dot(U,Sig)
    Ap     = np.dot(Atemps,V)

    #Assert if SVD failled
#    assert(np.allclose(A,Ap))
    return Modes, Atemps


#Schrodinger Approxiation for initial Data
class schrodinger2d(object):  

    """
    initial Values are The Meshgrid X,Y, the initial waveFunction, the Potential
    and dt.
    """      
   
    def __init__(self, X, psi0, V, dt
                 , hbar=1, m=1, t0=0.0):
                  
        self.dt     = dt
        self.X      = X
        self.dx     = X[1]-X[0]
        self.V      = V
        self.hbar   = hbar
        self.m      = m
        self.t      = t0
        self.N      = len(X)
        self.dk    = 2 * np.pi / (self.N * self.dx)
        
        K           = (np.arange(self.N) - self.N/2) * self.dk
        self.K     = np.fft.fftshift(K)

        self.steps()         
        
        self.psi_x    = psi0
        
        self.real_psi = np.abs(self.psi_x*self.psi_x)
        
     
    def steps(self):
        self.evolve_x_half      = np.exp(-(1j * (self.dt / 2) * V / hbar))
        
        self.evolve_x           = self.evolve_x_half * self.evolve_x_half
        
        self.evolve_k           = np.exp(-(1j * self.hbar * (self.K * self.K) * self.dt / (2 * self.m)))
                                
                            
    def x_to_k(self):
        self.psi_k   = np.fft.fft(self.psi_x)
        
    def k_to_x(self):
        self.psi_x    = np.fft.ifft(self.psi_k)
    
        
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
    N          = pow(2,9)
    dx          = 0.1
    X           = (np.arange(N) - N/2) * dx

#    Initiating TimeValues    
    t0          = 0.0
    dt          = 0.005
 
#    Calculating the initail wavefunction psi0 
    x0          = 0
    y0          = 0
    varx        = 1
    vary        = 1
    k0x         = 1
    k0y         = 1
    Coef        = 1
    psi0        = gaussian1d(X,Y,Coef,x0,y0,varx,vary,k0x,k0y)
    
    hbar        = 1
    m           = 1
#   Initiating the Potential   
    boundary    = abs(0.7*X[0])
    pot         = 10E8  
    V           = X**2
   
 #  Creating Schrodinger Object - this should finnaly be in a loop with several Inital wavefunction
    s = schrodinger2d(X = X,
                        psi0 = psi0,
                        V = V,
                        dt = dt,
                        hbar = hbar,
                        m = 1,
                        t0 = t0)
     
# Performing the Calculation of the Wavefunction - A will be the DataMatrix - as above this had to bigger in the end
    snapshots = 10  
    betweenSnp = 50
    A           = np.zeros((snapshots,len(psi0)), dtype=np.complex)
    Movie1      = np.zeros((len(psi0),snapshots))
    
    print("Start with Leapfrog")     
    for i in range(snapshots):
        Movie1[:,i] = s.real_psi
        A[i,:] = s.psi_x
        s.snapshot(betweenSnp)
      
    ani_frame(X,Movie1,"test.mp4", "leapfrog")
        
    dt = dt*betweenSnp

# Calculating the Modes    
    Mpath = "Modes/"
    if not os.path.exists(Mpath):
        os.mkdir(Mpath)
    print("Start with SVD")
    Modes, Atemps = SVD_Modes(A,Mpath,N, nModes=10)
    
# For simplicity and calculating time just taking the first "cut" Modes
# NOTE: if cut>2 ode will NOT converge
    cut = 5
    
# Inserting the same initial Value for the galerkin Eq as above the calculated Wave with leapfrog!
    y0 = Atemps[0,:cut]
 
# Calculating the Laplacian of every POD-Mode   
    DModes = np.zeros((cut,N), np.complex)   
    Coef = np.zeros((cut,N), np.complex) 
    
    for i in range(cut):
        DModes[i,:] = laplacian(Modes[i],dx)
        Coef[i,:]   = (DModes[i]*hbar/2-V*Modes[i]/hbar)
    
    H = np.zeros((cut,cut), dtype=np.complex)
    for k in range(cut):
        for l in range(cut):
            H[k,l] = (np.conjugate(Modes[k])*Coef[l]).sum()
        
    assert(np.allclose(H,conjugate(H.T)))


    def fun(t,y):
        """
        returns:        da_k/dt =  Σ_i a_i H[k,i] (here with hbar=1)
        the summation of i is performed in the middle loop
        the integral is performed in the last line before returning alpha
        alpha in the return is a vektor of length cut
        """
        print(t)
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
    r = ode(fun).set_integrator('zvode', method='bdf', with_jacobian=False)
    r.set_initial_value(y0, t0)
    
    steps = 4
    Tf = dt*steps
    print("Start with Ode")
    while r.successful() and r.t < Tf:
        POD_A.append(r.y)
        r.integrate(r.t+dt)
        solution.append([r.t,r.y])
        
    for i in range(steps):
        print("Leapfrog: ", Atemps[i,0:cut])
        print("POD: ",POD_A[i])
        print("Allclose: ", allclose(Atemps[i,0:cut],POD_A[i]),"\n")
        