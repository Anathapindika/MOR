# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 13:33:39 2015

@author: TinoValentin
"""

import numpy as np
import matplotlib.pylab as plt
from movie import*
import os

def laplacian(Array2d,dx,dy):

    x  = np.gradient(Array2d,dx)[1]    
    xx = np.gradient(x,dx)[1]
 
    y  = np.gradient(Array2d,dy)[0]
    yy = np.gradient(y,dy)[0]

    return xx+yy
    
    
def gaussian2d(x, y, Coef, x0, y0, varx, vary, k0x, k0y):
    
    return (Coef * np.exp(-0.5 * (((x - x0) * 1. / varx) ** 2 
    + ((y - y0) * 1. / vary) ** 2) + 1j * x * k0x + 1j * y *k0y))
    

    
class schrodinger2d(object):        
    pi = np.pi    
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
        self.dkx    = 2 * pi / (self.Nx * self.dx)
        self.dky    = 2 * pi / (self.Ny * self.dy)
        
        kx                  = (np.arange(self.Nx) - self.Nx/2) * self.dkx
        ky                  = (np.arange(self.Ny) - self.Ny/2) * self.dky
        self.Kx, self.Ky    = np.meshgrid(kx,ky) 

        self.steps()         
        
        self.psi_x    = psi0
        
        self.real_psi = np.abs(self.psi_x*self.psi_x)
        
        
    def steps(self):
        self.evolve_x_half      = np.exp(-(1j * (self.dt / 2) * self.V / self.hbar))
        
        self.evolve_x           = self.evolve_x_half * self.evolve_x_half
        
        self.evolve_k           = np.exp(-(1j * self.hbar * np.fft.fftshift(self.Kx * self.Kx
                                + self.Ky * self.Ky) * self.dt / (2 * self.m)))
                                
    def x_to_k(self):
        self.psi_k   = np.fft.fft2(self.psi_x)
        
    def k_to_x(self):
        self.psi_x    = np.fft.ifft2(self.psi_k)
        
    def snapshot(self, Steps):
        
        self.psi_x    *= self.evolve_x_half
        
        for i in range(Steps):
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
        self.t      += self.dt*Steps
    
    
def SVD_Modes(A, Mpath, Nx, Ny, PlotModes = True, nModes = 100):
    

    At      = np.conjugate(A.T)   
    AAt     = dot(A,At)
    sig,U   = linalg.eigh(AAt)

    sorter  = sig.argsort()[::-1] 
    sig     = sqrt(sig)
    sig     = sig[sorter]
    U       = U[:,sorter]


    path = Mpath + "sig.jpg"
    i = sig_print(sig,path)
    
    Ut      = np.conjugate(U.T)
    Sig     = diag(sig)
   
    UtA     = dot(Ut,A)
    V       = np.zeros_like(UtA, dtype=np.complex)
    Modes = np.zeros((len(V[:,0]),Ny,Nx), dtype=np.complex)
    
    for i in range(len(V[:,0])):
        V[i,:] = UtA[i,:]/linalg.norm(UtA[i,:])
        Modes[i,:,:] = np.reshape(V[i,:],(Ny,Nx))
   

    if PlotModes:        
        Mpath = Mpath + "Modes/"
        if not os.path.exists(Mpath):
            os.mkdir(Mpath)
        
        for i in range(nModes):
            imtemp = np.real(Modes[i,:,:])
            print_frame(imtemp,Mpath+"mode_"+str(i)+".png")
            
    Atemps = np.dot(U,Sig)

    return Modes, Atemps
    
def potential(X, Y, boundary, pot):
    
    V = np.zeros_like(X)
    boundx = abs(X[0,0]*boundary)
    boundy = abs(Y[0,0]*boundary)
    V[X>boundx]  = pot
    V[X<-boundx] = pot
    V[Y>boundy]  = pot
    V[Y<-boundy] = pot
    
    return V