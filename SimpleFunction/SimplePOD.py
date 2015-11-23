# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 08:53:17 2015

@author: TinoValentin
"""

import numpy as np
import matplotlib.pylab as plt
from numpy import*
from matplotlib.pylab import*
import matplotlib.animation as animation
import os
from movie import*
from scipy.integrate import ode 


#Singular Value Decomposition - Output POD-Modes   
def SVD_Modes(A, Nx, Ny, nModes = 10):
    

    At      = np.conjugate(A.T)   
    AAt     = dot(A,At)
    sig,U   = linalg.eigh(AAt)

    sorter  = sig.argsort()[::-1] 
    sig     = sqrt(sig)
    sig     = sig[sorter]
    U       = U[:,sorter]


    path =  "sig.jpg"
    i = sig_print(sig,path)
    
    Ut      = np.conjugate(U.T)
    Sig     = diag(sig)
   
    UtA     = dot(Ut,A)
    V       = np.zeros_like(UtA, dtype=np.complex)
    Modes = np.zeros((len(V[:,0]),Ny,Nx), dtype=np.complex)
    
    for i in range(len(V[:,0])):
        
        if linalg.norm(UtA[i])!= 0:
            V[i,:] = UtA[i,:]/linalg.norm(UtA[i,:])
            Modes[i,:,:] = np.reshape(V[i,:],(Ny,Nx))
       
        else:
            
            V[i,:] = zeros_like(UtA[i,:])
            Modes[i,:,:] = zeros_like(Modes[i-1])

    if nModes:        
        Mpath =  "Modes/"
        if not os.path.exists(Mpath):
            os.mkdir(Mpath)
        
        for i in range(nModes):
            imtemp = np.real(Modes[i,:,:])
            print_frame(imtemp,Mpath+"mode_"+str(i)+".png")
            
    Atemps = np.dot(U,Sig)

    return Modes, Atemps
    

#Simple function 
psy = lambda t,x,y: exp(10*1j*t)*(x**2+y**2)

#Initiating the Grid and Time Values
t0 = 0
tf = 15
dt = 0.1
T = arange(t0,tf,dt)
x = arange(-5,5,0.1)
y = arange(-5,5,0.1)
X,Y = meshgrid(x,y)
Nx  = len(X)
Ny = len(Y)
Nt = len(T)

# DataMatrix
A = zeros((Nt,Nx*Ny), dtype=complex)
Movie = zeros((Ny,Nx,Nt))
frame = 0
for t in T:
    temp   = psy(t,X,Y)
    A[t,:] = temp.flatten()
    Movie[:,:,frame] = real(temp)
    frame += 1
    
    
ani_frame(Movie,"Original.mp4","Versuch")
ModesAll, Atemps = SVD_Modes(A,Nx,Ny,nModes =10)

#Take only the first [cut] Modes - because fast decay
cut = 1
ModesCut = ModesAll[:cut,:,:]


#Galerkin EQ Function
def fun(t,y):    
    alpha = zeros((cut), dtype=complex)
    for k in range(cut):
        temp = zeros((Ny,Nx), dtype=complex)
        
        for i in range(cut):
            temp[:,:] += y[i]*ModesCut[i]
            
        temp = temp*conjugate(ModesCut[k])
        alpha[k] = 10j*temp.sum()
    return alpha


# Solving the GalerkinEquation with the same InitailValue as above    
y0 = Atemps[0,:cut]
r = ode(fun).set_integrator('zvode', method='bdf', with_jacobian=False)
r.set_initial_value(y0, t0)

POD_A = []
time = []
while r.successful() and r.t < tf:
    time.append(r.t)
    POD_A.append(r.y)    
    r.integrate(r.t+dt)


#Plotting expected vs POD Values
tt = arange(len(POD_A))*dt 
plt.plot(tt,POD_A, 'rx', label="POD")
fyy = lambda t: y0*np.exp(10j*t) 
yy = (fyy(tt))
plt.plot(tt,yy, label="Expected (~exp(10it))")
plt.legend(loc = "upper right")
plt.show()
