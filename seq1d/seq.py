# First we redefine fft and ifft to automate wrap/unwrap.

#import scipy.fftpack as scf
#def fft(fx):
#    return scf.fftshift(scf.fft(scf.fftshift(fx))) 
#def ifft(fk):
#    return scf.fftshift(scf.ifft(scf.fftshift(fk)))
#

from scipy.fftpack import fft, ifft

# Set up the grid in x and k.

import numpy as np
from numpy import arange, pi, exp, array, matrix, zeros

N = 2**12             # Assumed a power of two.
R = (pi*N/2)**.5     # Arbitrary, but this choice makes x and k ranges equal.
dx = 2*R/N
X = (arange(N)-N/2)*dx
dk = 2*pi/(N*dx)
K = (arange(N)-N/2)*dk
K[:N/2],K[N/2:] = 1*K[N/2:],1*K[:N/2]

# Choose potential and initial conditions.

#V = X**2
sep = X[-1]/4
psi = exp(-(X-sep)**2/2) + exp(-(X+sep)**2/2) + 0j
psi_k = fft(psi)

def V(psi):
    rho = psi*psi.conj()
    phi_k = fft(rho)
    phi = abs(phi_k)
    phi_k[1:] /= -K[1:]**2
    phi = ifft(phi_k).real
    return phi

# The semi-spectral method for the time-dependent Schroedinger equation.

def evolve():
    delta_t = 0.05
    global psi,psi_k
    psi *= exp(1j*V(psi)*delta_t/2)
    psi_k = fft(psi)
    psi_k *= exp(1j*K*K*delta_t)
    psi = ifft(psi_k)
    psi *= exp(1j*V(psi)*delta_t/2)
    psi_k = fft(psi)

S = N
T = 100

A = matrix(zeros(dtype=complex,shape=(T,S)))
for t in range(T):
    A[t,:] = 1*psi
    evolve()

M = A*A.conj().T
U = np.linalg.eigh(M)[1]
sigVT = U.conj().T*A
VT = 0*sigVT
for r in range(len(sigVT)):
    f = sigVT[r]*sigVT[r].conj().T
    norm = f[0,0]**(-0.5)
    VT[r] = norm * sigVT[r]

B = 16
Phi = matrix(zeros(dtype=complex,shape=(B,S)))
for b in range(B):
    Phi[b,:] = VT[-1-b,:]

Gal = A*Phi.conj().T    
Ab = Gal*Phi
    
#Phitld = 0*Phi
#    lapl = (Phi[b,:-2] + Phi[b,2:] - 2*Phi[b,1:-1])/(X[1]-X[0])**2
#    Phitld[b,1:-1] = lapl - V[1:-1]*Phi[b,1:-1]
#H = zeros(dtype=complex,shape=(B,B))
#for k in range(B):
#    for l in range(B):
#        H[k,l] = np.sum(Phi[k].conj()*Phitld[l])
# print(H)

# Now animate everything!

from matplotlib.pyplot import figure, show
from matplotlib.animation import FuncAnimation

def grinit():
    global fig,panelx,panelk,curvx_re,curvx_im,curvV
    fig.subplots_adjust(hspace=0.3)
    panelx = fig.add_subplot(211,axisbg='black')
    panelx.set_xlabel('position')
    panelx.set_xlim(X[0],X[-1])
    curvx_re = panelx.plot(X,psi.real,color='white',lw=3)[0]
    curvx_im = panelx.plot(X,psi.imag,color='magenta',lw=3)[0]
    panelV = fig.add_subplot(212,axisbg='black')
    panelV.set_xlabel('potential')
    panelV.set_xlim(X[0],X[-1])
    curvV = panelV.plot(X,V(psi),color='red')[0]

def frame(t):
    psi = array(A[t%T]-Ab[t%T])[0,:]
    maxx = max(abs(psi))
    panelx.set_ylim(-maxx,maxx)
    curvx_re.set_ydata(psi.real)
    curvx_im.set_ydata(psi.imag)
    curvV.set_ydata(V(psi))

fig = figure()
dummy = FuncAnimation(fig, frame, init_func=grinit, interval=25)

#fig.subplots_adjust(hspace=0.3)
#panelx = fig.add_subplot(211)
#panelx.set_xlabel('position')
#panelx.set_xlim(X[0],X[-1])
#for b in range(B):
#    Y = Phi[b]
#    panelx.plot(X,Y.real)
#    panelx.plot(X,Y.imag)

show()

