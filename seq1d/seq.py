# First we redefine fft and ifft to automate wrap/unwrap.

import scipy.fftpack as scf

def fft(fx):
    return scf.fftshift(scf.fft(scf.fftshift(fx)))

def ifft(fk):
    return scf.fftshift(scf.ifft(scf.fftshift(fk)))

# Set up the grid in x and k.

import numpy as np
from numpy import arange, pi, exp, array, matrix, zeros

N = 2**6             # Assumed a power of two.
R = (pi*N/2)**.5     # Arbitrary, but this choice makes x and k ranges equal.
step = 2*R/N
X = (arange(N)-N/2)*step
K = 2*pi/N*(arange(N)-N/2)/step

# Choose potential and initial conditions.

V = X**2
psi = exp(-X**2/2)+0j
psi_k = fft(psi)

# The semi-spectral method for the time-dependent Schroedinger equation.

def evolve():
    delta_t = 0.05
    global psi,psi_k
    psi *= exp(1j*V*delta_t/2)
    psi_k = fft(psi)
    psi_k *= exp(1j*K*K*delta_t)
    psi = ifft(psi_k)
    psi *= exp(1j*V*delta_t/2)
    psi_k = fft(psi)

S = N
T = 30

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


# Now animate everything!

from pylab import figure, plot, show
from matplotlib.animation import FuncAnimation

def grinit():
    global fig,panelx,panelk,curvx_re,curvx_im,curvk_re,curvk_im
    fig.subplots_adjust(hspace=0.3)
    panelx = fig.add_subplot(211,axisbg='black')
    panelx.set_xlabel('position')
    panelx.set_xlim(X[0],X[-1])
    curvx_re = panelx.plot(X,psi.real,color='white',lw=3)[0]
    curvx_im = panelx.plot(X,psi.imag,color='magenta',lw=3)[0]
#    panelk = fig.add_subplot(212,axisbg='black')
#    panelk.set_xlabel('momentum')
#    panelk.set_xlim(K[0],K[-1])
#    curvk_re = panelk.plot(K,psi_k.real,color='white',lw=3)[0]
#    curvk_im = panelk.plot(K,psi_k.imag,color='magenta',lw=3)[0]

def frame(t):
    print(t)
    print(A.shape)
    psi = array(A[t%T])[0,:]
    print(psi.shape)
    maxx = max(abs(psi))
#    maxk = max(abs(psi_k))
    panelx.set_ylim(-maxx,maxx)
#    panelk.set_ylim(-maxk,maxk)
    curvx_re.set_ydata(psi.real)
    curvx_im.set_ydata(psi.imag)
#    curvk_re.set_ydata(psi_k.real)
#    curvk_im.set_ydata(psi_k.imag)

fig = figure()
#dummy = FuncAnimation(fig, frame, init_func=grinit, interval=25)

fig.subplots_adjust(hspace=0.3)
panelx = fig.add_subplot(211)
panelx.set_xlabel('position')
panelx.set_xlim(X[0],X[-1])
for r in range(1,4):
    Y = array(VT[-r])[0]
    panelx.plot(X,Y.real)
    panelx.plot(X,Y.imag)

show()

