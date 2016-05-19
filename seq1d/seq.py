
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

sep = X[-1]/4
# psi = exp(-(X-sep)**2/2) + exp(-(X+sep)**2/2) + 0j

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
T = 250
L = 4

A = matrix(zeros(dtype=complex,shape=(L*T,S)))
for l in range(L):
    fr = 0.997 + l*0.002
    print(fr)
    psi = exp(-(X-sep)**2/2) + fr*exp(-(X+sep)**2/2) + 0j
    for t in range(T):
        A[l*T+t,:] = 1*psi
        evolve()    

M = A*A.conj().T
U = np.linalg.eigh(M)[1]
sigVT = U.conj().T*A
VT = 0*sigVT
for r in range(len(sigVT)):
    f = sigVT[r]*sigVT[r].conj().T
    norm = f[0,0]**(-0.5)
    VT[r] = norm * sigVT[r]

B = 64
Phi = matrix(zeros(dtype=complex,shape=(B,S)))

for b in range(B):
    Phi[b,:] = VT[-1-b,:]

Gal = A*Phi.conj().T    
    
dGal = 0*Gal
for l in range(L):
    lo,hi = l*T,(l+1)*T
    dGal[lo,:] = Gal[lo+1,:] - Gal[lo,:]
    dGal[lo+1:hi-1,:] = (Gal[lo+2:hi,:]-Gal[lo:hi-2,:])/2
    dGal[hi-1,:] = Gal[hi-1,:] - Gal[hi-2,:]

# Interpolation

from scipy.interpolate import NearestNDInterpolator

def realise(gal):
    T,B = gal.shape
    rgal = zeros(shape=(T,2*B))
    for b in range(B):
        galb = array(gal[:,b])[:,0]
        rgal[:,2*b] = galb.real
        rgal[:,2*b+1] = galb.imag
    return rgal
    
rgal = realise(Gal)
tree = []
for b in range(B):
    dgb = array(dGal[:,b])[:,0]
    f = NearestNDInterpolator(rgal,dgb)
    tree.append(f)

    
def deriv(galt):
    f = matrix(zeros(dtype=complex,shape=(1,B)))
    rgalt = realise(galt)
    for b in range(B):
        f[0,b] = tree[b](rgalt)[0]
    return f
    
Aa = matrix(zeros(dtype=complex,shape=(T,S)))
sim = matrix(zeros(dtype=complex,shape=(T,B))) 

psi = exp(-(X-sep)**2/2) + exp(-(X+sep)**2/2) + 0j
phim = matrix(zeros(dtype=complex,shape=(1,S)))
phim[0,:] = psi
sim[0,:] = phim*Phi.conj().T

sim[1,:] = sim[0,:] + deriv(sim[0,:])
for t in range(1,T-1):
    sim[t+1,:] = sim[t-1,:] + 2*deriv(sim[t,:])

Ab = sim*Phi

for t in range(T):
    Aa[t,:] = 1*psi
    evolve()

# Now animate everything!

from matplotlib.pyplot import figure, show
from matplotlib.animation import FuncAnimation

def grinit():
    global fig,panelx,panelk,curvx_re,curvx_im,curv,curvV
    fig.subplots_adjust(hspace=0.3)
    panelx = fig.add_subplot(211,axisbg='black')
    panelx.set_xlabel('position')
    panelx.set_xlim(X[0],X[-1])
    curvx_re = panelx.plot(X,psi.real,color='cyan',lw=2)[0]
    curvx_im = panelx.plot(X,psi.imag,color='magenta',lw=2)[0]
    curv = panelx.plot(X,abs(psi)**2,color='white',lw=3)[0]
    panelV = fig.add_subplot(212,axisbg='black')
    panelV.set_xlabel('potential')
    panelV.set_xlim(X[0],X[-1])
    curvV = panelV.plot(X,V(psi),color='red')[0]

def frame(t):
    psi = array(Aa[t%T,:]-Ab[t%T,:])[0]
#    psi = array(A[t%T,:]-A[t%T+T,:])[0]
    maxx = max(abs(psi))
    panelx.set_ylim(-maxx,maxx)
    curvx_re.set_ydata(psi.real)
    curvx_im.set_ydata(psi.imag)
    curv.set_ydata(abs(psi)**2)
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

