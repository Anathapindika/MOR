import numpy as np
import numpy.random as npr
from scipy.interpolate import NearestNDInterpolator

N = 100
B = 3

g = np.ndarray(shape=(N,B))
v = np.ndarray(shape=(N))

def fun(x):
    return 100*x[0] + 10*x[1] + x[2]

for n in range(N):
    r = npr.sample(B)
    g[n,:] = r
    v[n] = fun(r)

print(g)
print(v)

tree = NearestNDInterpolator(g,v)

print(tree)

x = npr.sample(B)

print(x,fun(x),tree(x))



