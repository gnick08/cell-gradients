import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.special as sp

font = {'size'   : 20}

mpl.rc('font', **font)

R=100
Z=500
s=0.01
num=40
V=Z
zeros=np.array(sp.jn_zeros(0,num))

x1,c1=np.loadtxt("cyl_D100_z500_r100.txt",usecols=(0,1), unpack=True)
#print(zeros)

def coeff(x):
    return np.sinh(x*s/R)/(sp.j1(x)*x*x*np.sinh(x*Z/R))

A=coeff(zeros)
A=A/A[0]
#print(A)

X=np.arange(0.5,Z,1)

def f(x):
    sum=0
    for i in range(num):
        sum+=A[i]*np.sinh(zeros[i]*(Z-x)/R)
    return sum

Y=f(X)
Y=Y/np.sum(Y)

np.savetxt("cylnum_z500_r100.txt", np.c_[X,V*Y])

plt.xlabel(r"z ($\mu m$)")
plt.ylabel(r"$c(z)/c_0$")
plt.yscale("log")
plt.plot(x1/10,c1,label="simulation")
plt.plot(X/10,V*Y,label="analytic")
#plt.plot(x1/10,c1)
#plt.plot(X/10,V*Y,label=f"first {num} terms")
plt.legend()
plt.show()
