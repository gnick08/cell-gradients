import numpy as np
import math
import matplotlib.pyplot as plt


R=100
Q=4*R**3/3
rp=5
#x1,c=np.loadtxt("conc_D100_r100_v10_rp5.txt", unpack=True)
#x1,c=np.loadtxt("test_a600.txt", unpack=True)
#x2,c2=np.loadtxt("conc_s10_r500_c3.txt", unpack=True)
#c1=np.sum(c)
#print(c1)
dl=1

def f(x,eps):
    b=R-eps
    a=R*R/b
    q1=1
    q=-a*q1/R
    #return 2*(q1*(np.sqrt(R*R-x*x+(x+b)*(x+b))-np.abs(x+b))+q*(np.sqrt(R*R-x*x+(x+a)*(x+a))-np.abs(x+a)))/(R*R-x*x)
    return 2*(q1*(np.sqrt(R*R-x*x+(x+b)*(x+b))-np.abs(x+b))+q*(np.sqrt(R*R-x*x+(x+a)*(x+a))-np.abs(x+a)))
    #return 2*dl*(q1*(np.sqrt(R*R-x*x+(x+b)*(x+b))-np.abs(x+b))+q*(np.sqrt(R*R-x*x+(x+a)*(x+a))-np.abs(x+a)))/(math.pi*(R*R*dl-(z+0.5*dl)**3/3+(z-0.5*dl)**3/3))
    #return q1/np.abs(x+b) + q/np.abs(x+a)    


z=np.arange(-R+0.5*dl,R,dl)
l2=np.sum(f(z,rp))*dl
plt.yscale('log')
#print(l1)
#plt.plot(z,f(z,0.01)/(R*R-z*z)/l2,'-o')
plt.plot((z+R)/R,(Q*f(z,rp)/(R*R-z*z)/l2), label='analytic')
#plt.plot(x1-R,c/dl,'-o')
#plt.plot((x1)/R,(c/dl), '--')
np.savetxt('conc_D100_r100_v10_an.txt', np.c_[z+R, (Q*f(z,rp)/(R*R-z*z)/l2)])
#plt.plot((x1)/600,np.log(c))
#plt.plot(x2-R,c2)
#plt.plot(x2-R,np.log(c2))
plt.legend()
plt.show()
