from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt


# constant values
q = 2
wD = 2/3
g = 1.5


#function that returns the DE as list.
def DDP(y, t,q,g,wD):
    o,x,ph= y
    dydt = [x, -x/q -np.sin(o) + g*np.cos(wD*t),wD]
    return dydt

# for intial conditon i.e angle and velocity of oscillation
y0 = [ 0, 0, 0]


# poincare space constants
p = (np.pi)/wD
n= 5000
T= n*p  #p=2pi/w_d
h=p/n   #step size
t = np.arange(0,T,h) 

sol = odeint(DDP, y0, t,args=(q,g,wD))
x1 = [sol[n*i, 0] for i in range(n)]
x = np.array(list(map(lambda x: x, x1))) 
z=  np.arctan2(np.sin(x), np.cos(x)) 
y = [sol[n*j, 1] for j in range(n)]  



#graphics
plt.scatter(z, y, color='b',s=0.5)
plt.xlabel(r'$\theta$ (radian)')
plt.ylabel(r'$\omega$ (radian/sec)')
plt.legend(loc='best')
plt.title('The Poincare section')
plt.grid(linestyle='--')
plt.show()