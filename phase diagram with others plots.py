from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

# constant values and initial values
q = 0.1
b= 1/q
wD = 1/q
g = 1/q
x0=1
o0=0
ph0=0

#function that returns the DDP as list.
def DDP(y, t,g,b,wD):
    x,o,ph= y
    dydt = [o, - b*o - np.sin(x) + g * np.cos(wD*t),wD]
    return dydt

# time scale and numerical solution 

y0 = [ x0, o0, ph0]
t = np.linspace(0,10,1000)
sol = odeint(DDP, y0, t, args=(g,b,wD))
x_dot,o_dot ,ph_dot= sol[:, 0],sol[:, 1],sol[:,2]

# this limits angle in range -pi to pi 

p1= np.arctan2(np.sin(x_dot), np.cos(x_dot))
p2= np.arctan2(np.sin(o_dot), np.cos(o_dot))
p3= np.arctan2(np.sin(ph_dot), np.cos(ph_dot)) 


#graphics

figure, axis = plt.subplots(2, 2)
# For phase plane
axis[0, 0].scatter(p1, p2,color='b',s=1)
axis[0, 0].set_title("phase plane")
# For theta / t
axis[0, 1].plot(t, p1)
axis[0, 1].set_title("theta(radian) / t(s)")
# For omega / t
axis[1, 0].plot(t,p2)
axis[1, 0].set_title("omega(radian/s) / t(s)")
# For phi /t
axis[1, 1].plot(t,p3)
axis[1, 1].set_title("phi(radian/s) / t(s)") 
# Combine all the operations and display
plt.show()