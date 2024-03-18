from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# constant values
q = 0.1
b= 1/q
wD =1/q
g = 1/q
x0=1
o0=0
ph0=0

#function that returns the DE as list.
def DDP(y, t,g,b,wD):
     x,o,ph= y
     dydt = [o , - b*o - np.sin(x) + g * np.cos(wD*t),wD]
     return dydt

# intial conditon nad time scale
y0 = [ x0, o0, ph0]
t = np.linspace(0,100,1000)
sol = odeint(DDP, y0, t, args=(g,b,wD))
x_dot,o_dot ,ph_dot= sol[:, 0],sol[:, 1],sol[:,2]

# this limits angle in range -pi to pi
p = np.arctan2(np.sin(o_dot), np.cos(o_dot))
p2= np.arctan2(np.sin(x_dot), np.cos(x_dot))

x=[]
y=[]

fig, ax = plt.subplots()

def animate(i):
     x.append(np.sin(x_dot[i]))
     y.append(-np.cos(x_dot[i]))
     ax.clear()
     ax.scatter(x, y ,color='b', s=0.5)
     plt.xlim([-1.3,1.3])
     plt.ylim([-1.3,1.3])
     if i>5:
        del x[0]
        del y[0]


ani = FuncAnimation(fig, animate, frames=len(t), interval=.05, repeat=False)

plt.show()
