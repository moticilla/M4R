from random import randrange, random
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#Variables:
delta_t = 0.2
t_max = 6
n = 7
d_max = 4
angle_max = 1.57 #pi/2
k = 3.15
a = 9.2
w_var = 1.3
c = 3.61
t_set    = np.arange(0, t_max, delta_t)
position = np.zeros((len(t_set), n, 2))
phi      = np.zeros((len(t_set),n))
r        = np.zeros((len(t_set),n))
for i in range(n):
    position[0,i,0] = 0.5*random()
    position[0,i,1] = 0.5*random()
    phi[0,i] = 3* random()
    r[0,i]   = 0.3* random()
def d(i,p,t):
    return np.sqrt((position[t,p,0]-position[t,i,0])**2 + (position[t,p,1]-position[t,i,1])**2)
def angle(i,p,t):
    theta = np.arctan((position[t,i,1]-position[t,p,1])/ (position[t,i,0] - position[t,p,0]))
    return phi[t,p]-theta
def w(i,p,t):
    dis = d(i,p,t)
    return ( a/(np.exp(w_var * dis)+a) )

for t in range(1, len(t_set)):
    for p in range(n):
        r[t,p]   = r[t-1,p]
        phi[t,p] = phi[t-1,p]
        localN = 0
        phi_sum = 0
        r_sum   = 0
        for i in range(n):
            dis     = d(i,p,t-1)
            ang     = angle(i,p,t-1)
            if (dis< d_max) and (np.absolute(ang)<angle_max) and (i != p):#the angle condition is wrong
                localN  += 1
                w_i      = w(i,p,t-1)
                phi_sum += w_i * np.sin(phi[t-1,i] - phi[t-1,p])
                r_sum   += w_i *(r[t-1,i]-r[t-1,p])
        if localN > 0:
            r[t,p]   = r[t-1,p]   + delta_t * (c/localN) * r_sum
            phi[t,p] = phi[t-1,p] - delta_t * (k/localN) * phi_sum
        distance = r * delta_t
        position[t,p,0] = position[t-1,p,0] + distance[t-1,p]*np.cos(phi[t-1,p])
        position[t,p,1] = position[t-1,p,1] + distance[t-1,p]*np.sin(phi[t-1,p])

for i in range(n):
    plt.scatter(position[:,i,0],position[:,i,1])
my_quiver = np.zeros((len(t_set), n, 2))
for t in range(len(t_set)-1):
    for p in range(n):
        my_quiver[t,p,0] = -position[t,p,0] + position[t+1,p,0]
        my_quiver[t,p,1] = -position[t,p,1] + position[t+1,p,1]
my_quiver = my_quiver/200
for i in range(n):
    plt.quiver(position[::4,i,0], position[::4,i,1],my_quiver[::4,i,0], my_quiver[::4,i,1], units = 'width')

#plt.show()

#read in data from the paper
exp1_LateralDeviation_strip = pd.read_csv('/home/ellie/Documents/M4R/Exp1_LateralDeviation_strip.txt', sep="\s+", skip_blank_lines=True, 
                   header=None, skiprows = [0,1,12,15,16,17,18, 29,32,33,34,35,46])
exp1_LateralDeviation_strip.columns = ["Subject", "Near 0", "Far 0", "Near 3", "Near 6", "Near 9", "Near 12", "Far 3", "Far 6", "Far 9", "Far 12"]
print(exp1_LateralDeviation_strip)