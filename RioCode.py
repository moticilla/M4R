from random import random, choice
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import special
#Variables:
delta_t = 0.2
t_max = 60
n = 30
d_max = 4
pi = 3.14
angle_max = pi/2
k = 3.15
a = 9.2
w_var = 1.3
c = 3.61
t_set    = np.arange(0, t_max, delta_t)
position = np.zeros((len(t_set), n, 2))
phi      = np.zeros((len(t_set),n))
r        = np.zeros((len(t_set),n))
def d(i,p,t):
    return np.sqrt((position[t,p,0]-position[t,i,0])**2 + (position[t,p,1]-position[t,i,1])**2)
def angle(i,p,t):
    theta = np.arctan((position[t,i,1]-position[t,p,1])/ (position[t,i,0] - position[t,p,0]))
    return phi[t,p]-theta
def w(i,p,t):
    dis = d(i,p,t)
    return ( a/(np.exp(w_var * dis)+a) )
def ogive(mu, sigma,t):
    return 0.5 *(1+special.erf((t-mu)/sigma*np.sqrt(2)))
def perturb(S, t, Heading = True, Speed = False):
    #S is the set of virtual neighbours to eb perturbed either by speed or by direction
    sign = choice([1,-1])
    for i in range(int(0.5//delta_t)):
        if Heading == True:
            phi[t+i] = phi[t] + sign*(10/360)*2*pi*ogive(0,0.083,i*delta_t-0.25)
        if Speed == True:
            r[t+i] = r[t] + 0.3*sign*ogive(0,0.083,i*delta_t-0.25)

#set position of 12 neighbours for exp1
for i in range(14):
    position[0,i,0] = 1.5*np.sin(pi/4 + i*pi/7)
    position[0,i,1] = 1.5*np.cos(pi/4 + i*pi/7)
for i in range(14,30):    
    position[0,i,0] = 3.5*np.sin(pi/4 + (i-14)*pi/8)
    position[0,i,1] = 3.5*np.cos(pi/4 + (i-14)*pi/8)
#plt.scatter(position[0,:,0], position[0,:,1])
#plt.show()

#make the 18 neighbours not in the origional view jitter

#make all virtual people speed up to 1.3m/s in 3 s
for i in range(3):
    for t in range(int(3//delta_t)):
        r[t,i] = ogive(0,0.5,t*delta_t - 1.5)
'''
for i in range(n):
    position[0,i,0] = 0.5*random.random()
    position[0,i,1] = 0.5*random.random()
    phi[0,i] = 3* random.random()
    r[0,i]   = 0.3* random.random()
'''

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
            if (dis< d_max) and (np.absolute(ang)<angle_max) and (i != p):#need to check but I think the angle is right now
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

plt.show()

#read in data from the paper
exp1_LateralDeviation_strip = pd.read_csv('/home/ellie/Documents/M4R/exp1/Exp1_LateralDeviation_strip.txt', sep="\s+", skip_blank_lines=True, 
                   header=None, skiprows = [0,1,12,15,16,17,18, 29,32,33,34,35,46])
exp1_LateralDeviation_strip.columns = ["Subject", "Near 0", "Far 0", "Near 3", "Near 6", "Near 9", "Near 12", "Far 3", "Far 6", "Far 9", "Far 12"]
exp1_LateralDeviation_strip.drop([10,13,14,25,28,29,40], axis = 0)
#print(exp1_LateralDeviation_strip)
#multiply left turn by -1 and add to right turns

x = np.linspace(-1.5, 1.5)
plt.plot(x, ogive(0,0.5,x))

plt.show()