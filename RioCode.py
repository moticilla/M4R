from random import random, choice, normalvariate
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import special
#Variables:
delta_t = 0.05
t_max = 50
n = 31
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
r        = np.zeros((len(t_set),n)) #this is r dot not r
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
def perturb(S, t, sign = 1, Heading = True, Speed = False):
    #S is the set of virtual neighbours to be perturbed either by speed or by heading at time t
    for i in S:
        if Heading == True:
            for all_t in range(int(t//delta_t),len(t_set)):
                phi[all_t, i] = phi[int(t//delta_t),i] + sign*(10/360)*2*pi
            for t_step in range(int(0.5//delta_t)):
                phi[int(t//delta_t)+t_step,i] = phi[int(t//delta_t),i] + sign*(10/360)*2*pi*ogive(t,0.083,t_step*delta_t-0.25)
        if Speed   == True:
            for all_t in range(int(t//delta_t),len(t_set)):
                r[all_t, i] = r[int(t//delta_t),i]   + sign*0.3
            for t_step in range(int(0.5//delta_t)):
                r[int(t//delta_t)+t_step,i]   = r[int(t//delta_t),i]   + sign*0.3*ogive(t,0.083,t_step*delta_t-0.25)

#set position of 30 neighbours for experiment 1 in two circles with real person in the middle
for i in range(14):
    position[0,i,0] = normalvariate(1.5,0.15)*np.sin(normalvariate(i*pi/7,(8/360)*2*pi))
    position[0,i,1] = normalvariate(1.5,0.15)*np.cos(normalvariate(i*pi/7,(8/360)*2*pi))
for i in range(14,30):
    position[0,i,0] = normalvariate(3.5,0.15)*np.sin(normalvariate((i-14)*pi/8,(8/360)*2*pi))
    position[0,i,1] = normalvariate(3.5,0.15)*np.cos(normalvariate((i-14)*pi/8,(8/360)*2*pi))
position[0,30,0] = 0
position[0,30,1] = 0
#plt.scatter(position[0,:,0], position[0,:,1])
#for i in range(31):
#    plt.annotate(i, (position[0,i,0], position[0,i,1]))
#plt.show()

#make all virtual people speed up to 1.3m/s in 3 s
for i in range(30):
    for t in range(int(3//delta_t)):
        r[t,i] = ogive(0,0.5,t*delta_t - 1.5)

#define subset S for experiment 1
#either 0,3,6,9,12 neighbours in S
S = np.array([0,1,28,5,7,12])
#make the perturbed set speed up after total of 5 seconds
perturb(S, 5, Heading = True, Speed = False)
print(phi[:,28])
#print(phi[:,1])
#get position of virtual neighbours from the perturbed speed and heading
for t in range(1, len(t_set)):
    for p in range(n):
            distance        = r * delta_t
            position[t,p,0] = position[t-1,p,0] + distance[t-1,p]*np.cos(phi[t-1,p])
            position[t,p,1] = position[t-1,p,1] + distance[t-1,p]*np.sin(phi[t-1,p])

#end trial after participant has walked 12m condition

#move the "real" person
for t in range(1, len(t_set)):
    p = 30
    r[t,p]   = r[t-1,p]
    phi[t,p] = phi[t-1,p]
    localN   = 0
    phi_sum  = 0
    r_sum    = 0
    for i in range(30):
        dis = d(i,p,t-1)
        ang = angle(i,p,t-1)
        if (dis< d_max) and (np.absolute(ang)<angle_max) and (i != p):#need to check but I think the angle is right now
            localN  += 1
            w_i      = w(i,p,t-1)
            phi_sum += w_i * np.sin(phi[t-1,i] - phi[t-1,p])
            r_sum   += w_i *(r[t-1,i]-r[t-1,p])
    if localN > 0:
        r[t,p]   = r[t-1,p]   + delta_t * (c/localN) * r_sum
        phi[t,p] = phi[t-1,p] - delta_t * (k/localN) * phi_sum
    distance        = r * delta_t
    position[t,p,0] = position[t-1,p,0] + distance[t-1,p]*np.cos(phi[t-1,p])
    position[t,p,1] = position[t-1,p,1] + distance[t-1,p]*np.sin(phi[t-1,p])

#set all initial positions randomly
'''
for i in range(n):
    position[0,i,0] = 0.5*random.random()
    position[0,i,1] = 0.5*random.random()
    phi[0,i] = 3* random.random()
    r[0,i]   = 0.3* random.random()
'''

#move according to model
'''
for t in range(1, len(t_set)):
    for p in range(n):
        r[t,p]   = r[t-1,p]
        phi[t,p] = phi[t-1,p]
        localN   = 0
        phi_sum  = 0
        r_sum    = 0
        for i in range(n):
            dis = d(i,p,t-1)
            ang = angle(i,p,t-1)
            if (dis< d_max) and (np.absolute(ang)<angle_max) and (i != p):#need to check but I think the angle is right now
                localN  += 1
                w_i      = w(i,p,t-1)
                phi_sum += w_i * np.sin(phi[t-1,i] - phi[t-1,p])
                r_sum   += w_i *(r[t-1,i]-r[t-1,p])
        if localN > 0:
            r[t,p]   = r[t-1,p]   + delta_t * (c/localN) * r_sum
            phi[t,p] = phi[t-1,p] - delta_t * (k/localN) * phi_sum
        distance        = r * delta_t
        position[t,p,0] = position[t-1,p,0] + distance[t-1,p]*np.cos(phi[t-1,p])
        position[t,p,1] = position[t-1,p,1] + distance[t-1,p]*np.sin(phi[t-1,p])
'''
for i in range(n):
    plt.scatter(position[:,i,0],position[:,i,1])
#my_quiver = np.zeros((len(t_set), n, 2))
#for t in range(len(t_set)-1):
#    for p in range(n):
#        my_quiver[t,p,0] = -position[t,p,0] + position[t+1,p,0]
#        my_quiver[t,p,1] = -position[t,p,1] + position[t+1,p,1]
#my_quiver = my_quiver/200
#for i in range(n):
#    plt.quiver(position[::4,i,0], position[::4,i,1],my_quiver[::4,i,0], my_quiver[::4,i,1], units = 'width')

plt.show()

#read in data from the paper
exp1_LateralDeviation_strip = pd.read_csv('/home/ellie/Documents/M4R/exp1/Exp1_LateralDeviation_strip.txt', sep="\s+", skip_blank_lines=True, 
                   header=None, skiprows = [0,1,12,15,16,17,18, 29,32,33,34,35,46])
exp1_LateralDeviation_strip.columns = ["Subject", "Near 0", "Far 0", "Near 3", "Near 6", "Near 9", "Near 12", "Far 3", "Far 6", "Far 9", "Far 12"]
exp1_LateralDeviation_strip.drop([10,13,14,25,28,29,40], axis = 0)
#print(exp1_LateralDeviation_strip)
t = 5
x = [t_step*delta_t for t_step in range(int(0.5//delta_t))]
y = [(10/360)*2*pi*ogive(0,0.083,t_step*delta_t-0.25) for t_step in range(int(0.5//delta_t))]

x1 = np.arange(-5,5,0.01)
y1 = ogive(3,0.083,x1)
#plt.plot(x1,y1)
#plt.show()