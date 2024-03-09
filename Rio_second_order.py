#making changes
from random import random, choice, normalvariate
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import special
delta_t = 0.05
t_max = 15
n = 31
d_max = 5
pi = np.pi
angle_max = pi/2
k = 3.15
a = 9.2
w_var = 1.3
c = 3.61
t_set    = np.arange(0, t_max, delta_t)
position = np.zeros((len(t_set), n, 2))
phi      = np.zeros((len(t_set),n))
r        = np.zeros((len(t_set),n)) #this is r dot not r
def d(pos,position,p,t):
    return np.sqrt((position[t,p,0]-pos[t,:,0])**2 + (position[t,p,1]-pos[t,:,1])**2)
def angle(pos,position,p,t):
    theta = np.arctan(np.absolute(pos[t,:,1]-position[t,p,1])/ (pos[t,:,0] - position[t,p,0]))
    return (np.absolute(phi[t,p]-theta)-pi/2)*np.sign(theta) + pi/2
def w(pos,p,t):
    dis = d(pos,position,p,t)
    return ( a/(np.exp(w_var * dis)+a) )
def ogive(mu, sigma,t):
    return 0.5 *(1+special.erf((t-mu)/sigma*np.sqrt(2)))

def move_person():
    total_distance = 0
    t = 1
    global phi
    global r
    global position
    u = np.zeros((len(t_set),n)) #w = phi dot
    p = 30
    while total_distance<12 and t<len(t_set):
        r[t,p]   = r[t-1,p]
        phi[t,p] = phi[t-1,p]
        localN   = 0
        phi_sum  = 0
        r_sum    = 0
        reduced_position = position[:,:30,:]
        dis = d(reduced_position,position,p,t-1)
        ang = angle(reduced_position,position,p,t-1)
        w_i = w(reduced_position,p,t-1)
        for i in range(30):
            if (dis[i]< d_max) and (np.absolute(ang[i])<angle_max) and (i != p):
                localN  += 1
                phi_sum += w_i[i] * np.sin((phi[t-1,i] - phi[t-1,p]))
                r_sum   += w_i[i] *(r[t-1,i]-r[t-1,p])
        if localN > 0:
            u[t,p]   = u[t-1,p]   + delta_t * (k/localN) * phi_sum 
            r[t,p]   = r[t-1,p]   + delta_t * (c/localN) * r_sum
            phi[t,p] = phi[t-1,p] + delta_t * u[t-1,p]
        distance        = r * delta_t
        total_distance= np.sum(distance[:,30])
        position[t,p,0] = position[t-1,p,0] + distance[t-1,p]*np.cos(phi[t-1,p]*2*pi/360)
        position[t,p,1] = position[t-1,p,1] + distance[t-1,p]*np.sin(phi[t-1,p]*2*pi/360)
        t += 1
