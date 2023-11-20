from random import randrange
import numpy as np
#Variables:
delta_t = 0.2
t = 0
t_max = 5
d_max = 4
angle_max = 90
k = 3.15
a = 9.2
w = 1.3
n = 3
c = 3.61
t_set = np.arange(0, t_max, delta_t)
r = np.zeros(2, n, len(t_set))
r[:][:][0] = randrange(1,3)
phi = np.zeros(n, len(t_set))
def d(i,p,t):
    return np.sqrt((r[1][p][t]-r[1][i][t])**2 + (r[2][p][t]-r[2][i][t])**2)
def angle(i,p,t):
    return phi[p][t]-phi[i][t]
def w(i,p,t):
    return ( a/(np.exp(w*d(i,p,t)+a)) )

for t in range(len(t_set)):
    for p in range(n):
        for i in range(n):
            if d(i,p,t)< d_max and np.absolute(angle(i,p,t))<angle_max:
                w_i = w(i,p,t)
                phi[p][t+1] = phi[p][t] - delta_t* (k/n) * w_i * np.sin(phi[i][t] - phi[p][t])